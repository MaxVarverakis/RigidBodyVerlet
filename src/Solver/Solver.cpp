#include "Solver.hpp"

Solver::Solver(const uint window_size)
    : m_substeps    { 8 }
    , m_window_size { window_size }
    , m_particles(m_numParticles)
{
    m_circles.reserve(m_numParticles);

    m_grid_spacing = 2 * std::sqrt(m_mass);
    // m_grid_spacing = window_size;

    // for square domain
    // determine whether or not to pad grid with extra layer of cells
    double ratio { window_size / m_grid_spacing };
    if (ratio - static_cast<uint>(ratio) > 0.0)
    {
        std::cout << "Padding enabled!" << '\n';
        m_max_cell_idx = static_cast<std::size_t>(window_size / m_grid_spacing);
    }
    else
    {
        std::cout << "No padding!" << '\n';
        m_max_cell_idx = static_cast<std::size_t>(window_size / m_grid_spacing) - 1;
    }

    m_max_cell_idx_int = { static_cast<int>(m_max_cell_idx) };

    m_grid = std::vector<std::vector<std::size_t>>((m_max_cell_idx+1) * (m_max_cell_idx+1));
    m_neighbor_indices = std::vector<std::vector<std::size_t>>(m_grid.size());

    // maximum 4 particles can fit in cell if grid spacing = diameter and all particles have same diameter
    for (std::vector<std::size_t>& cell : m_grid) { cell.reserve(4); }

    precomputeNeighborIndices();

    for (uint i = 0; i < 10; ++i) { m_cached_fudge_factors[i] = generateFudgeFactor(1e-6); }
    
    m_workers.reserve(m_num_threads);
    m_cells_per_thread = static_cast<std::size_t>(m_grid.size() / m_num_threads);
    m_particles_per_thread = static_cast<std::size_t>(m_particles.size / m_num_threads);
    
}

void Solver::precomputeNeighborIndices()
{
    for (int cell_i = 0; cell_i <= m_max_cell_idx_int; ++cell_i)
    {
        for (int cell_j = 0; cell_j <= m_max_cell_idx_int; ++cell_j)
        {
            for (int k = 0; k < 5; ++k)
            {
                const auto& [i, j] = m_neighbor_offsets[k];
                
                int neighbor_i { (cell_i + i) };
                int neighbor_j { (cell_j + j) };
                
                // handle boundaries
                if
                (
                    neighbor_i < 0                                ||
                    neighbor_j < 0                                ||
                    neighbor_i > static_cast<int>(m_max_cell_idx) ||
                    neighbor_j > static_cast<int>(m_max_cell_idx)
                    ) { continue; }
                
                m_neighbor_indices
                [
                    flattenGridIndices(static_cast<std::size_t>(cell_i), static_cast<std::size_t>(cell_j))
                ].emplace_back
                (
                    flattenGridIndices(static_cast<std::size_t>(neighbor_i), static_cast<std::size_t>(neighbor_j))
                );
            }
        }
    }
}

double Solver::generateFudgeFactor(const double& max_fudgeness)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dis(-max_fudgeness, max_fudgeness);
    return dis(gen);
}

double Solver::getFudgeFactor(const double& hash_key)
{
    return m_cached_fudge_factors[static_cast<uint>(hash_key) % 10];
}

Point2D Solver::applyFudgeFactor(const Point2D& n_perp)
{
    const double& fudge { getFudgeFactor(n_perp.x()) };

    return
    {
        n_perp.x() * fudge,
        n_perp.y() * fudge
    };
}

std::size_t Solver::particlePositionToCellIndex(std::size_t particle_ID)
{
    return flattenGridIndices
        (
            static_cast<std::size_t>(m_particles.position[particle_ID].x() / m_grid_spacing),
            static_cast<std::size_t>(m_particles.position[particle_ID].y() / m_grid_spacing)
        );
}

void Solver::binNewParticle(std::size_t particle_ID)
{
    std::size_t cell_idx { particlePositionToCellIndex(particle_ID) };
    m_particles.cell_idx[particle_ID] = cell_idx;
    m_particles.cell_idx_idx[particle_ID] = m_grid[cell_idx].size();
    m_grid[cell_idx].emplace_back(particle_ID);
}

void Solver::removeParticleFromCell(std::size_t particle_ID)
{
    /*
    To avoid reshifting data after removal of particle,
    swap particle to be removed with last particle in
    cell array before removal.
    */

    std::size_t cell_idx = m_particles.cell_idx[particle_ID];
    std::size_t cell_idx_idx = m_particles.cell_idx_idx[particle_ID];
    
    if (cell_idx_idx != m_grid[cell_idx].size() - 1)
    {
        std::swap(m_grid[cell_idx][cell_idx_idx], m_grid[cell_idx].back());
        // Update the index of the particle that was swapped
        m_particles.cell_idx_idx[m_grid[cell_idx][cell_idx_idx]] = cell_idx_idx;
    }
    // Remove the last element (which was the particle that just moved)
    m_grid[cell_idx].pop_back();
    
    // m_grid[particle.cell_idx].erase(
    //     // remove moves all instances of particle.ID to the end of the vector and returns an iterator starting at the first move
    //     std::remove(m_grid[particle.cell_idx].begin(), m_grid[particle.cell_idx].end(), particle.ID),
    //     m_grid[particle.cell_idx].end()
    // );
}

void Solver::binParticles()
{
    for (std::size_t& particle_ID : m_particles.ID)
    {
        std::size_t cell_idx { particlePositionToCellIndex(particle_ID) };

        // only rebin particles that moved to a new cell
        if (cell_idx != m_particles.cell_idx[particle_ID])
        {
            // m_grid[cell_idx].erase(m_grid[cell_idx].begin() + particle.cell_idx_idx);
            removeParticleFromCell(particle_ID);

            m_particles.cell_idx[particle_ID] = cell_idx;
            m_particles.cell_idx_idx[particle_ID] = m_grid[cell_idx].size();
            m_grid[cell_idx].emplace_back(particle_ID);
        }
    }
}

std::array<std::size_t, 2> Solver::unflattenGridIndices(const std::size_t& idx)
{
    return {static_cast<std::size_t>(idx / (m_max_cell_idx+1)), idx % (m_max_cell_idx+1)};
}

std::size_t Solver::flattenGridIndices(const std::size_t& x_idx, const std::size_t& y_idx)
{
    return x_idx * (m_max_cell_idx + 1) + y_idx;
}

void Solver::resolveParticleCollisions(std::size_t particle_ID, std::size_t& other_particle_ID)
{
    double min_distance { m_particles.radius[particle_ID] + m_particles.radius[other_particle_ID] };
    if ( (m_particles.position[particle_ID] - m_particles.position[other_particle_ID]).mag2() < min_distance * min_distance )
    {
        Point2D r { m_particles.position[particle_ID] - m_particles.position[other_particle_ID] };
        Point2D n { r.normalized() };
        
        // adjust positions so there is no overlap
        double overlap { m_particles.radius[particle_ID] + m_particles.radius[other_particle_ID] - r.magnitude() };

        // introduce fudge factor so particles can't stack themselves out of the domain
        Point2D fudge_vec { applyFudgeFactor(n.perp()) };

        // particle.position       += overlap/2 * n + fudge_vec;
        // other_particle.position -= overlap/2 * n - fudge_vec;
        /*
        particle.position       += overlap/2 * n;
        other_particle.position -= overlap/2 * n;
        */

        double mass_ratio_i { 2 * m_particles.mass[other_particle_ID] / (m_particles.mass[particle_ID] + m_particles.mass[other_particle_ID]) };
        double mass_ratio_j { 2 *       m_particles.mass[particle_ID] / (m_particles.mass[particle_ID] + m_particles.mass[other_particle_ID]) };
        
        
        m_particles.position[particle_ID]       += overlap/2 * n * mass_ratio_i * m_damping + fudge_vec;
        m_particles.position[other_particle_ID] -= overlap/2 * n * mass_ratio_j * m_damping - fudge_vec;
        
        /*
        // for elastic collisions
        Point2D v { particle.velocity - other_particle.velocity };
        double v_normal { v.dot(n) / r.magnitude() }; // extra division by |r| for later formula
        
        // update velocities based off impulse
        Point2D v_i { mass_ratio_i * v_normal * (particle.radius + other_particle.radius) * n * m_damping };
        Point2D v_j { mass_ratio_j * v_normal * (particle.radius + other_particle.radius) * n * m_damping };
        
        particle.velocity       -= v_i;
        other_particle.velocity += v_j;
        */
    }
}

void Solver::resolveCollisions(std::size_t thread_ID)
{
    for (
            std::size_t cell_idx = (thread_ID * m_cells_per_thread);
            cell_idx < ( (thread_ID == m_num_threads - 1) ? m_grid.size() : (thread_ID + 1) * m_cells_per_thread );
            ++cell_idx
        )
    {
        for (std::size_t neighbor_idx : m_neighbor_indices[cell_idx])
        {
            if (!m_grid[cell_idx].empty() && !m_grid[neighbor_idx].empty())
            {
                // Recall: each cell contains particle ID's which correspond to their indices in the m_particles vector
                for (std::size_t i = 0; i < m_grid[cell_idx].size(); ++i)
                {
                    std::size_t particle_ID { m_particles.ID[m_grid[cell_idx][i]] };
                    for (std::size_t j = 0; j < m_grid[neighbor_idx].size(); ++j)
                    {
                        std::size_t other_particle_ID { m_particles.ID[m_grid[neighbor_idx][j]] };
                        if (particle_ID != other_particle_ID)
                        {
                            for (int i = 0; i < max_collision_resolution_count; ++i) { resolveParticleCollisions(particle_ID, other_particle_ID); }
                        }
                    }
                }
            }
        }
    }
}

void Solver::applyBoundaryConditions(std::size_t particle_ID)
{
    m_particles.position[particle_ID].setX(std::clamp(m_particles.position[particle_ID].x(), m_particles.radius[particle_ID], m_window_size - m_particles.radius[particle_ID]));
    m_particles.position[particle_ID].setY(std::clamp(m_particles.position[particle_ID].y(), m_particles.radius[particle_ID], m_window_size - m_particles.radius[particle_ID]));
}

void Solver::applyNoVelocityBC(std::size_t particle_ID)
{
    Point2D& position { m_particles.position[particle_ID] };
    Point2D& position_previous { m_particles.position_previous[particle_ID] };
    const double radius { m_particles.radius[particle_ID] };

    if (position.x() < radius)
    {
        position.setX(radius);
        position_previous.setX(radius - m_damping * (m_particles.position_previous[particle_ID].x() - radius));
    }
    else if (position.x() > m_window_size - radius)
    {
        position.setX(m_window_size - radius);
        position_previous.setX((m_window_size - radius) - m_damping * (m_particles.position_previous[particle_ID].x() - (m_window_size - radius)));
    }
    
    if (position.y() < radius)
    {
        position.setY(radius);
        position_previous.setY(radius - m_damping * (m_particles.position_previous[particle_ID].y() - radius));
    }
    else if (position.y() > m_window_size - radius)
    {
        position.setY(m_window_size - radius);
        position_previous.setY((m_window_size - radius) - m_damping * (m_particles.position_previous[particle_ID].y() - (m_window_size - radius)));
    }
}

void Solver::parallelizeCollisions()
{
    for (std::size_t thread_ID = 0; thread_ID < m_num_threads; ++thread_ID)
    {
        m_workers.emplace_back(&Solver::resolveCollisions, this, thread_ID);
    }

    for (auto& thread : m_workers)
    {
        thread.join();  // Join all threads to ensure they complete before moving on
    }

    m_workers.clear();
}

void Solver::updateParticlesPerThread(std::size_t thread_ID, const int substep, const double dt, Graphics& graphics, sf::RenderTarget& window_target)
{
    for (
            std::size_t particle_ID = (thread_ID * m_particles_per_thread);
            particle_ID < ( (thread_ID == m_num_threads - 1) ? m_particles.size : (thread_ID + 1) * m_particles_per_thread );
            ++particle_ID
        )
        {   
            Point2D position { m_particles.position[particle_ID] };

            // Point2D a_tot { m_a_global + m_config.mouse_force / particle.mass * graphics.m_mouse.particleInteraction(particle.position) };
            // here we clamp displacement so that "velocity" gets clamped: |dx/dt| < v_max ==> |dx| < v_max * dt = v_max / frame_rate
            Point2D a_tot { m_a_global + m_mouse_force * graphics.m_mouse.particleInteraction(position) };
            // Point2D a_tot { m_a_global };

            
            // Point2D next_r { 2 * particle.position - particle.velocity + a_tot*dt*dt };
            Point2D next_r { position + (position - m_particles.position_previous[particle_ID]).clamped(1.5) + a_tot*dt*dt };

            m_particles.position_previous[particle_ID] = position;
            m_particles.position[particle_ID] = next_r;
            applyNoVelocityBC(particle_ID);
            
            /*
            // update the position according to Verlet integration
            particle.position = particle.position + particle.velocity*dt + 0.5*a_tot*dt*dt;
            applyBoundaryConditions(particle);
            
            // Point2D next_r { particle.position + particle.velocity*dt + 0.5*a_tot*dt*dt };
            // set next position according to boundary conditions
            // particle.position.setX(applyBoundaryConditions(next_r.x(), particle.radius));
            // particle.position.setY(applyBoundaryConditions(next_r.y(), particle.radius));
            // trying to apply z-BC's in 2D simulation will cause damage!

            // update velocity
            particle.velocity += a_tot * dt;

            // clamp velocity so things don't explode... 200 per component seems reasonable but adjust accordingly
            particle.velocity.clamp(500);
            

            
            // From the position-setting block, now any particles
            // that were headed out of bounds will have their
            // appropriate position component equal to the domain bound
            

            // adjust velocity if particle hits boundary
            if (( particle.position.x() == m_window_size - particle.radius ) || ( particle.position.x() == particle.radius )) {
                particle.velocity.setX(-particle.velocity.x() * m_damping); // Reverse x-velocity if boundary hit
            }
            if (( particle.position.y() == m_window_size - particle.radius ) || ( particle.position.y() == particle.radius )) {
                particle.velocity.setY(-particle.velocity.y() * m_damping); // Reverse y-velocity if boundary hit
            }
            */

            std::size_t cell_idx { particlePositionToCellIndex(particle_ID) };
            // only rebin the particle if it moves to a new cell
            if (cell_idx != m_particles.cell_idx[particle_ID])
            {
                // m_grid[cell_idx].erase(m_grid[cell_idx].begin() + static_cast<std::vector<std::size_t>::difference_type>(particle.cell_idx_idx));

                removeParticleFromCell(particle_ID);

                m_particles.cell_idx[particle_ID] = cell_idx;
                m_particles.cell_idx_idx[particle_ID] = m_grid[cell_idx].size();
                m_grid[cell_idx].emplace_back(particle_ID);
            }

            if (substep + 1 == m_substeps)
            {
                // since we're already looping through particles, might as well render it now too!
                m_circles[particle_ID].setPosition(m_particles.position[particle_ID].toVf());
                renderParticle(particle_ID, window_target);

            }
        }
}

void Solver::parallelizeUpdates(const int substep, const double dt, Graphics& graphics, sf::RenderTarget& window_target)
{
    for (std::size_t thread_ID = 0; thread_ID < m_num_threads; ++thread_ID)
        {
            m_workers.emplace_back([this, thread_ID, substep, dt, &graphics, &window_target]
            {
                this->updateParticlesPerThread(thread_ID, substep, dt, graphics, window_target);
            });
        }

        for (auto& thread : m_workers)
        {
            thread.join();  // Join all threads to ensure they complete before moving on
        }
        
        m_workers.clear();
}

void Solver::update(Graphics& graphics, sf::RenderTarget& window_target)
{
    const double& dt { m_dt / m_substeps };
    
    for (int substep = 0; substep < m_substeps; ++substep)
    {
        parallelizeCollisions();
        
        parallelizeUpdates(substep, dt, graphics, window_target);

    }
}

sf::Text Solver::cellNumberText(const std::size_t& idx, Graphics& graphics)
{
    float character_size { 16 };
    sf::Text text(graphics.m_font, std::to_string(idx));
    text.setCharacterSize(static_cast<unsigned int>(character_size)); // in pixels, not points!
    text.setFillColor(sf::Color::Green); // set the color
    
    const auto& [idx_x, idx_y] = unflattenGridIndices(idx);
    text.setPosition(
    {
        (static_cast<float>(idx_x) + 0.5f) * static_cast<float>(m_grid_spacing) - character_size / 2,
        (static_cast<float>(idx_y) + 0.5f) * static_cast<float>(m_grid_spacing) - character_size / 2
    });

    return text;
}

void Solver::drawCells(sf::RenderTarget& target_window)
{
    for (std::size_t i = 0; i < m_grid.size(); ++i)
    {
        sf::RectangleShape boundingBox(
        {
            static_cast<float>(m_grid_spacing),
            static_cast<float>(m_grid_spacing)
        });
        
        const auto& [idx_x, idx_y] = unflattenGridIndices(i);

        boundingBox.setPosition(
        {
            (static_cast<float>(idx_x)) * static_cast<float>(m_grid_spacing),
            (static_cast<float>(idx_y)) * static_cast<float>(m_grid_spacing)
        });

        boundingBox.setOutlineColor(sf::Color::Green);
        boundingBox.setOutlineThickness(0.5f);
        boundingBox.setFillColor(sf::Color::Transparent);
        target_window.draw(boundingBox);
    }
}

void Solver::drawCells(sf::RenderTarget& target_window, Graphics& graphics)
{
    for (std::size_t i = 0; i < m_grid.size(); ++i)
    {
        target_window.draw(cellNumberText(i, graphics));

        sf::RectangleShape boundingBox(
        {
            static_cast<float>(m_grid_spacing),
            static_cast<float>(m_grid_spacing)
        });
        
        const auto& [idx_x, idx_y] = unflattenGridIndices(i);

        boundingBox.setPosition(
        {
            (static_cast<float>(idx_x)) * static_cast<float>(m_grid_spacing),
            (static_cast<float>(idx_y)) * static_cast<float>(m_grid_spacing)
        });

        boundingBox.setOutlineColor(sf::Color::Green);
        boundingBox.setOutlineThickness(0.5f);
        boundingBox.setFillColor(sf::Color::Transparent);
        target_window.draw(boundingBox);
    }
}

void Solver::particleToCircleShape(std::size_t particle_ID)
{
    const float radius { static_cast<float>(m_particles.radius[particle_ID]) };
    
    sf::CircleShape circle(radius);
    circle.setOrigin({radius, radius});
    circle.setPosition(m_particles.position[particle_ID].toVf());
    circle.setFillColor(m_particles.color[particle_ID]);

    m_circles.emplace_back(circle);
}

void Solver::renderParticle(std::size_t particle_ID, sf::RenderTarget& window_target)
{
    window_target.draw(m_circles[particle_ID]);
}

void Solver::renderParticles(sf::RenderTarget& window_target)
{
    for (std::size_t particle_ID : m_particles.ID)
    {
        renderParticle(particle_ID, window_target);
    }
}

void Solver::arrowKeyGravity()
{
    double a_mag { m_a_global.magnitude() };

    if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right))
    {
        m_a_global = {a_mag, 0.0};
    }
    else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left))
    {
        m_a_global = {-a_mag, 0.0};
    }
    else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Up))
    {
        m_a_global = {0.0, -a_mag};
    }
    else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Down))
    {
        m_a_global = {0.0, a_mag};
    }
}

void Solver::drawMouse(sf::RenderWindow& window,Graphics& graphics)
{
    if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
    {
        graphics.m_mouse.setCursorPosition(sf::Mouse::getPosition(window));
        
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::LShift))
        {
            graphics.m_mouse.setMouseRadiusColor(sf::Color::Green);
        }
        else
        {
            graphics.m_mouse.setMouseRadiusColor(sf::Color::Red);
        }
        
        window.draw(graphics.m_mouse.m_shape);
    }
}

void Solver::run()
{
    Graphics graphics(50.0f);
    
    ParticleGenerator generator
    (
        1.0f, // milliseconds
        m_numParticles, // particles
        m_mass,
        {static_cast<double>(m_window_size)/4.0, static_cast<double>(m_window_size)/8.0},
        {1.0, 1.0}
    );

    float elapsed_time { 0.0f };
    float fps_elapsed_time { 0.0f };
    float generator_elapsed_time { 0.0f };
    sf::Text fps_text { graphics.fpsText() };
    
    // sf::ContextSettings settings;
    // settings.antiAliasingLevel = 1;
    
    // sf::RenderWindow window(sf::VideoMode({m_window_size, m_window_size}), "Simulation :)", sf::Style::Default, sf::State::Windowed, settings);
    sf::RenderWindow window(sf::VideoMode({m_window_size, m_window_size}), "Simulation :)");

    window.setFramerateLimit(m_frame_rate);

    while (window.isOpen())
    {
        // check all the window's events that were triggered since the last iteration of the loop
        while (const std::optional event = window.pollEvent())
        {
            // "close requested" event: we close the window
            if (event->is<sf::Event::Closed>() || sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Escape)) { window.close(); }
        }
        
        if (generator.generateBool(generator_elapsed_time) && m_particles.size < generator.m_numParticles)
        {
            generator_elapsed_time = 0.0f;

            for (int counter = 0; counter < m_numGenerators; ++counter)
            {
                const std::size_t& ID { m_particles.size };
                generator.generate(m_particles, m_dt, {counter * 100.0, 0.0});
                particleToCircleShape(ID);
                renderParticle(ID, window);
                binNewParticle(ID);

                m_particles.size += 1;
            }
        }

        arrowKeyGravity();
        
        // clear the window with black color
        window.clear(sf::Color::Black);

        // must keep update after the window clear since update also renders the particles
        update(graphics, window);

        // draw cell borders (for debugging purposes). Pass `graphics` if also want cell indices displayed
        // drawCells(window);

        // renderParticles(window);

        if (fps_elapsed_time  >= 500.0f)
        {
            fps_elapsed_time = 0.0f;
            fps_text = graphics.fpsText();
        }
        elapsed_time = static_cast<float>(graphics.m_clock.restart().asMilliseconds());
        generator_elapsed_time += elapsed_time;
        fps_elapsed_time       += elapsed_time;
        
        window.draw(fps_text);
        window.draw(graphics.particleCount(m_particles.size));
        drawMouse(window, graphics);
        
        window.display();
    }
}
