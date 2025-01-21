#include "Solver.hpp"

Solver::Solver(const uint window_size)
    : m_substeps    { 8 }
    , m_window_size { window_size }
{
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
                
                std::size_t neighbor_i { static_cast<std::size_t>(cell_i + i) };
                std::size_t neighbor_j { static_cast<std::size_t>(cell_j + j) };
                
                // handle boundaries
                if
                (
                    neighbor_i < 0              ||
                    neighbor_j < 0              ||
                    neighbor_i > m_max_cell_idx ||
                    neighbor_j > m_max_cell_idx
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

Point2D Solver::applyFudgeFactor(const Point2D& n_perp, const double& max_fudgeness)
{
    double fudge { generateFudgeFactor(max_fudgeness) };

    return
    {
        n_perp.x() * fudge,
        n_perp.y() * fudge
    };
}

std::size_t Solver::particlePositionToCellIndex(Particle& particle)
{
    return flattenGridIndices
        (
            static_cast<std::size_t>(particle.position.x()/m_grid_spacing),
            static_cast<std::size_t>(particle.position.y()/m_grid_spacing)
        );
}

void Solver::binNewParticle(Particle& particle)
{
    std::size_t cell_idx { particlePositionToCellIndex(particle) };
    particle.cell_idx = cell_idx;
    particle.cell_idx_idx = m_grid[cell_idx].size();
    m_grid[cell_idx].emplace_back(particle.ID);
}

void Solver::removeParticleFromCell(Particle& particle)
{
    /*
    To avoid reshifting data after removal of particle,
    swap particle to be removed with last particle in
    cell array before removal.
    */

    std::size_t cell_idx = particle.cell_idx;
    std::size_t cell_idx_idx = particle.cell_idx_idx;
    
    if (cell_idx_idx != m_grid[cell_idx].size() - 1)
    {
        std::swap(m_grid[cell_idx][cell_idx_idx], m_grid[cell_idx].back());
        // Update the index of the particle that was swapped
        m_particles[m_grid[cell_idx][cell_idx_idx]].cell_idx_idx = cell_idx_idx;
    }
    // Remove the last element (which was the particle that just moved)
    m_grid[cell_idx].pop_back();
    
    // m_grid[particle.cell_idx].erase(
    //     // remove moves all instances of particle.ID to the end of the vector and returns an iterator starting at the first move
    //     std::remove(m_grid[particle.cell_idx].begin(), m_grid[particle.cell_idx].end(), particle.ID),
    //     m_grid[particle.cell_idx].end()
    // );
}

void Solver::binParticles(bool use_temp_grid)
{
    if (use_temp_grid)
    {
        std::vector<std::vector<std::size_t>> temp_grid(m_grid.size());
    
        for (Particle& particle : m_particles)
        {
            std::size_t cell_idx { particlePositionToCellIndex(particle) };
            temp_grid[cell_idx].emplace_back(particle.ID);
            particle.cell_idx = cell_idx;
            particle.cell_idx_idx = temp_grid[cell_idx].size();
        }

        m_grid = std::move(temp_grid);
    }
    else
    {
        for (Particle& particle : m_particles)
        {
            std::size_t cell_idx { particlePositionToCellIndex(particle) };

            // only rebin particles that moved to a new cell
            if (cell_idx != particle.cell_idx)
            {
                // m_grid[cell_idx].erase(m_grid[cell_idx].begin() + particle.cell_idx_idx);
                removeParticleFromCell(particle);

                particle.cell_idx = cell_idx;
                particle.cell_idx_idx = m_grid[cell_idx].size();
                m_grid[cell_idx].emplace_back(particle.ID);
            }
        }
    }
}

std::array<std::size_t, 2> Solver::unflattenGridIndices(const std::size_t& idx)
{
    return std::array<std::size_t, 2>{static_cast<std::size_t>(idx / (m_max_cell_idx+1)), idx % (m_max_cell_idx+1)};
}

std::size_t Solver::flattenGridIndices(const std::size_t& x_idx, const std::size_t& y_idx)
{
    return x_idx * (m_max_cell_idx + 1) + y_idx;
}

void Solver::resolveParticleCollisions(Particle& particle, Particle& other_particle)
{
    double min_distance { particle.radius + other_particle.radius };
    if ( (particle.position - other_particle.position).mag2() < min_distance * min_distance )
    {
        Point2D r { particle.position - other_particle.position };
        Point2D n { r.normalized() };
        
        // adjust positions so there is no overlap
        double overlap { particle.radius + other_particle.radius - r.magnitude() };

        // introduce fudge factor so particles can't stack themselves out of the domain
        Point2D fudge_vec { applyFudgeFactor(n.perp(), 1e-6) };

        // particle.position       += overlap/2 * n + fudge_vec;
        // other_particle.position -= overlap/2 * n - fudge_vec;
        /*
        particle.position       += overlap/2 * n;
        other_particle.position -= overlap/2 * n;
        */

        double mass_ratio_i { 2 * other_particle.mass / (particle.mass + other_particle.mass) };
        double mass_ratio_j { 2 *       particle.mass / (particle.mass + other_particle.mass) };
        
        
        particle.position       += overlap/2 * n * mass_ratio_i * m_damping + fudge_vec;
        other_particle.position -= overlap/2 * n * mass_ratio_j * m_damping - fudge_vec;
        
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

void Solver::resolveCollisions()
{
    for (std::size_t cell_idx = 0; cell_idx < m_grid.size(); ++cell_idx)
    {
        for (std::size_t neighbor_idx : m_neighbor_indices[cell_idx])
        {
            if (!m_grid[cell_idx].empty() && !m_grid[neighbor_idx].empty())
            {
                // Recall: each cell contains particle ID's which correspond to their indices in the m_particles vector
                for (std::size_t i = 0; i < m_grid[cell_idx].size(); ++i)
                {
                    Particle& particle { m_particles[m_grid[cell_idx][i]] };

                    for (std::size_t j = 0; j < m_grid[neighbor_idx].size(); ++j)
                    {
                        Particle& other_particle { m_particles[m_grid[neighbor_idx][j]] };
                        if (particle.ID != other_particle.ID)
                        {
                            for (int i = 0; i < max_collision_resolution_count; ++i) { resolveParticleCollisions(particle, other_particle); }
                        }
                    }
                }
            }
        }
    }
}

void Solver::applyBoundaryConditions(Particle& particle)
{
    particle.position.setX(std::clamp(particle.position.x(), particle.radius, m_window_size - particle.radius));
    particle.position.setY(std::clamp(particle.position.y(), particle.radius, m_window_size - particle.radius));
}

void Solver::applyNoVelocityBC(Particle& particle)
{
    if (particle.position.x() < particle.radius)
    {
        particle.position.setX(particle.radius);
        particle.velocity.setX(particle.radius - m_damping * (particle.velocity.x() - particle.radius));
    }
    else if (particle.position.x() > m_window_size - particle.radius)
    {
        particle.position.setX(m_window_size - particle.radius);
        particle.velocity.setX((m_window_size - particle.radius) - m_damping * (particle.velocity.x() - (m_window_size - particle.radius)));
    }

    if (particle.position.y() < particle.radius)
    {
        particle.position.setY(particle.radius);
        particle.velocity.setY(particle.radius - m_damping * (particle.velocity.y() - particle.radius));
    }
    else if (particle.position.y() > m_window_size - particle.radius)
    {
        particle.position.setY(m_window_size - particle.radius);
        particle.velocity.setY((m_window_size - particle.radius) - m_damping * (particle.velocity.y() - (m_window_size - particle.radius)));
    }
}

void Solver::updateAndRenderParticles(Graphics& graphics, sf::RenderTarget& window_target)
{
    const double& dt { m_dt / m_substeps };
    
    for (int substep = 0; substep < m_substeps; ++substep)
    {
        resolveCollisions();
        
        for (std::size_t i = 0; i < m_particles.size(); ++ i)
        {
            Particle& particle { m_particles[i] };

            // Point2D a_tot { m_a_global + m_config.mouse_force / particle.mass * graphics.m_mouse.particleInteraction(particle.position) };
            // here we clamp displacement so that "velocity" gets clamped: |dx/dt| < v_max ==> |dx| < v_max * dt = v_max / frame_rate
            Point2D a_tot { m_a_global + m_mouse_force * graphics.m_mouse.particleInteraction(particle.position) };
            // Point2D a_tot { m_a_global };

            
            // Point2D next_r { 2 * particle.position - particle.velocity + a_tot*dt*dt };
            Point2D next_r { particle.position + (particle.position - particle.velocity).clamped(1.5) + a_tot*dt*dt };
            particle.velocity = particle.position;
            particle.position = next_r;
            applyNoVelocityBC(particle);
            

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

            std::size_t cell_idx { particlePositionToCellIndex(particle) };
            // only rebin the particle if it moves to a new cell
            if (cell_idx != particle.cell_idx)
            {
                // m_grid[cell_idx].erase(m_grid[cell_idx].begin() + static_cast<std::vector<std::size_t>::difference_type>(particle.cell_idx_idx));

                removeParticleFromCell(particle);

                particle.cell_idx = cell_idx;
                particle.cell_idx_idx = m_grid[cell_idx].size();
                m_grid[cell_idx].emplace_back(particle.ID);
            }

            if (substep + 1 == m_substeps)
            {
                // since we're already looping through particles, might as well render it now too!
                renderParticle(particle, window_target);
            }
        }
    }
}

sf::Text Solver::cellNumberText(const std::size_t& idx, Graphics& graphics)
{
    float character_size { 16 };
    sf::Text text(graphics.m_font, std::to_string(idx));
    text.setCharacterSize(static_cast<unsigned int>(character_size)); // in pixels, not points!
    text.setFillColor(sf::Color::Green); // set the color
    
    const auto& [idx_x, idx_y] = unflattenGridIndices(idx);
    text.setPosition(sf::Vector2f
    (
        (static_cast<float>(idx_x) + 0.5f) * static_cast<float>(m_grid_spacing) - character_size / 2,
        (static_cast<float>(idx_y) + 0.5f) * static_cast<float>(m_grid_spacing) - character_size / 2
    ));

    return text;
}

void Solver::drawCells(sf::RenderTarget& target_window)
{
    for (std::size_t i = 0; i < m_grid.size(); ++i)
    {
        sf::RectangleShape boundingBox(sf::Vector2f
        (
            static_cast<float>(m_grid_spacing),
            static_cast<float>(m_grid_spacing)
        ));
        
        const auto& [idx_x, idx_y] = unflattenGridIndices(i);

        boundingBox.setPosition(sf::Vector2f
        (
            (static_cast<float>(idx_x)) * static_cast<float>(m_grid_spacing),
            (static_cast<float>(idx_y)) * static_cast<float>(m_grid_spacing)
        ));

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

        sf::RectangleShape boundingBox(sf::Vector2f
        (
            static_cast<float>(m_grid_spacing),
            static_cast<float>(m_grid_spacing)
        ));
        
        const auto& [idx_x, idx_y] = unflattenGridIndices(i);

        boundingBox.setPosition(sf::Vector2f
        (
            (static_cast<float>(idx_x)) * static_cast<float>(m_grid_spacing),
            (static_cast<float>(idx_y)) * static_cast<float>(m_grid_spacing)
        ));

        boundingBox.setOutlineColor(sf::Color::Green);
        boundingBox.setOutlineThickness(0.5f);
        boundingBox.setFillColor(sf::Color::Transparent);
        target_window.draw(boundingBox);
    }
}

void Solver::renderParticle(Particle& particle, sf::RenderTarget& window_target)
{
    sf::CircleShape circle(static_cast<float>(particle.radius));
    circle.setOrigin(sf::Vector2f(static_cast<float>(particle.radius), static_cast<float>(particle.radius)));
    circle.setPosition(sf::Vector2f(static_cast<float>(particle.position.x()), static_cast<float>(particle.position.y())));
    circle.setFillColor(particle.color);
    window_target.draw(circle);
}

void Solver::renderParticles(sf::RenderTarget& window_target)
{
    for (Particle& particle : m_particles)
    {
        renderParticle(particle, window_target);
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
        10.0f, // milliseconds
        m_numParticles, // particles
        m_mass,
        {static_cast<double>(m_window_size)/4.0, static_cast<double>(m_window_size)/8.0},
        {1.0, 1.0}
    );

    float elapsed_time { 0.0f };
    float fps_elapsed_time { 0.0f };
    float generator_elapsed_time { 0.0f };
    sf::Text fps_text { graphics.fpsText() };

    sf::RenderWindow window(sf::VideoMode({m_window_size, m_window_size}), "Verlet Simulation");
    window.setFramerateLimit(m_frame_rate);
    // sf::ContextSettings settings;
    // settings.antiAliasingLevel = 1;

    while (window.isOpen())
    {
        // check all the window's events that were triggered since the last iteration of the loop
        while (const std::optional event = window.pollEvent())
        {
            // "close requested" event: we close the window
            if (event->is<sf::Event::Closed>() || sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Escape)) { window.close(); }
        }
        
        if (generator.generateBool(generator_elapsed_time) && m_particles.size() < static_cast<std::size_t>(generator.m_numParticles))
        {
            generator_elapsed_time = 0.0f;

            for (int counter = 0; counter < m_numGenerators; ++counter)
            {
                Particle new_particle { generator.generate(m_particles.size(), {counter * 100.0, 0.0}) };
                // for simple Verlet
                new_particle.velocity = new_particle.position - generator.m_velocity * m_dt;
                binNewParticle(new_particle);
                m_particles.emplace_back(new_particle);
                renderParticle(new_particle, window);
            }
        }

        arrowKeyGravity();
        
        // clear the window with black color
        window.clear(sf::Color::Black);

        // must keep update after the window clear since update also renders the particles
        updateAndRenderParticles(graphics, window);

        // draw cell borders (for debugging purposes). Pass `graphics` if also want cell indices displayed
        // drawCells(window);

        if (fps_elapsed_time  >= 500.0f)
        {
            fps_elapsed_time = 0.0f;
            fps_text = graphics.fpsText();
        }
        elapsed_time = static_cast<float>(graphics.m_clock.restart().asMilliseconds());
        generator_elapsed_time += elapsed_time;
        fps_elapsed_time       += elapsed_time;
        
        renderParticles(window);
        window.draw(fps_text);
        window.draw(graphics.particleCount(m_particles.size()));
        drawMouse(window, graphics);
        
        window.display();
    }
}

