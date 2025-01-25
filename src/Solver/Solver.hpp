#pragma once

#include <iostream>
#include <random>
#include "../ParticleGenerator/ParticleGenerator.hpp"

class Solver
{
private:
    const int m_substeps;
    uint m_window_size;
    Particles m_particles;

    std::size_t m_max_cell_idx;
    int m_max_cell_idx_int;
    double m_grid_spacing;
    std::vector<std::vector<std::size_t>> m_grid;
    std::vector<std::vector<std::size_t>> m_neighbor_indices;
    
    static constexpr int m_neighbor_offsets[5][2] { {0,0}, {0,1}, {1,0}, {1,1}, {1,-1} };
    
    double m_cached_fudge_factors[10];

    const uint32_t m_frame_rate { 120 };
    const double m_dt { 1.0 / m_frame_rate };
    Point2D m_a_global { 0.0, 500.0 };
    const double m_damping { 0.75 };
    const double m_mass { 25.0 };
    const double m_mouse_force { 20.0 * m_a_global.magnitude() };
    const int max_collision_resolution_count { 1 };
    const std::size_t m_numParticles { 3000 };
    const int m_numGenerators { 5 };

public:
    Solver(const uint window_size);

    void precomputeNeighborIndices();

    double generateFudgeFactor(const double& max_fudgeness);

    double getFudgeFactor(const double& hash_key);

    Point2D applyFudgeFactor(const Point2D& n_perp);

    std::size_t flattenGridIndices(const std::size_t& x_idx, const std::size_t& y_idx);

    std::array<std::size_t, 2> unflattenGridIndices(const std::size_t& idx);

    std::size_t particlePositionToCellIndex(std::size_t particle_ID);

    void binNewParticle(std::size_t particle_ID);

    void removeParticleFromCell(std::size_t particle_ID);

    void binParticles();

    void resolveParticleCollisions(std::size_t particle_ID, std::size_t& other_particle_ID);

    void resolveCollisions();
    
    void applyBoundaryConditions(std::size_t particle_ID);

    void applyNoVelocityBC(std::size_t particle_ID);

    void updateAndRenderParticles(Graphics& graphics, sf::RenderTarget& window_target);

    sf::Text cellNumberText(const std::size_t& idx, Graphics& graphics);

    void drawCells(sf::RenderTarget& target_window);

    void drawCells(sf::RenderTarget& target_window, Graphics& graphics);

    void renderParticles(sf::RenderTarget& window_target);
    
    void renderParticle(std::size_t particle_ID, sf::RenderTarget& window_target);

    void arrowKeyGravity();

    void drawMouse(sf::RenderWindow& window,Graphics& graphics);

    void run();
};
