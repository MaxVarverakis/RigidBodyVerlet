#pragma once

#include <iostream>
#include <random>
#include "../ParticleGenerator/ParticleGenerator.hpp"

class Solver
{
private:
    const int m_substeps;
    uint m_window_size;
    std::vector<Particle> m_particles;

    std::size_t m_max_cell_idx;
    int m_max_cell_idx_int;
    double m_grid_spacing;
    std::vector<std::vector<std::size_t>> m_grid;
    std::vector<std::vector<std::size_t>> m_neighbor_indices;
    
    static constexpr int m_neighbor_offsets[5][2] { {0,0}, {0,1}, {1,0}, {1,1}, {1,-1} };
    
    const uint32_t m_frame_rate { 60 };
    const double m_dt { 1.0 / m_frame_rate };
    Point2D m_a_global { 0.0, 500.0 };
    const double m_damping { 0.75 };
    const double m_mass { 16.0 };
    const double m_mouse_force { 20.0 * m_a_global.magnitude() };
    const int max_collision_resolution_count { 1 };
    const int m_numParticles { 3000 };
    const int m_numGenerators { 5 };

public:
    Solver(const uint window_size);

    void precomputeNeighborIndices();

    double generateFudgeFactor(const double& max_fudgeness);

    Point2D applyFudgeFactor(const Point2D& n_perp, const double& max_fudgeness);

    std::size_t flattenGridIndices(const std::size_t& x_idx, const std::size_t& y_idx);

    std::array<std::size_t, 2> unflattenGridIndices(const std::size_t& idx);

    std::size_t particlePositionToCellIndex(Particle& particle);

    void binNewParticle(Particle& particle);

    void removeParticleFromCell(Particle& particle);

    void binParticles(bool use_temp_grid = false);

    void resolveParticleCollisions(Particle& particle, Particle& other_particle);
    
    // void resolveCellCollisions(const int cell_i, const int cell_j, const int rel_i, const int rel_j);

    void resolveCollisions();
    
    void applyBoundaryConditions(Particle& particle);

    void applyNoVelocityBC(Particle& particle);

    void updateAndRenderParticles(Graphics& graphics, sf::RenderTarget& window_target);

    sf::Text cellNumberText(const std::size_t& idx, Graphics& graphics);

    void drawCells(sf::RenderTarget& target_window);

    void drawCells(sf::RenderTarget& target_window, Graphics& graphics);

    void renderParticles(sf::RenderTarget& window_target);
    
    void renderParticle(Particle& particle, sf::RenderTarget& window_target);

    void arrowKeyGravity();

    void drawMouse(sf::RenderWindow& window,Graphics& graphics);

    void run();
};
