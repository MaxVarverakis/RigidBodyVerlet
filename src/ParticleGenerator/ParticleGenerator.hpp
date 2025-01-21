#pragma once

#include "../Points/Points.hpp"
#include "../Graphics/Graphics.hpp"

struct ParticleGenerator
{
    float m_generator_time_offset;
    int m_numParticles;
    const double m_mass;
    Point2D m_location;
    Point2D m_velocity;

    ParticleGenerator
    (
        const float generator_time_offset,
        const int numParticles,
        const double mass,
        const Point2D location,
        const Point2D velocity
    )
        : m_generator_time_offset {generator_time_offset}
        , m_numParticles {numParticles}
        , m_mass {mass}
        , m_location {location}
        , m_velocity {velocity}
    {}
    
    Particle generate(const std::size_t idx, const Point2D& location_offset = {0.0, 0.0});

    bool generateBool(const float& elapsed_time);

    sf::Color colorOnSpawnIterationLinear(const float& normalized_iteration);
};
