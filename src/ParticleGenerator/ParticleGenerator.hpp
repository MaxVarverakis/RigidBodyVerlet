#pragma once


#include "../Points/Points.hpp"
#include "../Particles/Particles.hpp"
#include "../Graphics/Graphics.hpp"

struct ParticleGenerator
{
    float m_generator_time_offset;
    std::size_t m_numParticles;
    const double m_mass;
    Point2D m_location;
    Point2D m_velocity;

    ParticleGenerator
    (
        const float generator_time_offset,
        const std::size_t numParticles,
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
    
    void generate(Particles& particles, const double dt, const Point2D& location_offset = {0.0, 0.0});

    bool generateBool(const float& elapsed_time);

    sf::Color colorOnSpawnIterationLinear(std::size_t idx);
};
