#include "ParticleGenerator.hpp"

Particle ParticleGenerator::generate(const std::size_t idx, const Point2D& location_offset)
{
    return
    {
        idx,
        sf::Color { colorOnSpawnIterationLinear(idx) },
        m_mass,
        m_location + location_offset,
        m_velocity,
    };
}

bool ParticleGenerator::generateBool(const float& elapsed_time)
{
    return elapsed_time >= m_generator_time_offset;
}

sf::Color ParticleGenerator::colorOnSpawnIterationLinear(const float& idx)
{
    const float normalized_idx { idx / m_numParticles };

    std::uint8_t red;
    std::uint8_t green;
    std::uint8_t blue;

    // if (normalized_idx == 0.0f)
    // {
    //     red   = 255;
    //     green = 255;
    //     blue  = 255;
    // }
    // First segment: Red to Green
    if (normalized_idx <= (1.0f / 3.0f))
    {
        float t = normalized_idx * 3.0f; // Scale to [0, 1]
        red = static_cast<std::uint8_t>((1.0f - t) * 255);
        green = static_cast<std::uint8_t>(t * 255);
        blue = 0;
    }
    // Second segment: Green to Blue
    else if (normalized_idx <= (2.0f / 3.0f))
    {
        float t = (normalized_idx - 1.0f / 3.0f) * 3.0f; // Scale to [0, 1]
        red = 0;
        green = static_cast<std::uint8_t>((1.0f - t) * 255);
        blue = static_cast<std::uint8_t>(t * 255);
    }
    // Third segment: Blue to Red
    else
    {
        float t = (normalized_idx - 2.0f / 3.0f) * 3.0f; // Scale to [0, 1]
        red = static_cast<std::uint8_t>(t * 255);
        green = 0;
        blue = static_cast<std::uint8_t>((1.0f - t) * 255);
    }
    
    return sf::Color(red, green, blue);
}
