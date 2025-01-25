#include "ParticleGenerator.hpp"

void ParticleGenerator::generate(Particles& particles, const double dt, const Point2D& location_offset)
{    
    const Point2D& location          { m_location + location_offset };
    const Point2D& previous_location { location - m_velocity * dt };

    particles.position.emplace_back(         location);
    particles.position_previous.emplace_back(previous_location);
    particles.mass.emplace_back(             m_mass);
    particles.radius.emplace_back(           std::sqrt(m_mass));
    particles.color.emplace_back(            colorOnSpawnIterationLinear(particles.size));
    particles.ID.emplace_back(               particles.size);
    particles.cell_idx.emplace_back(         0);
    particles.cell_idx_idx.emplace_back(     0);

}

bool ParticleGenerator::generateBool(const float& elapsed_time)
{
    return elapsed_time >= m_generator_time_offset;
}

sf::Color ParticleGenerator::colorOnSpawnIterationLinear(std::size_t idx)
{
    const float normalized_idx { static_cast<float>(idx) / m_numParticles };

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
