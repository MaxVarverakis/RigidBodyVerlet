#pragma once

#include <SFML/Graphics.hpp>
#include "../Points/Points.hpp"
#include "../Mouse/Mouse.hpp"

struct Graphics
{
    Mouse m_mouse;
    sf::Clock m_clock;
    sf::Font m_font;

    Graphics(float mouse_radius)
        : m_mouse(mouse_radius)
    {
        // initialize text font
        if (!m_font.openFromFile("Bitstream-Vera-Sans-Mono/VeraMono.ttf")) {
            std::cerr << "Failed to load font!" << std::endl;
        }
        
        m_clock.start();
    };

    sf::Text fpsText();

    sf::Text particleCount(const std::size_t& count);
};
