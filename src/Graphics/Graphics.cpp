#include "Graphics.hpp"

sf::Text Graphics::fpsText()
{
    // Measure frame time
    float fps { 1.0f / m_clock.getElapsedTime().asSeconds() };
    sf::Text text(m_font, "FPS: " + std::to_string(static_cast<int>(fps)));

    text.setCharacterSize(16); // in pixels, not points!
    text.setFillColor(sf::Color::Red); // set the color
    text.setPosition({10.0f, 10.0f}); // Position at top-left corner

    return text;
}

sf::Text Graphics::particleCount(const std::size_t& count)
{
    sf::Text text(m_font, "Particle Count: " + std::to_string(count));

    text.setCharacterSize(16); // in pixels, not points!
    text.setFillColor(sf::Color::Red); // set the color
    text.setPosition({10.0f, 26.0f}); // Position at top-left corner (offset vertically by FPS text size)

    return text;
}
