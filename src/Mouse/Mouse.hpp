#pragma once

#include <SFML/Graphics.hpp>
#include "../Points/Points.hpp"

struct Mouse
{
    sf::CircleShape m_shape;
    
    Mouse(const float& mouse_radius)
    {
        m_shape.setRadius(mouse_radius);
        m_shape.setOrigin({mouse_radius, mouse_radius});
        m_shape.setFillColor(sf::Color::Transparent);
        m_shape.setOutlineThickness(1);
        m_shape.setPointCount(32);
    }

    double x     () { return static_cast<double>(m_shape.getPosition().x); }
    double y     () { return static_cast<double>(m_shape.getPosition().y); }
    double radius() { return static_cast<double>(m_shape.getRadius()); }
    Point2D positionP2D() { return Point2D{this->x(), this->y()}; }

    void setCursorPosition(sf::Vector2i cursor_location) { m_shape.setPosition(sf::Vector2f(cursor_location)); }
    
    void setMouseRadiusColor(sf::Color color) { m_shape.setOutlineColor(color); }

    Point2D particleInteraction(const Point2D& position)
    {
        Point2D cursor_location { this->positionP2D() };
        double distance { position.distanceTo(cursor_location) };
        if
        (
            sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)
            && distance < this->radius()
        )
        {
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::LShift))
            {
                return (position - cursor_location).normalized(); // * distance / this->radius();
            }
            else
            {
                return (cursor_location - position).normalized(); // * distance / this->radius();
            }
        }
        else
        {
            return {0.0, 0.0};
        }
    }
};
