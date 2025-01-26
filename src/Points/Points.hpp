#pragma once

#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>

class Point3D
{
private:
    double m_x;
    double m_y;
    double m_z;

public:
    Point3D(double x, double y, double z) : m_x{x}, m_y{y}, m_z{z} {};

    // Getters
    double x() const { return m_x; }
    double y() const { return m_y; }
    double z() const { return m_z; }

    // Setters
    void setX(double x) { m_x = x; }
    void setY(double y) { m_y = y; }
    void setZ(double z) { m_z = z; }

    // Methods
    double distanceTo(const Point3D& other) const
    {
        return std::sqrt(
              (m_x - other.x()) * (m_x - other.x())
            + (m_y - other.y()) * (m_y - other.y())
            + (m_z - other.z()) * (m_z - other.z())
        );
    }
    
    Point3D cross(const Point3D& other) const
    {
        return {
            m_y * other.z() - m_z * other.y(),
            m_z * other.x() - m_x * other.z(),
            m_x * other.y() - m_y * other.x()
        };
    }

    double dot(const Point3D& other) const
    {
        return (m_x * other.x() + m_y * other.y() + m_z * other.z());
    }

    Point3D operator-(const Point3D& other) const
    {
        return { m_x - other.x(), m_y - other.y(), m_z - other.z() };
    }

    Point3D operator-(const double& value) const
    {
        return { m_x - value, m_y - value, m_z - value };
    }
    
    Point3D operator+(const double& value) const
    {
        return { m_x + value, m_y + value, m_z + value };
    }

    Point3D operator+(const Point3D& other) const
    {
        return { m_x + other.x(), m_y + other.y(), m_z + other.z() };
    }

    Point3D operator*(const double scalar) const
    {
        return { scalar * m_x, scalar * m_y, scalar * m_z };
    }

    Point3D operator/(const double scalar) const
    {
        return { m_x / scalar, m_y / scalar, m_z / scalar };
    }

    bool operator==(const Point3D& other) const
    {
        return m_x == other.x() && m_y == other.y() && m_z == other.z();
    }
    
    bool operator!=(const double& value) const
    {
        return m_x != value || m_y != value || m_z != value;
    }

    void operator+=(const Point3D& other)
    {
        m_x += other.x();
        m_y += other.y();
        m_z += other.z();
    }

    void operator*=(const double& scalar)
    {
        m_x *= scalar;
        m_y *= scalar;
        m_z *= scalar;
    }

    void operator-=(const Point3D& other)
    {
        m_x -= other.x();
        m_y -= other.y();
        m_z -= other.z();
    }

    const Point3D normalized() const
    {
        return *this / (this->magnitude());
    }

    void normalize()
    {
        // modifies `this` Point3D by normalizing it

        double magnitude { this->magnitude() };
        m_x /= magnitude;
        m_y /= magnitude;
        m_z /= magnitude;
    }

    double mag2() const { return m_x*m_x + m_y*m_y + m_z*m_z; }

    double magnitude() const { return std::sqrt(m_x*m_x + m_y*m_y + m_z*m_z);}

    // Overload the << operator to output the Point3D in the format (x, y, z)
    friend std::ostream& operator<<(std::ostream& os, const Point3D& point) {
        os << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")";
        return os;
    }
};

// overload scalar multiplication, but to make it commutative we define it outside the Point2D class
Point3D operator*(const double scalar, const Point3D& point);

class Point2D
{
private:
    double m_x;
    double m_y;

public:
    Point2D(double x, double y) : m_x{x}, m_y{y} {};

    // Getters
    double x() const { return m_x; }
    double y() const { return m_y; }

    // Setters
    void setX(double x) { m_x = x; }
    void setY(double y) { m_y = y; }

    // Methods
    sf::Vector2f toVf()
    {
        return { static_cast<float>(m_x), static_cast<float>(m_y) };
    }

    double distanceTo(const Point2D& other) const
    {
        return std::sqrt(
              (m_x - other.x()) * (m_x - other.x())
            + (m_y - other.y()) * (m_y - other.y())
        );
    }

    void clamp(const double max_value)
    {
        if (m_x > max_value)
        {
            m_x = max_value;
        }
        else if (m_x < -max_value)
        {
            m_x = -max_value;
        }
        
        if (m_y > max_value)
        {
            m_y = max_value;
        }
        else if (m_y < -max_value)
        {
            m_y = -max_value;
        }
    }
    
    Point2D clamped(const double max_value) const
    {
        double new_x { m_x };
        double new_y { m_y };

        if (m_x > max_value)
        {
            new_x = max_value;
        }
        else if (m_x < -max_value)
        {
            new_x = -max_value;
        }
        
        if (m_y > max_value)
        {
            new_y = max_value;
        }
        else if (m_y < -max_value)
        {
            new_y = -max_value;
        }

        return { new_x, new_y };
    }

    double dot(const Point2D& other) const
    {
        return (m_x * other.x() + m_y * other.y());
    }

    Point2D perp() const
    {
        return {-m_y, m_x};
    }

    Point2D operator-(const Point2D& other) const
    {
        return { m_x - other.x(), m_y - other.y() };
    }

    Point2D operator-(const double& value) const
    {
        return { m_x - value, m_y - value };
    }
    
    Point2D operator+(const double& value) const
    {
        return { m_x + value, m_y + value };
    }

    Point2D operator+(const Point2D& other) const
    {
        return { m_x + other.x(), m_y + other.y() };
    }

    Point2D operator*(const double scalar) const
    {
        return { scalar * m_x, scalar * m_y };
    }

    Point2D operator/(const double scalar) const
    {
        return { m_x / scalar, m_y / scalar };
    }

    bool operator==(const Point2D& other) const
    {
        return m_x == other.x() && m_y == other.y();
    }
    
    bool operator!=(const double& value) const
    {
        return m_x != value || m_y != value;
    }

    void operator+=(const Point2D& other)
    {
        m_x += other.x();
        m_y += other.y();
    }

    void operator*=(const double& scalar)
    {
        m_x *= scalar;
        m_y *= scalar;
    }

    void operator-=(const Point2D& other)
    {
        m_x -= other.x();
        m_y -= other.y();
    }

    const Point2D normalized() const
    {
        return *this / (this->magnitude());
    }

    void normalize()
    {
        // modifies `this` Point2D by normalizing it

        double magnitude { this->magnitude() };
        m_x /= magnitude;
        m_y /= magnitude;
    }

    double mag2() const { return m_x*m_x + m_y*m_y; }

    double magnitude() const { return std::sqrt(m_x*m_x + m_y*m_y);}

    // Overload the << operator to output the Point2D in the format (x, y, z)
    friend std::ostream& operator<<(std::ostream& os, const Point2D& point) {
        os << "(" << point.x() << ", " << point.y() << ")";
        return os;
    }
};

// overload scalar multiplication, but to make it commutative we define it outside the Point2D class
Point2D operator*(const double scalar, const Point2D& point);

struct Particle
{
    const std::size_t ID;
    sf::Color color;
    const double mass;
    Point2D position;
    Point2D velocity;
    
    std::size_t cell_idx { 0 };
    std::size_t cell_idx_idx { 0 };
    
    // const double radius { (mass > 1.0) ? std::sqrt(mass) : mass };
    const double radius { std::sqrt(mass) };

    bool operator==(const Particle& other) const
    {
        return ID == other.ID;
    }
};
