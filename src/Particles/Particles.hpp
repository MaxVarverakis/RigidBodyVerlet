#pragma once

#include <SFML/Graphics.hpp>
#include "../Points/Points.hpp"

struct alignas(64) Particles
{
    std::size_t size { 0 };
    
    std::vector<Point2D>      position;
    std::vector<Point2D>      position_previous;

    std::vector<double>      mass;
    std::vector<double>      radius;
    std::vector<sf::Color>   color;
    std::vector<std::size_t> ID;
    std::vector<std::size_t> cell_idx;
    std::vector<std::size_t> cell_idx_idx;

    Particles(const std::size_t& max_particles)
    {
        position.reserve(max_particles);
        position_previous.reserve(max_particles);
        mass.reserve(max_particles);
        radius.reserve(max_particles);
        color.reserve(max_particles);
        ID.reserve(max_particles);
        cell_idx.reserve(max_particles);
        cell_idx_idx.reserve(max_particles);
    }

    // Point2D position(std::size_t index) const
    // {
    //     return { x[index], y[index] };
    // }
    
    // Point2D previousPosition(std::size_t index) const
    // {
    //     return { x_previous[index], y_previous[index] };
    // }

    // void setPosition(std::size_t index, const Point2D& position)
    // {
    //     x[index] = position.x();
    //     y[index] = position.y();
    // }
    
    // void addPosition(std::size_t index, const Point2D& position)
    // {
    //     x[index] += position.x();
    //     y[index] += position.y();
    // }
    
    // void subtractPosition(std::size_t index, const Point2D& position)
    // {
    //     x[index] -= position.x();
    //     y[index] -= position.y();
    // }
    
    // void setPreviousPosition(std::size_t index, const Point2D& position)
    // {
    //     x_previous[index] = position.x();
    //     y_previous[index] = position.y();
    // }

    // void addPreviousPosition(std::size_t index, const Point2D& position)
    // {
    //     x_previous[index] += position.x();
    //     y_previous[index] += position.y();
    // }
    
    // void subtractPreviousPosition(std::size_t index, const Point2D& position)
    // {
    //     x_previous[index] -= position.x();
    //     y_previous[index] -= position.y();
    // }
};
