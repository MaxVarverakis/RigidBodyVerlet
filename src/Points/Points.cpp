#include "Points.hpp"

Point2D operator*(const double scalar, const Point2D& point)
{
    return Point2D(scalar * point.x(), scalar * point.y());
}

Point3D operator*(const double scalar, const Point3D& point)
{
    return Point3D(scalar * point.x(), scalar * point.y(), scalar * point.z());
}
