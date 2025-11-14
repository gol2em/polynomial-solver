#ifndef GEOMETRY_H
#define GEOMETRY_H

/**
 * @file geometry.h
 * @brief Geometry tools for polynomial solver
 * 
 * This module provides geometric utilities for working with polynomial curves,
 * including point representations, intervals, and geometric operations.
 */

namespace polynomial_solver {

/**
 * @class Point2D
 * @brief Represents a 2D point
 */
class Point2D {
public:
    /**
     * @brief Default constructor
     */
    Point2D();
    
    /**
     * @brief Constructor with coordinates
     * @param x X-coordinate
     * @param y Y-coordinate
     */
    Point2D(double x, double y);
    
    /**
     * @brief Destructor
     */
    ~Point2D();
    
    // TODO: Add geometric operations
    // - Distance calculations
    // - Vector operations
    
private:
    double x_;
    double y_;
};

/**
 * @class Interval
 * @brief Represents a 1D interval
 */
class Interval {
public:
    /**
     * @brief Default constructor
     */
    Interval();
    
    /**
     * @brief Constructor with bounds
     * @param min Lower bound
     * @param max Upper bound
     */
    Interval(double min, double max);
    
    /**
     * @brief Destructor
     */
    ~Interval();
    
    // TODO: Add interval operations
    // - Intersection
    // - Union
    // - Contains check
    
private:
    double min_;
    double max_;
};

} // namespace polynomial_solver

#endif // GEOMETRY_H

