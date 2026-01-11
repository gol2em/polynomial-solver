#ifndef GEOMETRY_BASE_H
#define GEOMETRY_BASE_H

/**
 * @file geometry_base.h
 * @brief Templated geometry tools for polynomial solver
 *
 * This module provides templated geometric utilities that work with any
 * scalar type (double, mpfr::mpreal, etc.). The tolerance for geometric
 * operations is automatically adjusted based on the scalar type's precision.
 *
 * For dimensions > 2, only bounding box operations are fully supported.
 * Convex hull algorithms are implemented for 2D and 3D only.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

// Include Boost BEFORE opening polynomial_solver namespace
#ifdef USE_MPFR_BACKEND
#include <boost/multiprecision/mpfr.hpp>
#endif

#ifdef ENABLE_QUADMATH
extern "C" {
    #include <quadmath.h>
}
#endif

namespace polynomial_solver {

//=============================================================================
// Tolerance traits - compute appropriate epsilon for each scalar type
//=============================================================================

/**
 * @brief Traits class for computing geometric tolerance based on scalar type.
 *
 * For double: uses ~1e-15 (machine epsilon * small factor)
 * For mpfr::mpreal: uses 10^(-precision_bits * log10(2) + 5)
 * For __float128: uses ~1e-30
 */
template<typename Scalar>
struct GeometryToleranceTraits {
    static Scalar epsilon() {
        // Default: use numeric_limits if available
        return Scalar(100) * std::numeric_limits<Scalar>::epsilon();
    }
};

// Specialization for double
template<>
struct GeometryToleranceTraits<double> {
    static double epsilon() {
        return 1e-15;  // Same as original GEOMETRY_EPSILON
    }
};

// Specialization for float
template<>
struct GeometryToleranceTraits<float> {
    static float epsilon() {
        return 1e-6f;
    }
};

// Specialization for long double
template<>
struct GeometryToleranceTraits<long double> {
    static long double epsilon() {
        return 1e-18L;
    }
};

#ifdef ENABLE_QUADMATH
template<>
struct GeometryToleranceTraits<__float128> {
    static __float128 epsilon() {
        // Use strtoflt128 instead of Q suffix for portability
        return strtoflt128("1e-30", nullptr);
    }
};
#endif

// MPFR specialization (if available)
#ifdef USE_MPFR_BACKEND
template<>
struct GeometryToleranceTraits<boost::multiprecision::mpfr_float> {
    static boost::multiprecision::mpfr_float epsilon() {
        // Get current precision in bits, convert to decimal digits
        // eps ≈ 2^(-precision) ≈ 10^(-precision * log10(2))
        // We use a slightly larger tolerance (add 5 digits margin)
        unsigned int bits = boost::multiprecision::mpfr_float::default_precision();
        unsigned int decimal_digits = static_cast<unsigned int>(bits * 0.30103);  // log10(2)
        unsigned int eps_digits = (decimal_digits > 10) ? decimal_digits - 5 : 5;

        boost::multiprecision::mpfr_float eps("1e-" + std::to_string(eps_digits));
        return eps;
    }
};
#endif

//=============================================================================
// Machine epsilon traits - compute appropriate machine epsilon for each type
//=============================================================================

/**
 * @brief Traits class for computing machine epsilon based on scalar type.
 *
 * This is the smallest value such that 1 + epsilon != 1.
 * Used for clamping box widths to avoid reporting zero error.
 *
 * For double: uses std::numeric_limits<double>::epsilon() ≈ 2.2e-16
 * For mpfr::mpreal: uses 2^(-precision)
 * For __float128: uses ~1.9e-34
 */
template<typename Scalar>
struct MachineEpsilonTraits {
    static Scalar epsilon() {
        // Default: use numeric_limits if available
        return std::numeric_limits<Scalar>::epsilon();
    }
};

// Specialization for double (explicit for clarity)
template<>
struct MachineEpsilonTraits<double> {
    static double epsilon() {
        return std::numeric_limits<double>::epsilon();
    }
};

// Specialization for float
template<>
struct MachineEpsilonTraits<float> {
    static float epsilon() {
        return std::numeric_limits<float>::epsilon();
    }
};

// Specialization for long double
template<>
struct MachineEpsilonTraits<long double> {
    static long double epsilon() {
        return std::numeric_limits<long double>::epsilon();
    }
};

#ifdef ENABLE_QUADMATH
template<>
struct MachineEpsilonTraits<__float128> {
    static __float128 epsilon() {
        // FLT128_EPSILON ≈ 1.93e-34, but the Q suffix requires -fext-numeric-literals
        // Use strtoflt128 for portability
        return strtoflt128("1.92592994438723585305597794258492732e-34", nullptr);
    }
};
#endif

// MPFR specialization (if available)
#ifdef USE_MPFR_BACKEND
template<>
struct MachineEpsilonTraits<boost::multiprecision::mpfr_float> {
    static boost::multiprecision::mpfr_float epsilon() {
        // Machine epsilon for MPFR is 2^(1-precision)
        unsigned int bits = boost::multiprecision::mpfr_float::default_precision();
        return boost::multiprecision::ldexp(
            boost::multiprecision::mpfr_float(1),
            static_cast<int>(1 - bits));
    }
};
#endif

//=============================================================================
// Basic templated geometry structures
//=============================================================================

/**
 * @struct ConvexPolyhedronBase
 * @brief Generic convex polyhedron in R^{d}, represented by its vertices.
 *
 * @tparam Scalar The scalar type for coordinates (double, mpfr::mpreal, etc.)
 */
template<typename Scalar>
struct ConvexPolyhedronBase {
    std::vector<std::vector<Scalar>> vertices;
    std::size_t intrinsic_dim;

    ConvexPolyhedronBase() : intrinsic_dim(0) {}

    /// Returns the ambient dimension (dimension of the space)
    std::size_t ambient_dimension() const {
        return vertices.empty() ? 0u : vertices.front().size();
    }

    /// Returns the intrinsic dimension (dimension of the polytope itself)
    std::size_t dimension() const {
        return intrinsic_dim;
    }
};

/**
 * @struct ConvexPolyhedronBoxBase
 * @brief Simple axis-aligned bounding box in R^{d}.
 *
 * @tparam Scalar The scalar type for coordinates
 */
template<typename Scalar>
struct ConvexPolyhedronBoxBase {
    std::vector<Scalar> min_coords;
    std::vector<Scalar> max_coords;

    std::size_t dimension() const { return min_coords.size(); }
};

//=============================================================================
// Type aliases for convenience
//=============================================================================

using ConvexPolyhedronDouble = ConvexPolyhedronBase<double>;
using ConvexPolyhedronBoxDouble = ConvexPolyhedronBoxBase<double>;

//=============================================================================
// Forward declarations of templated algorithms
//=============================================================================

template<typename Scalar>
std::size_t compute_intrinsic_dimension_impl(
    const std::vector<std::vector<Scalar>>& points,
    Scalar eps);

template<typename Scalar>
ConvexPolyhedronBase<Scalar>
convex_hull_impl(const std::vector<std::vector<Scalar>>& points);

template<typename Scalar>
bool intersect_convex_polyhedra_impl(
    const std::vector<ConvexPolyhedronBase<Scalar>>& polyhedra,
    ConvexPolyhedronBase<Scalar>& intersection);

template<typename Scalar>
bool intersect_convex_polyhedron_with_last_coordinate_zero_impl(
    const ConvexPolyhedronBase<Scalar>& polyhedron,
    ConvexPolyhedronBase<Scalar>& intersection);

template<typename Scalar>
ConvexPolyhedronBoxBase<Scalar>
bounding_box_impl(const ConvexPolyhedronBase<Scalar>& polyhedron);

} // namespace polynomial_solver

// Include template implementation
#include "core/geometry_impl.h"

#endif // GEOMETRY_BASE_H

