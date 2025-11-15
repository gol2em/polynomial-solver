#include "de_casteljau.h"

/**
 * @file de_casteljau.cpp
 * @brief Implementation of De Casteljau subdivision algorithm
 */

#include <vector>

namespace polynomial_solver {

DeCasteljau::DeCasteljau() {
    // TODO: Implement constructor if state is needed.
}

DeCasteljau::~DeCasteljau() {
    // TODO: Implement destructor if state is added.
}

double DeCasteljau::evaluate1D(const std::vector<double>& bernstein_coeffs,
                               double t)
{
    const std::size_t n = bernstein_coeffs.size();
    if (n == 0u) {
        return 0.0;
    }

    // Standard De Casteljau algorithm (univariate).
    std::vector<double> tmp(bernstein_coeffs.begin(), bernstein_coeffs.end());

    for (std::size_t r = 1; r < n; ++r) {
        const std::size_t upper = n - r;
        for (std::size_t i = 0; i < upper; ++i) {
            tmp[i] = (1.0 - t) * tmp[i] + t * tmp[i + 1];
        }
    }

    return tmp[0];
}

void DeCasteljau::subdivide1D(const std::vector<double>& bernstein_coeffs,
                              double t,
                              std::vector<double>& left,
                              std::vector<double>& right)
{
    const std::size_t n = bernstein_coeffs.size();
    if (n == 0u) {
        left.clear();
        right.clear();
        return;
    }

    const std::size_t degree = n - 1u;

    // Build the De Casteljau triangle.
    std::vector<std::vector<double>> b(degree + 1u);
    b[0] = bernstein_coeffs;

    for (std::size_t r = 1; r <= degree; ++r) {
        const std::size_t upper = degree - r + 1u;
        b[r].resize(upper);
        for (std::size_t i = 0; i < upper; ++i) {
            b[r][i] = (1.0 - t) * b[r - 1u][i] + t * b[r - 1u][i + 1u];
        }
    }

    left.resize(n);
    right.resize(n);

    // Left segment: control points b[j][0], j = 0..degree.
    // Right segment: control points b[degree-j][j], j = 0..degree.
    for (std::size_t j = 0; j <= degree; ++j) {
        left[j]  = b[j][0];
        right[j] = b[degree - j][j];
    }
}

double DeCasteljau::evaluateTensorProduct(const std::vector<unsigned int>& degrees,
                                          const std::vector<double>& bernstein_coeffs,
                                          const std::vector<double>& parameters)
{
    const std::size_t dim = degrees.size();

    // Univariate case: directly delegate to evaluate1D.
    if (dim == 1u && parameters.size() == 1u) {
        const unsigned int deg = degrees[0];
        const std::size_t count = static_cast<std::size_t>(deg + 1u);
        if (bernstein_coeffs.size() != count) {
            return 0.0;
        }
        return evaluate1D(bernstein_coeffs, parameters[0]);
    }

    // Bivariate tensor-product case.
    if (dim == 2u && parameters.size() == 2u) {
        const unsigned int deg_x = degrees[0];
        const unsigned int deg_y = degrees[1];
        const std::size_t nx = static_cast<std::size_t>(deg_x + 1u);
        const std::size_t ny = static_cast<std::size_t>(deg_y + 1u);

        if (bernstein_coeffs.size() != nx * ny) {
            return 0.0;
        }

        const double x = parameters[0];
        const double y = parameters[1];

        // First evaluate along the second dimension (y) for each fixed x index.
        std::vector<double> row_values(nx, 0.0);

        for (std::size_t ix = 0; ix < nx; ++ix) {
            const std::size_t offset = ix * ny;
            std::vector<double> row(bernstein_coeffs.begin() + static_cast<std::ptrdiff_t>(offset),
                                    bernstein_coeffs.begin() + static_cast<std::ptrdiff_t>(offset + ny));
            row_values[ix] = evaluate1D(row, y);
        }

        // Then evaluate the resulting 1D Bernstein polynomial in x.
        return evaluate1D(row_values, x);
    }

    // Higher-dimensional tensor-product evaluation is not implemented yet.
    (void)bernstein_coeffs;
    (void)parameters;
    return 0.0;
}

} // namespace polynomial_solver

