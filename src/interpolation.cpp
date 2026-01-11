/**
 * @file interpolation.cpp
 * @brief Implementation of polynomial interpolation from function samples.
 */

#include "core/interpolation.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace polynomial_solver {

// ========== Abscissae Generation ==========

std::vector<double> Interpolation::generateAbscissae(unsigned int n, AbscissaeType type) {
    std::vector<double> points(n);

    switch (type) {
        case AbscissaeType::UNIFORM:
            // Uniform points: x_i = i / (n - 1) for i = 0, ..., n-1
            for (unsigned int i = 0; i < n; ++i) {
                points[i] = (n > 1) ? static_cast<double>(i) / (n - 1) : 0.5;
            }
            break;

        case AbscissaeType::CHEBYSHEV:
            // Chebyshev nodes (zeros of T_n): x_i = (1 + cos(π(2i+1)/(2n))) / 2
            // Mapped to [0, 1]
            for (unsigned int i = 0; i < n; ++i) {
                double theta = M_PI * (2.0 * i + 1.0) / (2.0 * n);
                points[i] = (1.0 + std::cos(theta)) / 2.0;
            }
            // Sort in ascending order
            std::sort(points.begin(), points.end());
            break;

        case AbscissaeType::CHEBYSHEV_EXTREMA:
            // Chebyshev extrema (Clenshaw-Curtis): x_i = (1 + cos(πi/(n-1))) / 2
            // Mapped to [0, 1]
            for (unsigned int i = 0; i < n; ++i) {
                double theta = (n > 1) ? M_PI * i / (n - 1) : 0.0;
                points[i] = (1.0 + std::cos(theta)) / 2.0;
            }
            // Sort in ascending order
            std::sort(points.begin(), points.end());
            break;
    }

    return points;
}

// ========== Least Squares Solver ==========

std::vector<double> Interpolation::solveLeastSquares(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& b) {

    std::size_t m = A.size();      // Number of samples
    std::size_t n = A[0].size();   // Number of coefficients

    // Build A^T A and A^T b
    std::vector<std::vector<double>> ATA(n, std::vector<double>(n, 0.0));
    std::vector<double> ATb(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < m; ++k) {
                ATA[i][j] += A[k][i] * A[k][j];
            }
        }
        for (std::size_t k = 0; k < m; ++k) {
            ATb[i] += A[k][i] * b[k];
        }
    }

    // Gaussian elimination with partial pivoting
    std::vector<double> x(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        // Find pivot
        std::size_t pivot = i;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (std::abs(ATA[j][i]) > std::abs(ATA[pivot][i])) {
                pivot = j;
            }
        }

        if (std::abs(ATA[pivot][i]) < 1e-14) {
            continue;  // Singular, skip
        }

        // Swap rows
        if (pivot != i) {
            std::swap(ATA[i], ATA[pivot]);
            std::swap(ATb[i], ATb[pivot]);
        }

        // Eliminate
        for (std::size_t j = i + 1; j < n; ++j) {
            double factor = ATA[j][i] / ATA[i][i];
            for (std::size_t k = i; k < n; ++k) {
                ATA[j][k] -= factor * ATA[i][k];
            }
            ATb[j] -= factor * ATb[i];
        }
    }

    // Back substitution
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        if (std::abs(ATA[i][i]) < 1e-14) {
            x[i] = 0.0;
            continue;
        }
        x[i] = ATb[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            x[i] -= ATA[i][j] * x[j];
        }
        x[i] /= ATA[i][i];
    }

    return x;
}

// ========== 1D Interpolation ==========

Polynomial Interpolation::interpolate1D(
    const std::function<double(double)>& func,
    unsigned int degree,
    double x_min,
    double x_max,
    AbscissaeType abscissae) {

    unsigned int n = degree + 1;  // Number of sample points
    std::vector<double> nodes = generateAbscissae(n, abscissae);

    // Sample function
    std::vector<double> values(n);
    for (unsigned int i = 0; i < n; ++i) {
        double x = transformFromUnit(nodes[i], x_min, x_max);
        values[i] = func(x);
    }

    // Build Vandermonde matrix for power basis
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    for (unsigned int i = 0; i < n; ++i) {
        double s = nodes[i];
        double s_pow = 1.0;
        for (unsigned int j = 0; j <= degree; ++j) {
            A[i][j] = s_pow;
            s_pow *= s;
        }
    }

    // Solve for power coefficients
    std::vector<double> power_coeffs = solveLeastSquares(A, values);

    return Polynomial::fromPower({degree}, power_coeffs);
}

// ========== 2D Interpolation ==========

Polynomial Interpolation::interpolate2D(
    const std::function<double(double, double)>& func,
    unsigned int degree_x,
    unsigned int degree_y,
    double x_min,
    double x_max,
    double y_min,
    double y_max,
    AbscissaeType abscissae) {

    unsigned int nx = degree_x + 1;
    unsigned int ny = degree_y + 1;
    unsigned int n_coeffs = nx * ny;

    std::vector<double> nodes_x = generateAbscissae(nx, abscissae);
    std::vector<double> nodes_y = generateAbscissae(ny, abscissae);

    // Sample function at tensor product grid
    std::vector<double> values;
    std::vector<std::pair<double, double>> sample_points;

    for (unsigned int j = 0; j < ny; ++j) {
        for (unsigned int i = 0; i < nx; ++i) {
            double x = transformFromUnit(nodes_x[i], x_min, x_max);
            double y = transformFromUnit(nodes_y[j], y_min, y_max);
            sample_points.push_back({nodes_x[i], nodes_y[j]});
            values.push_back(func(x, y));
        }
    }

    // Build Vandermonde matrix for tensor-product power basis
    // Coefficient ordering: (0,0), (0,1), ..., (0,dy), (1,0), ..., (dx, dy)
    std::vector<std::vector<double>> A(values.size(), std::vector<double>(n_coeffs, 0.0));

    for (std::size_t k = 0; k < values.size(); ++k) {
        double s = sample_points[k].first;
        double t = sample_points[k].second;

        std::size_t coeff_idx = 0;
        for (unsigned int i = 0; i <= degree_x; ++i) {
            double s_pow = std::pow(s, static_cast<double>(i));
            for (unsigned int j = 0; j <= degree_y; ++j) {
                double t_pow = std::pow(t, static_cast<double>(j));
                A[k][coeff_idx] = s_pow * t_pow;
                coeff_idx++;
            }
        }
    }

    // Solve for power coefficients
    std::vector<double> power_coeffs = solveLeastSquares(A, values);

    return Polynomial::fromPower({degree_x, degree_y}, power_coeffs);
}

// ========== General ND Interpolation ==========

Polynomial Interpolation::interpolateND(
    const std::function<double(const std::vector<double>&)>& func,
    const std::vector<unsigned int>& degrees,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    AbscissaeType abscissae) {

    std::size_t dim = degrees.size();
    if (dim == 0 || lower_bounds.size() != dim || upper_bounds.size() != dim) {
        throw std::invalid_argument("Interpolation: dimension mismatch in bounds/degrees");
    }

    // Generate abscissae for each dimension
    std::vector<std::vector<double>> nodes(dim);
    for (std::size_t d = 0; d < dim; ++d) {
        nodes[d] = generateAbscissae(degrees[d] + 1, abscissae);
    }

    // Compute total number of coefficients and samples
    unsigned int n_coeffs = 1;
    unsigned int n_samples = 1;
    for (std::size_t d = 0; d < dim; ++d) {
        n_coeffs *= (degrees[d] + 1);
        n_samples *= (degrees[d] + 1);
    }

    // Generate tensor product sample points and evaluate function
    std::vector<std::vector<double>> sample_points;  // Normalized coordinates
    std::vector<double> values;

    // Use multi-index enumeration
    std::vector<unsigned int> idx(dim, 0);
    for (unsigned int k = 0; k < n_samples; ++k) {
        // Compute sample point
        std::vector<double> s(dim);  // Normalized
        std::vector<double> x(dim);  // Original domain
        for (std::size_t d = 0; d < dim; ++d) {
            s[d] = nodes[d][idx[d]];
            x[d] = transformFromUnit(s[d], lower_bounds[d], upper_bounds[d]);
        }

        sample_points.push_back(s);
        values.push_back(func(x));

        // Increment multi-index (little-endian style)
        for (std::size_t d = 0; d < dim; ++d) {
            idx[d]++;
            if (idx[d] <= degrees[d]) break;
            idx[d] = 0;
        }
    }

    // Build Vandermonde matrix for tensor-product power basis
    std::vector<std::vector<double>> A(n_samples, std::vector<double>(n_coeffs, 0.0));

    for (unsigned int k = 0; k < n_samples; ++k) {
        const std::vector<double>& s = sample_points[k];

        // Enumerate coefficient multi-index
        std::vector<unsigned int> coeff_idx(dim, 0);
        for (unsigned int c = 0; c < n_coeffs; ++c) {
            // Compute monomial value: s[0]^coeff_idx[0] * s[1]^coeff_idx[1] * ...
            double monomial = 1.0;
            for (std::size_t d = 0; d < dim; ++d) {
                monomial *= std::pow(s[d], static_cast<double>(coeff_idx[d]));
            }
            A[k][c] = monomial;

            // Increment coefficient multi-index
            for (std::size_t d = 0; d < dim; ++d) {
                coeff_idx[d]++;
                if (coeff_idx[d] <= degrees[d]) break;
                coeff_idx[d] = 0;
            }
        }
    }

    // Solve for power coefficients
    std::vector<double> power_coeffs = solveLeastSquares(A, values);

    return Polynomial::fromPower(degrees, power_coeffs);
}

} // namespace polynomial_solver

