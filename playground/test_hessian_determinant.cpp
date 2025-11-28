// Test: Find zero set of Hessian determinant to determine convexity regions
// Function: f(x,y) = (10*x*y*x*y + sqrt(x*x*y)) + atan2(0.001*(sin(5*y)-2*x))
//                    + (10*y*y*y + x*x*x) + atan2(0.01*(sin(5*x)-2*y))
// Domain: [-1,1]^2 -> transform to [0,1]^2
// Goal: Find where det(Hessian) = 0 (boundary between convex/concave/saddle regions)

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

using namespace polynomial_solver;

// The original function on [-1,1]^2
double f_original(double x, double y) {
    // f1 = 10*x*y*x*y + sqrt(x*x*y)
    double f1 = 10.0 * x * y * x * y + std::sqrt(std::abs(x * x * y));
    
    // f2 = atan2(0.001, sin(5*y) - 2*x)
    double f2 = std::atan2(0.001, std::sin(5.0 * y) - 2.0 * x);
    
    // f3 = 10*y*y*y + x*x*x
    double f3 = 10.0 * y * y * y + x * x * x;
    
    // f4 = atan2(0.01, sin(5*x) - 2*y)
    double f4 = std::atan2(0.01, std::sin(5.0 * x) - 2.0 * y);
    
    return f1 + f2 + f3 + f4;
}

// Transform from [0,1]^2 to [-1,1]^2: x_orig = 2*u - 1, y_orig = 2*v - 1
double f_transformed(double u, double v) {
    double x = 2.0 * u - 1.0;
    double y = 2.0 * v - 1.0;
    return f_original(x, y);
}

// Numerical derivative (central difference)
double derivative(double (*func)(double, double), double u, double v, int var, double h = 1e-6) {
    if (var == 0) {
        return (func(u + h, v) - func(u - h, v)) / (2.0 * h);
    } else {
        return (func(u, v + h) - func(u, v - h)) / (2.0 * h);
    }
}

// Numerical second derivative
double second_derivative(double (*func)(double, double), double u, double v, int var1, int var2, double h = 1e-6) {
    if (var1 == 0 && var2 == 0) {
        // d²f/du²
        return (func(u + h, v) - 2.0 * func(u, v) + func(u - h, v)) / (h * h);
    } else if (var1 == 1 && var2 == 1) {
        // d²f/dv²
        return (func(u, v + h) - 2.0 * func(u, v) + func(u, v - h)) / (h * h);
    } else {
        // d²f/dudv (mixed derivative)
        return (func(u + h, v + h) - func(u + h, v - h) - func(u - h, v + h) + func(u - h, v - h)) / (4.0 * h * h);
    }
}

// Compute Hessian determinant at a point
double hessian_determinant(double u, double v) {
    double fuu = second_derivative(f_transformed, u, v, 0, 0);
    double fvv = second_derivative(f_transformed, u, v, 1, 1);
    double fuv = second_derivative(f_transformed, u, v, 0, 1);

    return fuu * fvv - fuv * fuv;
}

// Verify Hessian computation by checking at a few test points
void verify_hessian_computation() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Verifying Hessian Computation" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Test at a few points
    double test_points[][2] = {
        {0.5, 0.5},  // center
        {0.25, 0.25},
        {0.75, 0.75},
        {0.1, 0.9}
    };

    for (int i = 0; i < 4; ++i) {
        double u = test_points[i][0];
        double v = test_points[i][1];

        double f_val = f_transformed(u, v);
        double fuu = second_derivative(f_transformed, u, v, 0, 0);
        double fvv = second_derivative(f_transformed, u, v, 1, 1);
        double fuv = second_derivative(f_transformed, u, v, 0, 1);
        double det_H = hessian_determinant(u, v);

        std::cout << "Point (u=" << u << ", v=" << v << "):" << std::endl;
        std::cout << "  f = " << f_val << std::endl;
        std::cout << "  f_uu = " << fuu << std::endl;
        std::cout << "  f_vv = " << fvv << std::endl;
        std::cout << "  f_uv = " << fuv << std::endl;
        std::cout << "  det(H) = f_uu * f_vv - f_uv² = " << det_H << std::endl;
        std::cout << "  Verification: " << fuu << " * " << fvv << " - " << fuv << "² = "
                  << (fuu * fvv) << " - " << (fuv * fuv) << " = " << det_H << std::endl;
        std::cout << std::endl;
    }
}

// Proper polynomial interpolation using least-squares fit
// Sample at many points and fit a polynomial in power basis
Polynomial interpolate_proper(double (*func)(double, double),
                              unsigned int degree_u, unsigned int degree_v,
                              double u_min, double u_max, double v_min, double v_max) {
    // For proper interpolation, we need to solve a linear system
    // We'll use oversampling: sample at (degree+1)^2 points and use least-squares

    unsigned int n_samples_u = degree_u + 1;
    unsigned int n_samples_v = degree_v + 1;
    unsigned int n_coeffs = (degree_u + 1) * (degree_v + 1);

    // Sample function at uniform grid
    std::vector<double> sample_values;
    std::vector<std::pair<double, double>> sample_points;

    for (unsigned int j = 0; j < n_samples_v; ++j) {
        for (unsigned int i = 0; i < n_samples_u; ++i) {
            double s = static_cast<double>(i) / (n_samples_u - 1);
            double t = static_cast<double>(j) / (n_samples_v - 1);
            double u = u_min + (u_max - u_min) * s;
            double v = v_min + (v_max - v_min) * t;

            sample_points.push_back({s, t});
            sample_values.push_back(func(u, v));
        }
    }

    // Build Vandermonde-like matrix for 2D polynomial
    // For simplicity, we'll use a direct approach: fit in power basis
    std::vector<std::vector<double>> A(sample_values.size(), std::vector<double>(n_coeffs, 0.0));

    for (std::size_t k = 0; k < sample_values.size(); ++k) {
        double s = sample_points[k].first;
        double t = sample_points[k].second;

        std::size_t coeff_idx = 0;
        for (unsigned int i = 0; i <= degree_u; ++i) {
            for (unsigned int j = 0; j <= degree_v; ++j) {
                double s_pow = 1.0;
                for (unsigned int p = 0; p < i; ++p) s_pow *= s;
                double t_pow = 1.0;
                for (unsigned int q = 0; q < j; ++q) t_pow *= t;

                A[k][coeff_idx] = s_pow * t_pow;
                coeff_idx++;
            }
        }
    }

    // Solve least-squares: A^T A x = A^T b
    std::vector<std::vector<double>> ATA(n_coeffs, std::vector<double>(n_coeffs, 0.0));
    std::vector<double> ATb(n_coeffs, 0.0);

    for (std::size_t i = 0; i < n_coeffs; ++i) {
        for (std::size_t j = 0; j < n_coeffs; ++j) {
            for (std::size_t k = 0; k < sample_values.size(); ++k) {
                ATA[i][j] += A[k][i] * A[k][j];
            }
        }
        for (std::size_t k = 0; k < sample_values.size(); ++k) {
            ATb[i] += A[k][i] * sample_values[k];
        }
    }

    // Simple Gaussian elimination (for small systems)
    std::vector<double> power_coeffs(n_coeffs, 0.0);

    for (std::size_t i = 0; i < n_coeffs; ++i) {
        // Find pivot
        std::size_t pivot = i;
        for (std::size_t j = i + 1; j < n_coeffs; ++j) {
            if (std::abs(ATA[j][i]) > std::abs(ATA[pivot][i])) {
                pivot = j;
            }
        }

        if (std::abs(ATA[pivot][i]) < 1e-12) {
            continue;  // Singular, skip
        }

        // Swap rows
        if (pivot != i) {
            std::swap(ATA[i], ATA[pivot]);
            std::swap(ATb[i], ATb[pivot]);
        }

        // Eliminate
        for (std::size_t j = i + 1; j < n_coeffs; ++j) {
            double factor = ATA[j][i] / ATA[i][i];
            for (std::size_t k = i; k < n_coeffs; ++k) {
                ATA[j][k] -= factor * ATA[i][k];
            }
            ATb[j] -= factor * ATb[i];
        }
    }

    // Back substitution
    for (int i = n_coeffs - 1; i >= 0; --i) {
        if (std::abs(ATA[i][i]) < 1e-12) {
            power_coeffs[i] = 0.0;
            continue;
        }

        power_coeffs[i] = ATb[i];
        for (std::size_t j = i + 1; j < n_coeffs; ++j) {
            power_coeffs[i] -= ATA[i][j] * power_coeffs[j];
        }
        power_coeffs[i] /= ATA[i][i];
    }

    std::vector<unsigned int> degrees = {degree_u, degree_v};
    return Polynomial::fromPower(degrees, power_coeffs);
}

// Simple approximation: just use function values as Bernstein coefficients
Polynomial interpolate_bernstein_approx(double (*func)(double, double),
                                       unsigned int degree_u, unsigned int degree_v,
                                       double u_min, double u_max, double v_min, double v_max) {
    std::vector<double> bernstein_coeffs;

    for (unsigned int j = 0; j <= degree_v; ++j) {
        for (unsigned int i = 0; i <= degree_u; ++i) {
            double s = static_cast<double>(i) / degree_u;
            double t = static_cast<double>(j) / degree_v;
            double u = u_min + (u_max - u_min) * s;
            double v = v_min + (v_max - v_min) * t;

            double value = func(u, v);
            bernstein_coeffs.push_back(value);
        }
    }

    std::vector<unsigned int> degrees = {degree_u, degree_v};
    return Polynomial::fromBernstein(degrees, bernstein_coeffs);
}

int main(int argc, char* argv[]) {
    std::cout << "========================================" << std::endl;
    std::cout << "Hessian Determinant Zero Set Test" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Parse command line arguments
    unsigned int degree = 10;
    double degeneracy_multiplier = 10.0;
    unsigned int subdivisions = 4;  // Subdivide domain into k x k regions

    if (argc > 1) {
        degree = std::atoi(argv[1]);
    }
    if (argc > 2) {
        degeneracy_multiplier = std::atof(argv[2]);
    }
    if (argc > 3) {
        subdivisions = std::atoi(argv[3]);
    }

    std::cout << "Function: f(x,y) on [-1,1]^2" << std::endl;
    std::cout << "  f1 = 10*x*y*x*y + sqrt(|x*x*y|)" << std::endl;
    std::cout << "  f2 = atan2(0.001, sin(5*y) - 2*x)" << std::endl;
    std::cout << "  f3 = 10*y^3 + x^3" << std::endl;
    std::cout << "  f4 = atan2(0.01, sin(5*x) - 2*y)" << std::endl;
    std::cout << "  f = f1 + f2 + f3 + f4\n" << std::endl;
    
    std::cout << "Goal: Find zero set of det(Hessian(f))" << std::endl;
    std::cout << "  det(H) = f_uu * f_vv - f_uv^2\n" << std::endl;
    
    std::cout << "Bernstein interpolation degree: " << degree << "x" << degree << std::endl;
    std::cout << "Degeneracy multiplier: " << degeneracy_multiplier << std::endl;
    std::cout << "Domain subdivisions: " << subdivisions << "x" << subdivisions << " = "
              << (subdivisions * subdivisions) << " regions\n" << std::endl;

    // Test interpolation quality on a small region first
    if (subdivisions > 1) {
        std::cout << "========================================" << std::endl;
        std::cout << "Testing Interpolation Quality on Sample Region" << std::endl;
        std::cout << "========================================\n" << std::endl;

        // Pick first region: [0, 1/subdivisions] x [0, 1/subdivisions]
        double test_u_min = 0.0, test_u_max = 1.0 / subdivisions;
        double test_v_min = 0.0, test_v_max = 1.0 / subdivisions;

        std::cout << "Test region: u=[" << test_u_min << "," << test_u_max
                  << "], v=[" << test_v_min << "," << test_v_max << "]" << std::endl;

        // Compute local range
        double test_min = INFINITY, test_max = -INFINITY;
        for (unsigned int i = 0; i <= 20; ++i) {
            for (unsigned int j = 0; j <= 20; ++j) {
                double u = test_u_min + (test_u_max - test_u_min) * i / 20.0;
                double v = test_v_min + (test_v_max - test_v_min) * j / 20.0;
                double val = hessian_determinant(u, v);
                test_min = std::min(test_min, val);
                test_max = std::max(test_max, val);
            }
        }
        std::cout << "Local range: [" << test_min << ", " << test_max
                  << "], magnitude: " << (test_max - test_min) << "\n" << std::endl;

        // Method 1: Direct Bernstein approximation (current approach)
        std::cout << "Method 1: Direct Bernstein approximation (degree " << degree << "x" << degree << ")..." << std::endl;
        Polynomial test_poly_approx = interpolate_bernstein_approx(
            hessian_determinant, degree, degree, test_u_min, test_u_max, test_v_min, test_v_max);

        // Check interpolation error at many test points
        double max_abs_error_approx = 0.0;
        double sum_abs_error_approx = 0.0;
        int num_test_points = 0;

        for (unsigned int i = 0; i <= 9; ++i) {
            for (unsigned int j = 0; j <= 9; ++j) {
                // Use slightly offset points (not exactly on grid)
                double s = (i + 0.5) / 10.0;
                double t = (j + 0.5) / 10.0;
                double u = test_u_min + (test_u_max - test_u_min) * s;
                double v = test_v_min + (test_v_max - test_v_min) * t;

                // Evaluate both original and polynomial
                double orig_val = hessian_determinant(u, v);
                std::vector<double> local_coords = {s, t};
                double poly_val = test_poly_approx.evaluate(local_coords);

                double abs_error = std::abs(poly_val - orig_val);
                max_abs_error_approx = std::max(max_abs_error_approx, abs_error);
                sum_abs_error_approx += abs_error;
                num_test_points++;
            }
        }

        double avg_abs_error_approx = sum_abs_error_approx / num_test_points;
        double error_ratio_approx = max_abs_error_approx / (test_max - test_min);

        std::cout << "  Max absolute error: " << max_abs_error_approx << std::endl;
        std::cout << "  Avg absolute error: " << avg_abs_error_approx << std::endl;
        std::cout << "  Error/Range ratio: " << error_ratio_approx << std::endl;


        // Method 2: Proper polynomial interpolation
        std::cout << "\nMethod 2: Proper polynomial interpolation (degree " << degree << "x" << degree << ")..." << std::endl;
        Polynomial test_poly_proper = interpolate_proper(
            hessian_determinant, degree, degree, test_u_min, test_u_max, test_v_min, test_v_max);

        double max_abs_error_proper = 0.0;
        double sum_abs_error_proper = 0.0;
        num_test_points = 0;

        for (unsigned int i = 0; i <= 9; ++i) {
            for (unsigned int j = 0; j <= 9; ++j) {
                double s = (i + 0.5) / 10.0;
                double t = (j + 0.5) / 10.0;
                double u = test_u_min + (test_u_max - test_u_min) * s;
                double v = test_v_min + (test_v_max - test_v_min) * t;

                double orig_val = hessian_determinant(u, v);
                std::vector<double> local_coords = {s, t};
                double poly_val = test_poly_proper.evaluate(local_coords);

                double abs_error = std::abs(poly_val - orig_val);
                max_abs_error_proper = std::max(max_abs_error_proper, abs_error);
                sum_abs_error_proper += abs_error;
                num_test_points++;
            }
        }

        double avg_abs_error_proper = sum_abs_error_proper / num_test_points;
        double error_ratio_proper = max_abs_error_proper / (test_max - test_min);

        std::cout << "  Max absolute error: " << max_abs_error_proper << std::endl;
        std::cout << "  Avg absolute error: " << avg_abs_error_proper << std::endl;
        std::cout << "  Error/Range ratio: " << error_ratio_proper << std::endl;

        // Compare methods
        std::cout << "\nComparison:" << std::endl;
        std::cout << "  Approximation error ratio: " << error_ratio_approx << std::endl;
        std::cout << "  Proper interpolation error ratio: " << error_ratio_proper << std::endl;

        if (error_ratio_proper < error_ratio_approx * 0.5) {
            std::cout << "  ✅ Proper interpolation is significantly better!" << std::endl;
        } else if (error_ratio_proper < error_ratio_approx) {
            std::cout << "  ✓ Proper interpolation is somewhat better" << std::endl;
        } else {
            std::cout << "  ⚠️  Proper interpolation doesn't help much" << std::endl;
        }

        if (error_ratio_proper < 0.01) {
            std::cout << "  ✅ Proper interpolation quality is excellent (<1% error)\n" << std::endl;
        } else if (error_ratio_proper < 0.1) {
            std::cout << "  ✓ Proper interpolation quality is acceptable (<10% error)\n" << std::endl;
        } else {
            std::cout << "  ⚠️  Function is too complex for polynomial interpolation\n" << std::endl;
        }
    }

    // Verify Hessian computation first
    verify_hessian_computation();

    // Step 0: Analyze the Hessian determinant range (global)
    std::cout << "Step 0: Analyzing global Hessian determinant range..." << std::endl;
    double global_min = INFINITY, global_max = -INFINITY;
    for (unsigned int i = 0; i <= 100; ++i) {
        for (unsigned int j = 0; j <= 100; ++j) {
            double u = i / 100.0;
            double v = j / 100.0;
            double val = hessian_determinant(u, v);
            global_min = std::min(global_min, val);
            global_max = std::max(global_max, val);
        }
    }
    std::cout << "  Global Hessian det range: [" << global_min << ", " << global_max << "]" << std::endl;
    std::cout << "  Global range magnitude: " << (global_max - global_min) << "\n" << std::endl;

    // Step 1: Subdivide domain and solve in each region
    std::cout << "Step 1: Subdividing domain and solving in each region..." << std::endl;

    std::vector<SubdivisionBoxResult> all_converged;
    std::vector<SubdivisionBoxResult> all_unresolved;
    unsigned int total_resolved = 0;
    unsigned int total_unresolved = 0;
    unsigned int pruned_converged = 0;
    unsigned int pruned_unresolved = 0;
    bool any_degeneracy = false;

    double max_poly_residual = 0.0;
    double max_orig_residual = 0.0;

    // Residual pruning threshold for original Hessian determinant
    // Can be adjusted: 1e2 (strict), 1e3 (moderate), 1e4 (loose)
    const double residual_threshold = 1e2;  // Absolute threshold for |det(H)|

    for (unsigned int i = 0; i < subdivisions; ++i) {
        for (unsigned int j = 0; j < subdivisions; ++j) {
            // Define subdomain [u_min, u_max] x [v_min, v_max] in [0,1]^2
            double u_min = static_cast<double>(i) / subdivisions;
            double u_max = static_cast<double>(i + 1) / subdivisions;
            double v_min = static_cast<double>(j) / subdivisions;
            double v_max = static_cast<double>(j + 1) / subdivisions;

            std::cout << "\n  Region [" << i << "," << j << "]: u=[" << u_min << "," << u_max
                      << "], v=[" << v_min << "," << v_max << "]" << std::endl;

            // Analyze local range
            double local_min = INFINITY, local_max = -INFINITY;
            for (unsigned int ii = 0; ii <= 20; ++ii) {
                for (unsigned int jj = 0; jj <= 20; ++jj) {
                    double u = u_min + (u_max - u_min) * ii / 20.0;
                    double v = v_min + (v_max - v_min) * jj / 20.0;
                    double val = hessian_determinant(u, v);
                    local_min = std::min(local_min, val);
                    local_max = std::max(local_max, val);
                }
            }
            std::cout << "    Local range: [" << local_min << ", " << local_max
                      << "], magnitude: " << (local_max - local_min) << std::endl;

            // Use proper polynomial interpolation
            Polynomial local_poly = interpolate_proper(
                hessian_determinant, degree, degree, u_min, u_max, v_min, v_max);

            // Solve in this subdomain
            PolynomialSystem local_system({local_poly});

            SubdivisionConfig config = defaultSolverConfig();
            config.tolerance = 1e-4;
            config.max_depth = 12;
            config.degeneracy_multiplier = degeneracy_multiplier;

#ifdef ENABLE_GEOMETRY_DUMP
            // Geometry dumping is available in debug mode but disabled by default
            // config.dump_geometry = true;
            // config.dump_prefix = "dumps/hessian_det_region";
#endif

            Solver solver;
            SubdivisionSolverResult result = solver.subdivisionSolve(
                local_system, config, RootBoundingMethod::ProjectedPolyhedral);

            std::cout << "    Boxes: " << result.boxes.size()
                      << " (converged: " << result.num_resolved
                      << ", unresolved: " << (result.boxes.size() - result.num_resolved) << ")" << std::endl;

            if (result.degeneracy_detected) {
                any_degeneracy = true;
            }

            // Transform boxes back to global [0,1]^2 coordinates and compute residuals
            for (std::size_t k = 0; k < result.boxes.size(); ++k) {
                SubdivisionBoxResult box = result.boxes[k];

                // Transform from local [0,1]^2 to global [0,1]^2
                box.lower[0] = u_min + (u_max - u_min) * box.lower[0];
                box.lower[1] = v_min + (v_max - v_min) * box.lower[1];
                box.upper[0] = u_min + (u_max - u_min) * box.upper[0];
                box.upper[1] = v_min + (v_max - v_min) * box.upper[1];
                box.center[0] = u_min + (u_max - u_min) * box.center[0];
                box.center[1] = v_min + (v_max - v_min) * box.center[1];
                box.max_error[0] = (u_max - u_min) * box.max_error[0];
                box.max_error[1] = (v_max - v_min) * box.max_error[1];

                // Evaluate original Hessian determinant at box center
                double u = box.center[0];
                double v = box.center[1];
                double orig_val = hessian_determinant(u, v);

                if (k < result.num_resolved) {
                    // Converged box - check residual before accepting

                    // Evaluate local polynomial at local coordinates
                    std::vector<double> local_coords = {
                        (u - u_min) / (u_max - u_min),
                        (v - v_min) / (v_max - v_min)
                    };
                    double poly_val = local_poly.evaluate(local_coords);

                    max_poly_residual = std::max(max_poly_residual, std::abs(poly_val));
                    max_orig_residual = std::max(max_orig_residual, std::abs(orig_val));

                    // Prune based on original function residual
                    if (std::abs(orig_val) <= residual_threshold) {
                        all_converged.push_back(box);
                    } else {
                        pruned_converged++;
                    }
                } else {
                    // Unresolved box - also check if it's worth keeping
                    // Check residual at box center and corners
                    double min_residual = std::abs(orig_val);

                    // Check corners
                    double corners[4][2] = {
                        {box.lower[0], box.lower[1]},
                        {box.upper[0], box.lower[1]},
                        {box.lower[0], box.upper[1]},
                        {box.upper[0], box.upper[1]}
                    };

                    for (int c = 0; c < 4; ++c) {
                        double corner_val = hessian_determinant(corners[c][0], corners[c][1]);
                        min_residual = std::min(min_residual, std::abs(corner_val));
                    }

                    // Keep unresolved box if minimum residual is below threshold
                    if (min_residual <= residual_threshold) {
                        all_unresolved.push_back(box);
                    } else {
                        pruned_unresolved++;
                    }
                }
            }

            total_resolved += result.num_resolved;
            total_unresolved += (result.boxes.size() - result.num_resolved);
        }
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Overall Results (before pruning):" << std::endl;
    std::cout << "  Total converged boxes: " << total_resolved << std::endl;
    std::cout << "  Total unresolved boxes: " << total_unresolved << std::endl;
    std::cout << "  Degeneracy detected: " << (any_degeneracy ? "yes" : "no") << std::endl;

    std::cout << "\nResidual-based pruning (threshold = " << residual_threshold << "):" << std::endl;
    std::cout << "  Pruned converged boxes: " << pruned_converged
              << " (" << (100.0 * pruned_converged / std::max(1u, total_resolved)) << "%)" << std::endl;
    std::cout << "  Pruned unresolved boxes: " << pruned_unresolved
              << " (" << (100.0 * pruned_unresolved / std::max(1u, total_unresolved)) << "%)" << std::endl;

    std::cout << "\nFinal results (after pruning):" << std::endl;
    std::cout << "  Kept converged boxes: " << all_converged.size() << std::endl;
    std::cout << "  Kept unresolved boxes: " << all_unresolved.size() << std::endl;

    // Recompute residuals for kept boxes only
    if (!all_converged.empty()) {
        double max_kept_residual = 0.0;
        for (const auto& box : all_converged) {
            double u = box.center[0];
            double v = box.center[1];
            double orig_val = hessian_determinant(u, v);
            max_kept_residual = std::max(max_kept_residual, std::abs(orig_val));
        }

        std::cout << "\nKept converged box residuals:" << std::endl;
        std::cout << "  Max original function residual: " << max_kept_residual << std::endl;
        std::cout << "  (All kept boxes have residual <= " << residual_threshold << ")" << std::endl;
    }

    // Step 2: Write converged boxes (points on the curve)
    std::ofstream converged_file("dumps/hessian_det_converged.txt");
    converged_file << "# Hessian determinant zero set - converged boxes (points)\n";
    converged_file << "# Format: u_center v_center\n";
    converged_file << "# Note: (u,v) in [0,1]^2, corresponds to (x,y) = (2u-1, 2v-1) in [-1,1]^2\n";

    for (const auto& box : all_converged) {
        converged_file << box.center[0] << " " << box.center[1] << "\n";
    }
    converged_file.close();

    // Step 3: Write unresolved boxes to file for visualization
    std::ofstream boxes_file("dumps/hessian_det_boxes.txt");
    boxes_file << "# Hessian determinant zero set boxes (unresolved only)\n";
    boxes_file << "# Format: u_min u_max v_min v_max\n";
    boxes_file << "# Note: (u,v) in [0,1]^2, corresponds to (x,y) = (2u-1, 2v-1) in [-1,1]^2\n";

    for (const auto& box : all_unresolved) {
        boxes_file << box.lower[0] << " " << box.upper[0] << " "
                  << box.lower[1] << " " << box.upper[1] << "\n";
    }
    boxes_file.close();

    std::cout << "\n✅ Converged points written to: dumps/hessian_det_converged.txt" << std::endl;
    std::cout << "✅ Unresolved boxes written to: dumps/hessian_det_boxes.txt" << std::endl;

    // Step 4: Sample the Hessian determinant to understand the zero set
    std::cout << "\nStep 4: Sampling Hessian determinant for visualization..." << std::endl;

    std::ofstream sample_file("dumps/hessian_det_samples.txt");
    sample_file << "# Hessian determinant samples on [0,1]^2\n";
    sample_file << "# Format: u v det(H) x y\n";
    sample_file << "# where (x,y) = (2u-1, 2v-1)\n";

    const int n_samples = 100;
    for (int j = 0; j <= n_samples; ++j) {
        for (int i = 0; i <= n_samples; ++i) {
            double u = static_cast<double>(i) / n_samples;
            double v = static_cast<double>(j) / n_samples;
            double det_h = hessian_determinant(u, v);
            double x = 2.0 * u - 1.0;
            double y = 2.0 * v - 1.0;

            sample_file << u << " " << v << " " << det_h << " " << x << " " << y << "\n";
        }
    }
    sample_file.close();

    std::cout << "✅ Samples written to: dumps/hessian_det_samples.txt" << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary" << std::endl;
    std::cout << "========================================\n" << std::endl;

    std::cout << "The zero set of det(Hessian) separates regions where:" << std::endl;
    std::cout << "  • det(H) > 0: Locally convex or concave (check trace)" << std::endl;
    std::cout << "  • det(H) < 0: Saddle point (neither convex nor concave)" << std::endl;
    std::cout << "  • det(H) = 0: Boundary between regions\n" << std::endl;

    std::cout << "To visualize:" << std::endl;
    std::cout << "  python3 visualize_hessian_det.py\n" << std::endl;

    return 0;
}

