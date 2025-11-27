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

// Interpolate a function with Bernstein polynomial of degree k
Polynomial interpolate_bernstein(double (*func)(double, double), unsigned int degree_u, unsigned int degree_v) {
    std::cout << "Interpolating function with Bernstein polynomial of degree (" 
              << degree_u << ", " << degree_v << ")..." << std::endl;
    
    // Sample the function at Bernstein interpolation points (Chebyshev-like)
    std::vector<double> bernstein_coeffs;
    
    for (unsigned int j = 0; j <= degree_v; ++j) {
        for (unsigned int i = 0; i <= degree_u; ++i) {
            // Use uniform sampling for simplicity
            double u = static_cast<double>(i) / degree_u;
            double v = static_cast<double>(j) / degree_v;
            
            double value = func(u, v);
            bernstein_coeffs.push_back(value);
        }
    }
    
    std::vector<unsigned int> degrees = {degree_u, degree_v};
    return Polynomial(degrees, bernstein_coeffs);
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
    bool any_degeneracy = false;

    double max_poly_residual = 0.0;
    double max_orig_residual = 0.0;

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

            // Interpolate in this subdomain
            // We sample the function at uniform grid points and use them as Bernstein coefficients
            // This is an approximation, but works well for piecewise interpolation
            std::vector<unsigned int> degrees_local = {degree, degree};
            std::vector<double> coeffs_local;

            for (unsigned int jj = 0; jj <= degree; ++jj) {
                for (unsigned int ii = 0; ii <= degree; ++ii) {
                    // Sample at uniform grid in local [0,1]^2
                    double s = static_cast<double>(ii) / degree;
                    double t = static_cast<double>(jj) / degree;

                    // Transform to global [0,1]^2
                    double u = u_min + (u_max - u_min) * s;
                    double v = v_min + (v_max - v_min) * t;

                    // Evaluate Hessian determinant
                    double val = hessian_determinant(u, v);
                    coeffs_local.push_back(val);
                }
            }

            // Use Bernstein coefficients directly (approximation)
            Polynomial local_poly = Polynomial::fromBernstein(degrees_local, coeffs_local);

            // Solve in this subdomain
            PolynomialSystem local_system({local_poly});

            SubdivisionConfig config = defaultSolverConfig();
            config.tolerance = 1e-4;
            config.max_depth = 12;
            config.degeneracy_multiplier = degeneracy_multiplier;

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

                if (k < result.num_resolved) {
                    // Converged box - compute residuals
                    double u = box.center[0];
                    double v = box.center[1];

                    // Evaluate local polynomial at local coordinates
                    std::vector<double> local_coords = {
                        (u - u_min) / (u_max - u_min),
                        (v - v_min) / (v_max - v_min)
                    };
                    double poly_val = local_poly.evaluate(local_coords);
                    double orig_val = hessian_determinant(u, v);

                    max_poly_residual = std::max(max_poly_residual, std::abs(poly_val));
                    max_orig_residual = std::max(max_orig_residual, std::abs(orig_val));

                    all_converged.push_back(box);
                } else {
                    all_unresolved.push_back(box);
                }
            }

            total_resolved += result.num_resolved;
            total_unresolved += (result.boxes.size() - result.num_resolved);
        }
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Overall Results:" << std::endl;
    std::cout << "  Total converged boxes: " << total_resolved << std::endl;
    std::cout << "  Total unresolved boxes: " << total_unresolved << std::endl;
    std::cout << "  Degeneracy detected: " << (any_degeneracy ? "yes" : "no") << std::endl;

    if (total_resolved > 0) {
        std::cout << "\nConverged box residuals:" << std::endl;
        std::cout << "  Max polynomial residual: " << max_poly_residual << std::endl;
        std::cout << "  Max original function residual: " << max_orig_residual << std::endl;
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

