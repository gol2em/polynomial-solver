/**
 * @file hessian_zero_set.cpp
 * @brief Find the zero set of the Hessian determinant of a nonlinear function
 *
 * WORKFLOW:
 *   1. Divide domain into subregions for better local approximation
 *   2. In each subregion, interpolate f(x,y) as a polynomial
 *   3. Compute symbolic Hessian using Differentiation::hessian()
 *   4. Compute det(H) using polynomial arithmetic
 *   5. Find zero set using subdivision solver
 *   6. Refine box centers onto the curve using Newton's method with numerical gradient
 *
 * USAGE: ./hessian_zero_set [options]
 *   -r <half_width>   Region = [-r, r]², default: 1.5
 *   -s <subdivisions> Subdivisions per axis, default: 4
 *   -d <degree>       Polynomial degree for interpolation, default: 10
 *   -n <boxes>        Target boxes per subregion, default: 2000
 *   -t <tolerance>    Solver box tolerance, default: 1e-6
 *   -hp               Use high-precision refinement
 *   -b <bits>         Precision bits for HP mode, default: 128 (~38 digits)
 *   -q                Quiet mode (machine-readable output)
 *
 * EXAMPLE: f(x,y) = exp(-(x² + y²)), zero set is circle r = 1/√2
 */

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstring>

#ifdef ENABLE_HIGH_PRECISION
#include <result_refiner_hp.h>
#include <precision_context.h>
#endif

using namespace polynomial_solver;

// ============================================================================
// USER-DEFINED FUNCTION: Modify this section for your application
// ============================================================================

// The function f(x,y) whose Hessian zero set we want to find
double f_user(double x, double y) {
    return std::exp(-(x * x + y * y));
}

// Expected result for verification (set to NaN if unknown)
double expected_radius() { return 1.0 / std::sqrt(2.0); }

// ============================================================================
// CORE ALGORITHM
// ============================================================================

struct Config {
    double half_width = 1.5;
    unsigned int subdivisions = 4;
    unsigned int degree = 10;
    unsigned int target_boxes = 2000;  // Target boxes per subregion
    double tolerance = 1e-6;
    bool quiet = false;
    bool high_precision = false;
    unsigned int precision_bits = 128;  // ~38 decimal digits
    std::string output_file;            // Output file for refined points (empty = no output)
};

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-r") == 0 && i+1 < argc) cfg.half_width = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "-s") == 0 && i+1 < argc) cfg.subdivisions = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "-d") == 0 && i+1 < argc) cfg.degree = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "-n") == 0 && i+1 < argc) cfg.target_boxes = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "-t") == 0 && i+1 < argc) cfg.tolerance = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "-hp") == 0) cfg.high_precision = true;
        else if (std::strcmp(argv[i], "-b") == 0 && i+1 < argc) cfg.precision_bits = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "-o") == 0 && i+1 < argc) cfg.output_file = argv[++i];
        else if (std::strcmp(argv[i], "-q") == 0) cfg.quiet = true;
        else if (std::strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: " << argv[0] << " [-r hw] [-s sub] [-d deg] [-n boxes] [-t tol] [-hp] [-b bits] [-o file] [-q]\n"
                      << "  -o file  Output refined points to file (format: x y per line)\n";
            std::exit(0);
        }
    }
    return cfg;
}

/// Compute Hessian det polynomial on subregion [u0,u1]×[v0,v1] of unit domain
Polynomial compute_hessian_det(const Domain2D& domain, double u0, double u1, double v0, double v1, unsigned int deg) {
    // Create subdomain that maps [0,1]² to [u0,u1]×[v0,v1] in unit domain, then to user domain
    auto local_f = [&domain, u0, u1, v0, v1](double s, double t) {
        double u = u0 + (u1 - u0) * s;
        double v = v0 + (v1 - v0) * t;
        double x, y;
        domain.fromUnit(u, v, x, y);
        return f_user(x, y);
    };
    Polynomial f = Interpolation::interpolate2D(local_f, deg, deg, 0.0, 1.0, 0.0, 1.0, AbscissaeType::CHEBYSHEV);
    auto H = Differentiation::hessian(f);
    return H[0][0] * H[1][1] - H[0][1] * H[0][1];
}

// Compute Hessian determinant numerically from f_user
// det(H) = f_xx * f_yy - f_xy^2
double hessian_det_numerical(double x, double y) {
    const double h = 1e-5;
    double f00 = f_user(x, y);
    double f_xx = (f_user(x+h, y) - 2*f00 + f_user(x-h, y)) / (h*h);
    double f_yy = (f_user(x, y+h) - 2*f00 + f_user(x, y-h)) / (h*h);
    double f_xy = (f_user(x+h, y+h) - f_user(x+h, y-h) - f_user(x-h, y+h) + f_user(x-h, y-h)) / (4*h*h);
    return f_xx * f_yy - f_xy * f_xy;
}

struct Region { double u0, u1, v0, v1; Polynomial det_H; };

int main(int argc, char* argv[]) {
    Config cfg = parse_args(argc, argv);
    const unsigned int nsub = cfg.subdivisions;
    const double exp_r = expected_radius();

    // Create domain for coordinate transformations
    Domain2D domain = Domain2D::symmetric(cfg.half_width);

    // Compute degeneracy_multiplier from target_boxes
    // The Hessian det polynomial has degree 2*(d-2) in each variable
    // expected_max_roots = (2*(d-2))^2 for a single curve
    // degeneracy_threshold = degeneracy_multiplier * expected_max_roots
    // So: degeneracy_multiplier = target_boxes / expected_max_roots
    unsigned int hess_deg = 2 * (cfg.degree - 2);
    unsigned int expected_max_roots = hess_deg * hess_deg;
    double degeneracy_multiplier = static_cast<double>(cfg.target_boxes) / expected_max_roots;

    if (!cfg.quiet) {
        std::cout << "Hessian Zero Set Finder\n"
                  << "=======================\n"
                  << "Region: [" << domain.x_min << ", " << domain.x_max << "]^2\n"
                  << "Subdivisions: " << nsub << "x" << nsub << "\n"
                  << "Polynomial degree: " << cfg.degree << " (Hessian det degree: " << hess_deg << ")\n"
                  << "Target boxes/subregion: " << cfg.target_boxes << "\n"
                  << "Solver tolerance: " << cfg.tolerance << "\n"
                  << "Expected radius: " << exp_r << "\n\n";
    }

    // Step 1-5: Solve in each subregion
    std::vector<SubdivisionBoxResult> all_boxes;
    std::vector<Region> regions;
    std::vector<std::size_t> box_region;

    for (unsigned int i = 0; i < nsub; ++i) {
        for (unsigned int j = 0; j < nsub; ++j) {
            double u_min = static_cast<double>(i) / nsub;
            double u_max = static_cast<double>(i + 1) / nsub;
            double v_min = static_cast<double>(j) / nsub;
            double v_max = static_cast<double>(j + 1) / nsub;

            Polynomial det_H = compute_hessian_det(domain, u_min, u_max, v_min, v_max, cfg.degree);
            std::size_t ridx = regions.size();
            regions.push_back({u_min, u_max, v_min, v_max, det_H});

            SubdivisionConfig scfg = defaultSolverConfig();
            scfg.tolerance = cfg.tolerance;
            scfg.max_depth = 100;  // High max_depth, rely on degeneracy_multiplier for stopping
            scfg.degeneracy_multiplier = degeneracy_multiplier;

            Solver solver;
            auto result = solver.subdivisionSolve(PolynomialSystem({det_H}), scfg, RootBoundingMethod::ProjectedPolyhedral);

            for (auto& box : result.boxes) {
                all_boxes.push_back(box);
                box_region.push_back(ridx);
            }

            if (!cfg.quiet) {
                std::cout << "Region [" << i << "," << j << "]: " << result.boxes.size() << " boxes\n";
            }
        }
    }

    // Step 6: Refine using numerical Hessian determinant (function evaluations only)
    double max_box_err = 0.0, max_refined_err = 0.0;
    unsigned int n_refined = 0;
    std::vector<std::pair<double, double>> refined_points;  // Store refined (x, y) points

#ifdef ENABLE_HIGH_PRECISION
    if (cfg.high_precision) {
        // High-precision refinement with auto-computed optimal parameters
        PrecisionContext ctx(cfg.precision_bits);
        auto refine_cfg_hp = CurveRefinementConfigHP::fromPrecisionBits(cfg.precision_bits);

        if (!cfg.quiet) {
            std::cout << "\nTotal boxes: " << all_boxes.size() << "\n";
            std::cout << "Refinement: HIGH PRECISION (" << cfg.precision_bits << " bits, "
                      << "h=" << refine_cfg_hp.step_size_str << ", "
                      << "tol=" << refine_cfg_hp.residual_tolerance_str << ")\n";
        }

        // High-precision f(x,y) and its Hessian determinant
        auto f_hp = [](const mpreal& x, const mpreal& y) -> mpreal {
            return exp(-(x * x + y * y));
        };
        auto hessian_det_hp = makeHessianDetFunctionHP(f_hp, refine_cfg_hp.step_size_str);

        mpreal exp_r_hp = mpreal(1) / sqrt(mpreal(2));

        for (std::size_t k = 0; k < all_boxes.size(); ++k) {
            const auto& box = all_boxes[k];
            const auto& reg = regions[box_region[k]];

            // Map box center: local [0,1]² → subregion → global unit → user domain
            double s = (box.lower[0] + box.upper[0]) / 2.0;
            double t = (box.lower[1] + box.upper[1]) / 2.0;
            double u = reg.u0 + (reg.u1 - reg.u0) * s;
            double v = reg.v0 + (reg.v1 - reg.v0) * t;
            double xc, yc;
            domain.fromUnit(u, v, xc, yc);

            double r = std::sqrt(xc * xc + yc * yc);
            max_box_err = std::max(max_box_err, std::abs(r - exp_r));

            auto result = refineCurveNumericalHP(hessian_det_hp, xc, yc, refine_cfg_hp);
            if (result.converged) {
                double rx = static_cast<double>(result.x);
                double ry = static_cast<double>(result.y);
                refined_points.emplace_back(rx, ry);

                mpreal rr = sqrt(result.x * result.x + result.y * result.y);
                double err = static_cast<double>(abs(rr - exp_r_hp));
                max_refined_err = std::max(max_refined_err, err);
                n_refined++;
            }
        }
    } else
#endif
    {
        // Double-precision refinement
        if (!cfg.quiet) {
            std::cout << "\nTotal boxes: " << all_boxes.size() << "\n";
            std::cout << "Refinement: double precision (h=1e-5, tol=1e-5)\n";
        }

        // Note: Numerical second derivatives have ~h^2 error where h=1e-5, so ~1e-6 residual
        CurveRefinementConfig refine_cfg;
        refine_cfg.residual_tolerance = 1e-5;  // Limited by numerical derivative accuracy
        refine_cfg.max_iterations = 50;

        for (std::size_t k = 0; k < all_boxes.size(); ++k) {
            const auto& box = all_boxes[k];
            const auto& reg = regions[box_region[k]];

            // Map box center: local [0,1]² → subregion → global unit → user domain
            double s = (box.lower[0] + box.upper[0]) / 2.0;
            double t = (box.lower[1] + box.upper[1]) / 2.0;
            double u = reg.u0 + (reg.u1 - reg.u0) * s;
            double v = reg.v0 + (reg.v1 - reg.v0) * t;
            double xc, yc;
            domain.fromUnit(u, v, xc, yc);

            double r = std::sqrt(xc * xc + yc * yc);
            max_box_err = std::max(max_box_err, std::abs(r - exp_r));

            auto result = refineCurveNumerical(hessian_det_numerical, xc, yc, refine_cfg);
            if (result.converged) {
                refined_points.emplace_back(result.x, result.y);
                double rr = std::sqrt(result.x * result.x + result.y * result.y);
                max_refined_err = std::max(max_refined_err, std::abs(rr - exp_r));
                n_refined++;
            }
        }
    }

    // Write refined points to file if requested
    if (!cfg.output_file.empty()) {
        std::ofstream ofs(cfg.output_file);
        if (ofs) {
            ofs << std::setprecision(17);
            for (const auto& pt : refined_points) {
                ofs << pt.first << " " << pt.second << "\n";
            }
            if (!cfg.quiet) {
                std::cout << "Wrote " << refined_points.size() << " points to " << cfg.output_file << "\n";
            }
        } else {
            std::cerr << "Error: Could not open " << cfg.output_file << " for writing\n";
        }
    }

    // Output summary
    if (cfg.quiet) {
        std::cout << all_boxes.size() << " " << max_box_err << " " << max_refined_err << "\n";
    } else {
        std::cout << "\n=== Results ===\n"
                  << "Refined: " << n_refined << "/" << all_boxes.size() << "\n"
                  << "Max box error:     " << std::scientific << max_box_err << "\n"
                  << "Max refined error: " << max_refined_err << "\n";
    }

    return 0;
}
