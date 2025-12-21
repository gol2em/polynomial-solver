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
 *   6. Refine box centers onto the curve using Newton's method
 *
 * USAGE: ./hessian_zero_set [options]
 *   -r <half_width>   Region = [-r, r]², default: 1.5
 *   -s <subdivisions> Subdivisions per axis, default: 4
 *   -d <degree>       Polynomial degree, default: 10
 *   -t <tolerance>    Solver tolerance, default: 1e-6
 *   -m <max_depth>    Max subdivision depth, default: 15
 *   -q                Quiet mode (machine-readable output)
 *
 * EXAMPLE: f(x,y) = exp(-(x² + y²)), zero set is circle r = 1/√2
 */

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstring>

using namespace polynomial_solver;

// ============================================================================
// USER-DEFINED FUNCTION: Modify this section for your application
// ============================================================================

// The function f(x,y) whose Hessian zero set we want to find
double f_user(double x, double y) {
    return std::exp(-(x * x + y * y));
}

// The Hessian determinant (for refinement). If unknown, use polynomial refinement.
double hessian_det_user(double x, double y) {
    double r2 = x * x + y * y;
    return 4.0 * (1.0 - 2.0 * r2) * std::exp(-2.0 * r2);
}

// Expected result for verification
double expected_radius() { return 1.0 / std::sqrt(2.0); }

// ============================================================================
// CORE ALGORITHM
// ============================================================================

struct Config {
    double half_width = 1.5;
    unsigned int subdivisions = 4;
    unsigned int degree = 10;
    double tolerance = 1e-6;
    unsigned int max_depth = 15;
    bool quiet = false;
};

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-r") == 0 && i+1 < argc) cfg.half_width = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "-s") == 0 && i+1 < argc) cfg.subdivisions = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "-d") == 0 && i+1 < argc) cfg.degree = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "-t") == 0 && i+1 < argc) cfg.tolerance = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "-m") == 0 && i+1 < argc) cfg.max_depth = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "-q") == 0) cfg.quiet = true;
        else if (std::strcmp(argv[i], "-h") == 0) {
            std::cerr << "Usage: " << argv[0] << " [-r hw] [-s sub] [-d deg] [-t tol] [-m depth] [-q]\n";
            std::exit(0);
        }
    }
    return cfg;
}

inline void to_original(double u, double v, double hw, double& x, double& y) {
    x = (2.0 * u - 1.0) * hw;
    y = (2.0 * v - 1.0) * hw;
}

Polynomial compute_hessian_det(double u0, double u1, double v0, double v1, double hw, unsigned int deg) {
    auto local_f = [=](double s, double t) {
        double u = u0 + (u1 - u0) * s, v = v0 + (v1 - v0) * t, x, y;
        to_original(u, v, hw, x, y);
        return f_user(x, y);
    };
    Polynomial f = Interpolation::interpolate2D(local_f, deg, deg, 0.0, 1.0, 0.0, 1.0, AbscissaeType::CHEBYSHEV);
    auto H = Differentiation::hessian(f);
    return H[0][0] * H[1][1] - H[0][1] * H[0][1];
}

bool refine_numerical(const std::function<double(double,double)>& g, double x0, double y0,
                      double tol, unsigned int max_iters, double& xo, double& yo) {
    double x = x0, y = y0;
    const double h = 1e-8;
    for (unsigned int i = 0; i < max_iters; ++i) {
        double val = g(x, y);
        if (std::abs(val) < tol) { xo = x; yo = y; return true; }
        double gx = (g(x+h, y) - g(x-h, y)) / (2*h);
        double gy = (g(x, y+h) - g(x, y-h)) / (2*h);
        double gs = gx*gx + gy*gy;
        if (gs < 1e-30) { xo = x; yo = y; return false; }
        double t = -val / gs;
        x += t * gx; y += t * gy;
    }
    xo = x; yo = y;
    return false;
}

struct Region { double u0, u1, v0, v1; Polynomial det_H; };

int main(int argc, char* argv[]) {
    Config cfg = parse_args(argc, argv);
    const double hw = cfg.half_width;
    const unsigned int nsub = cfg.subdivisions;
    const double exp_r = expected_radius();

    if (!cfg.quiet) {
        std::cout << "Hessian Zero Set Finder\n"
                  << "=======================\n"
                  << "Region: [-" << hw << ", " << hw << "]^2\n"
                  << "Subdivisions: " << nsub << "x" << nsub << "\n"
                  << "Polynomial degree: " << cfg.degree << "\n"
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

            Polynomial det_H = compute_hessian_det(u_min, u_max, v_min, v_max, hw, cfg.degree);
            std::size_t ridx = regions.size();
            regions.push_back({u_min, u_max, v_min, v_max, det_H});

            SubdivisionConfig scfg = defaultSolverConfig();
            scfg.tolerance = cfg.tolerance;
            scfg.max_depth = cfg.max_depth;
            scfg.degeneracy_multiplier = 5.0;

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

    if (!cfg.quiet) {
        std::cout << "\nTotal boxes: " << all_boxes.size() << "\n";
    }

    // Step 6: Refine using original function
    double max_box_err = 0.0, max_refined_err = 0.0;
    unsigned int n_refined = 0;

    for (std::size_t k = 0; k < all_boxes.size(); ++k) {
        const auto& box = all_boxes[k];
        const auto& reg = regions[box_region[k]];

        double s = (box.lower[0] + box.upper[0]) / 2.0;
        double t = (box.lower[1] + box.upper[1]) / 2.0;
        double u = reg.u0 + (reg.u1 - reg.u0) * s;
        double v = reg.v0 + (reg.v1 - reg.v0) * t;
        double xc, yc;
        to_original(u, v, hw, xc, yc);

        double r = std::sqrt(xc * xc + yc * yc);
        max_box_err = std::max(max_box_err, std::abs(r - exp_r));

        double xr, yr;
        if (refine_numerical(hessian_det_user, xc, yc, 1e-14, 50, xr, yr)) {
            double rr = std::sqrt(xr * xr + yr * yr);
            max_refined_err = std::max(max_refined_err, std::abs(rr - exp_r));
            n_refined++;
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
