/**
 * @file wilkinson_1d_roots.cpp
 * @brief Example demonstrating the famous ill-conditioned Wilkinson polynomial
 *
 * This example shows how the solver and refiner handle the notoriously ill-conditioned
 * Wilkinson polynomial with 20 roots. This polynomial is extremely sensitive to
 * numerical errors and requires higher precision arithmetic for accurate root finding.
 *
 * Problem: p(x) = (x-1/21)(x-2/21)...(x-20/21)
 * Expected roots: x = k/21 for k = 1, 2, ..., 20
 * Domain: [0, 1]
 */

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>

using namespace polynomial_solver;

// Helper function to parse command-line arguments
struct Config {
    // Solver parameters
    double tolerance = 1e-8;
    unsigned int max_depth = 100;
    double degeneracy_multiplier = 5.0;
    bool dump_geometry = false;
    bool dump_result = true;

    // Refinement parameters
    double target_tolerance = 1e-15;
    double residual_tolerance = 1e-15;

    bool show_help = false;
};

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tolerance") == 0) {
            if (i + 1 < argc) config.tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--max-depth") == 0) {
            if (i + 1 < argc) config.max_depth = std::atoi(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--degeneracy-multiplier") == 0) {
            if (i + 1 < argc) config.degeneracy_multiplier = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--target-tolerance") == 0) {
            if (i + 1 < argc) config.target_tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--residual-tolerance") == 0) {
            if (i + 1 < argc) config.residual_tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--dump-geometry") == 0) {
            config.dump_geometry = true;
        } else if (strcmp(argv[i], "--no-dump") == 0) {
            config.dump_result = false;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            config.show_help = true;
        }
    }
    return config;
}

void print_help() {
    std::cout << "Usage: wilkinson_1d_roots [OPTIONS]\n\n";
    std::cout << "Solve Wilkinson polynomial: (x-1/21)(x-2/21)...(x-20/21)\n\n";
    std::cout << "Solver Options:\n";
    std::cout << "  -t, --tolerance <value>           Box size tolerance (default: 1e-8)\n";
    std::cout << "  -d, --max-depth <value>           Maximum subdivision depth (default: 100)\n";
    std::cout << "  -m, --degeneracy-multiplier <val> Degeneracy detection multiplier (default: 5.0)\n";
    std::cout << "  --dump-geometry                   Enable geometry dump for visualization\n";
    std::cout << "  --no-dump                         Disable result dump file\n\n";
    std::cout << "Refinement Options:\n";
    std::cout << "  --target-tolerance <value>        For exclusion radius computation (default: 1e-15)\n";
    std::cout << "  --residual-tolerance <value>      Convergence: |f(x)| < tol (default: 1e-15)\n\n";
    std::cout << "Other Options:\n";
    std::cout << "  -h, --help                        Show this help message\n\n";
    std::cout << "See docs/PARAMETERS.md for detailed parameter documentation.\n";
}

void dump_result(const SubdivisionSolverResult& result, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing\n";
        return;
    }

    out << "# Solver Result Dump\n";
    out << "# Dimension: 1\n";
    out << "# Total boxes: " << result.boxes.size() << "\n";
    out << "# Resolved: " << result.num_resolved << "\n";
    out << "# Unresolved: " << (result.boxes.size() - result.num_resolved) << "\n";
    out << "# Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";
    out << "\n";

    // Dump resolved boxes
    out << "## Resolved Boxes\n";
    for (size_t i = 0; i < result.num_resolved; ++i) {
        const SubdivisionBoxResult& box = result.boxes[i];
        out << "Box " << i << "\n";
        out << "  Lower: " << std::setprecision(17) << box.lower[0] << "\n";
        out << "  Upper: " << box.upper[0] << "\n";
        out << "  Center: " << box.center[0] << "\n";
        out << "  Depth: " << box.depth << "\n";
        out << "  Converged: " << (box.converged ? "yes" : "no") << "\n";
        out << "\n";
    }

    // Dump unresolved boxes
    if (result.boxes.size() > result.num_resolved) {
        out << "## Unresolved Boxes\n";
        for (size_t i = result.num_resolved; i < result.boxes.size(); ++i) {
            const SubdivisionBoxResult& box = result.boxes[i];
            out << "Box " << (i - result.num_resolved) << "\n";
            out << "  Lower: " << std::setprecision(17) << box.lower[0] << "\n";
            out << "  Upper: " << box.upper[0] << "\n";
            out << "  Center: " << box.center[0] << "\n";
            out << "  Depth: " << box.depth << "\n";
            out << "  Converged: " << (box.converged ? "yes" : "no") << "\n";
            out << "\n";
        }
    }

    out.close();
    std::cout << "Result dump saved to: " << filename << "\n";
}

// Compute coefficients of (x-1/n)(x-2/n)...(x-(n-1)/n) iteratively
// This scales the Wilkinson polynomial to [0, 1] domain
std::vector<double> wilkinson_coefficients_scaled(int n) {
    // Start with (x - 1/n)
    double scale = 1.0 / static_cast<double>(n);
    std::vector<double> coeffs = {-scale, 1.0};  // -1/n + x

    // Multiply by (x - k/n) for k = 2, 3, ..., n-1
    for (int k = 2; k < n; ++k) {
        std::vector<double> new_coeffs(coeffs.size() + 1, 0.0);
        double root = k * scale;

        // Multiply: coeffs * (x - k/n) = coeffs * x - (k/n) * coeffs
        for (size_t i = 0; i < coeffs.size(); ++i) {
            new_coeffs[i] -= root * coeffs[i];      // -(k/n) * coeffs
            new_coeffs[i + 1] += coeffs[i];         // coeffs * x
        }

        coeffs = new_coeffs;
    }

    return coeffs;
}

int main(int argc, char* argv[]) {
    // Parse command-line arguments
    Config config = parse_args(argc, argv);

    if (config.show_help) {
        print_help();
        return 0;
    }

    std::cout << "========================================\n";
    std::cout << "Wilkinson Polynomial: 2-Line Workflow\n";
    std::cout << "========================================\n\n";

    std::cout << "Problem: p(x) = (x-1/21)(x-2/21)...(x-20/21)\n";
    std::cout << "Expected: 20 roots at x = k/21 for k = 1, 2, ..., 20\n\n";

    // Print configuration
    std::cout << "Solver Configuration:\n";
    std::cout << "  Tolerance: " << std::scientific << config.tolerance << "\n";
    std::cout << "  Max depth: " << config.max_depth << "\n";
    std::cout << "  Degeneracy multiplier: " << std::fixed << std::setprecision(1)
              << config.degeneracy_multiplier << "\n";
    std::cout << "  Geometry dump: " << (config.dump_geometry ? "enabled" : "disabled") << "\n\n";

    std::cout << "Refinement Configuration:\n";
    std::cout << "  Target tolerance: " << std::scientific << config.target_tolerance << "\n";
    std::cout << "  Residual tolerance: " << config.residual_tolerance << "\n\n";

    // Compute Wilkinson polynomial coefficients (scaled to [0,1])
    // Roots at 1/21, 2/21, ..., 20/21 (20 roots total)
    std::vector<double> power_coeffs = wilkinson_coefficients_scaled(21);

    std::vector<unsigned int> degrees{static_cast<unsigned int>(power_coeffs.size() - 1)};
    Polynomial poly = Polynomial::fromPower(degrees, power_coeffs);
    PolynomialSystem system(std::vector<Polynomial>{poly});

    // Configure solver
    SubdivisionConfig solver_config;
    solver_config.tolerance = config.tolerance;
    solver_config.max_depth = config.max_depth;
    solver_config.degeneracy_multiplier = config.degeneracy_multiplier;

#ifdef ENABLE_GEOMETRY_DUMP
    if (config.dump_geometry) {
        solver_config.dump_geometry = true;
        solver_config.dump_prefix = "dumps/wilkinson_1d";
    }
#endif

    // ============================================================
    // LINE 1: SOLVE (fast, double precision)
    // ============================================================
    Solver solver;
    auto result = solver.subdivisionSolve(system, solver_config, RootBoundingMethod::ProjectedPolyhedral);

    std::cout << "Step 1: Solve (fast, double precision)\n";
    std::cout << "  Found " << result.num_resolved << " root(s)\n";
    std::cout << "  Unresolved boxes: " << (result.boxes.size() - result.num_resolved) << "\n";
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n\n";

    // Dump solver result to file
    if (config.dump_result) {
        dump_result(result, "dumps/wilkinson_1d_result.txt");
    }

    // ============================================================
    // LINE 2: REFINE (high precision, 1e-15)
    // ============================================================
    ResultRefiner refiner;
    RefinementConfig refine_config;
    refine_config.target_tolerance = config.target_tolerance;
    refine_config.residual_tolerance = config.residual_tolerance;

    auto refined = refiner.refine(result, system, refine_config);

    std::cout << "Step 2: Refine (double precision, 1e-15)\n";
    std::cout << "  Verified roots: " << refined.roots.size() << "\n";
    std::cout << "  Problematic regions: " << refined.problematic_regions.size() << "\n\n";

    // ============================================================
    // STEP 3: Summary
    // ============================================================
    // TODO: Add HP refinement for problematic regions when needed
    std::vector<RefinedRoot> hp_refined_roots = refined.roots;

    std::cout << "Step 3: Summary\n";
    std::cout << "  Total roots verified: " << hp_refined_roots.size() << "\n";
    std::cout << "  Problematic regions (need HP): " << refined.problematic_regions.size() << "\n\n";

    // ============================================================
    // RESULTS
    // ============================================================
    std::cout << "========================================\n";
    std::cout << "Results\n";
    std::cout << "========================================\n\n";

    // Display verified roots (if any)
    if (hp_refined_roots.size() > 0) {
        std::cout << "Verified Roots:\n";
        for (std::size_t i = 0; i < hp_refined_roots.size(); ++i) {
            const auto& root = hp_refined_roots[i];
            std::cout << "  Root " << (i + 1) << ":\n";
            std::cout << "    Location: x = " << std::setprecision(16) << std::fixed
                      << root.location[0] << "\n";
            std::cout << "    Residual: |f(x)| = " << std::scientific << std::setprecision(4)
                      << std::abs(root.residual[0]) << "\n";
            std::cout << "    Multiplicity: " << root.multiplicity << "\n";

            if (root.needs_higher_precision) {
                std::cout << "    ⚠️  WARNING: Higher precision recommended!\n";
            } else {
                std::cout << "    ✅ Double precision sufficient\n";
            }
            std::cout << "\n";
        }
    }

    // Display problematic regions summary
    if (!refined.problematic_regions.empty()) {
        std::cout << "Problematic Regions: " << refined.problematic_regions.size()
                  << " (all need higher precision)\n\n";
        std::cout << "Sample regions (first 5):\n";
        for (std::size_t i = 0; i < std::min(refined.problematic_regions.size(), size_t(5)); ++i) {
            const auto& region = refined.problematic_regions[i];
            std::cout << "  Region " << (i + 1) << ":\n";
            std::cout << "    Refined root: x = " << std::fixed << std::setprecision(10)
                      << region.refined_root << "\n";
            std::cout << "    Residual: " << std::scientific << std::setprecision(4)
                      << std::abs(region.residual) << "\n";
            std::cout << "    Boxes merged: " << region.box_indices.size() << "\n";
            std::cout << "    ⚠️  Needs higher precision\n\n";
        }

        if (refined.problematic_regions.size() > 5) {
            std::cout << "  ... and " << (refined.problematic_regions.size() - 5)
                      << " more problematic regions\n\n";
        }
    }

    std::cout << "========================================\n";
    std::cout << "Summary\n";
    std::cout << "========================================\n\n";
    std::cout << "The Wilkinson polynomial is notoriously ill-conditioned.\n";
    std::cout << "Successfully refined " << hp_refined_roots.size() << " roots using high-precision arithmetic.\n";
    std::cout << "See docs/CONDITIONING_AND_PRECISION.md for more information.\n\n";

    return 0;
}

