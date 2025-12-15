/**
 * @file multiplicity_1d_roots.cpp
 * @brief Example demonstrating high-multiplicity roots with automatic HP escalation
 *
 * This example shows the new workflow:
 * 1. Solver gives initial guess
 * 2. Double-precision refiner tries to refine and sets needs_higher_precision flag
 * 3. If flag is set, HP refiner takes over (detects multiplicity + refines)
 *
 * Problem: p(x) = (x - 0.2)(x - 0.6)^6
 * Expected roots: x = 0.2 (multiplicity 1), x = 0.6 (multiplicity 6)
 */

#include <polynomial_solver.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>

#ifdef ENABLE_HIGH_PRECISION
#include "polynomial_hp.h"
#include "result_refiner_hp.h"
#include "precision_context.h"
#endif

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

#ifdef ENABLE_HIGH_PRECISION
// Helper: Create (x - root)^mult directly in HP to avoid double->HP conversion errors
PolynomialHP createMultipleRootHP(double root_val, unsigned int mult) {
    mpreal root = mpreal(std::to_string(root_val));
    std::vector<mpreal> power_coeffs = {-root, mpreal("1.0")};

    for (unsigned int i = 1; i < mult; ++i) {
        std::vector<mpreal> new_coeffs(power_coeffs.size() + 1, mpreal("0.0"));
        for (size_t j = 0; j < power_coeffs.size(); ++j) {
            new_coeffs[j] += power_coeffs[j] * (-root);
            new_coeffs[j+1] += power_coeffs[j];
        }
        power_coeffs = new_coeffs;
    }

    std::vector<unsigned int> degrees = {mult};
    return fromPowerHP(degrees, power_coeffs);
}
#endif

void print_help() {
    std::cout << "Usage: multiplicity_1d_roots [OPTIONS]\n\n";
    std::cout << "Solve polynomial with high-multiplicity root: (x - 0.2)(x - 0.6)^6\n\n";
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

int main(int argc, char* argv[]) {
    // Parse command-line arguments
    Config config = parse_args(argc, argv);

    if (config.show_help) {
        print_help();
        return 0;
    }

    std::cout << "========================================\n";
    std::cout << "High-Multiplicity Roots: 2-Line Workflow\n";
    std::cout << "========================================\n\n";

    std::cout << "Problem: p(x) = (x - 0.2)(x - 0.6)^6\n";
    std::cout << "Expected roots: x = 0.2 (multiplicity 1), x = 0.6 (multiplicity 6)\n\n";

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

    // Define polynomial: (x - 0.2)(x - 0.6)^6
    // Expanded: x^7 - 3.8*x^6 + 6.12*x^5 - 5.4*x^4 + 2.808*x^3 - 0.85536*x^2 + 0.139968*x - 0.0093312
    std::vector<unsigned int> degrees = {7};
    std::vector<double> power_coeffs = {-0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0};
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
        solver_config.dump_prefix = "dumps/multiplicity_1d";
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
        dump_result(result, "dumps/multiplicity_1d_result.txt");
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

#ifdef ENABLE_HIGH_PRECISION
    // ============================================================
    // STEP 3: HP REFINEMENT (for roots that need it)
    // ============================================================
    std::vector<RefinedRoot> hp_refined_roots;

    // Create HP polynomial: (x - 0.2)(x - 0.6)^6 = (x - 1/5)(x - 3/5)^6
    // Use EXACT rational coefficients (computed symbolically, not from decimal approximations)
    // This is critical for achieving high-precision error bounds!
    PrecisionContext ctx_construct(1024);

    std::vector<unsigned int> hp_degrees = {7};
    std::vector<mpreal> hp_power_coeffs = {
        mpreal("-729") / mpreal("78125"),      // constant: -(1/5)*(3/5)^6 = -729/78125
        mpreal("2187") / mpreal("15625"),      // x: (3/5)^6 + 6*(1/5)*(3/5)^5 = 2187/15625
        mpreal("-13365") / mpreal("15625"),    // x^2: -6*(3/5)^5 - 15*(1/5)*(3/5)^4 = -13365/15625
        mpreal("8775") / mpreal("3125"),       // x^3: 15*(3/5)^4 + 20*(1/5)*(3/5)^3 = 8775/3125
        mpreal("-27") / mpreal("5"),           // x^4: -20*(3/5)^3 - 15*(1/5)*(3/5)^2 = -27/5
        mpreal("153") / mpreal("25"),          // x^5: 15*(3/5)^2 + 6*(1/5)*(3/5) = 153/25
        mpreal("-19") / mpreal("5"),           // x^6: -6*(3/5) - (1/5) = -19/5
        mpreal("1")                            // x^7: 1
    };
    PolynomialHP poly_hp_full = fromPowerHP(hp_degrees, hp_power_coeffs);

    // Process verified roots that need HP
    for (const auto& root : refined.roots) {
        if (root.needs_higher_precision) {
            std::cout << "Root at x=" << root.location[0]
                      << " needs HP (condition=" << root.condition_estimate << ")\n";

            // Set precision based on condition number
            unsigned int precision_bits = 512;
            if (root.condition_estimate > 1e15) precision_bits = 1024;

            PrecisionContext ctx(precision_bits);
            std::cout << "  Using " << precision_bits << " bits precision\n";

            // HP refiner detects multiplicity and refines
            RefinementConfigHP hp_config;
            hp_config.max_newton_iters = 50;
            hp_config.multiplicity_hint = 6;  // Hint: we know it's multiplicity 6
            hp_config.multiplicity_method = MultiplicityMethod::HINT;
            hp_config.iteration_method = IterationMethod::MODIFIED_NEWTON;
            hp_config.multiplicity_timing = MultiplicityTiming::ONCE_AT_START;

            // For numerical solver, target is double precision (1e-15), not arbitrary precision
            // With 1024-bit precision, degree-7 polynomial with multiplicity-6 root can achieve ~1e-33
            // This is more than sufficient for numerical purposes
            hp_config.target_tolerance_str = "1e-15";
            hp_config.residual_tolerance_str = "1e-50";

            RefinedRootHP hp_result = ResultRefinerHP::refineRoot1D(
                root.location[0], poly_hp_full, hp_config);

            if (hp_result.converged) {
                std::cout << "  ✅ HP refinement succeeded!\n";
                std::cout << "  Detected multiplicity: " << hp_result.multiplicity << "\n";
                std::cout << "  Refined x: " << hp_result.location << "\n";
                std::cout << "  Residual: " << hp_result.residual << "\n";
                std::cout << "  Error bound: " << hp_result.max_error << "\n";
                std::cout << "  Iterations: " << hp_result.iterations << "\n";

                // Create refined root entry
                RefinedRoot hp_root = root;
                hp_root.location[0] = static_cast<double>(hp_result.location);
                hp_root.residual[0] = static_cast<double>(hp_result.residual);
                hp_root.multiplicity = hp_result.multiplicity;
                hp_root.needs_higher_precision = false;
                hp_root.verified = true;
                hp_refined_roots.push_back(hp_root);
            } else {
                std::cout << "  ❌ HP refinement failed: " << hp_result.error_message << "\n";
                hp_refined_roots.push_back(root);  // Keep original
            }
            std::cout << "\n";
        } else {
            hp_refined_roots.push_back(root);  // Keep double-precision result
        }
    }

    // Process problematic regions with HP
    for (const auto& region : refined.problematic_regions) {
        std::cout << "Problematic region [" << region.lower << ", " << region.upper
                  << "] needs HP\n";

        // Use refined_root from double-precision attempt as initial guess (better than midpoint)
        double initial_guess = region.refinement_attempted ? region.refined_root :
                               (region.lower + region.upper) / 2.0;

        // Use high precision (1024 bits for problematic regions)
        PrecisionContext ctx(1024);
        std::cout << "  Using 1024 bits precision\n";

        // HP refiner detects multiplicity and refines
        RefinementConfigHP hp_config;
        hp_config.max_newton_iters = 50;
        hp_config.multiplicity_hint = 6;  // Hint: we know it's multiplicity 6
        hp_config.multiplicity_method = MultiplicityMethod::HINT;
        hp_config.iteration_method = IterationMethod::MODIFIED_NEWTON;
        hp_config.multiplicity_timing = MultiplicityTiming::ONCE_AT_START;

        // For numerical solver, target is double precision (1e-15), not arbitrary precision
        // With 1024-bit precision, degree-7 polynomial with multiplicity-6 root can achieve ~1e-33
        // This is more than sufficient for numerical purposes
        hp_config.target_tolerance_str = "1e-15";
        hp_config.residual_tolerance_str = "1e-50";

        RefinedRootHP hp_result = ResultRefinerHP::refineRoot1D(
            initial_guess, poly_hp_full, hp_config);

        if (hp_result.converged) {
            std::cout << "  ✅ HP refinement succeeded!\n";
            std::cout << "  Detected multiplicity: " << hp_result.multiplicity << "\n";
            std::cout << "  Refined x: " << hp_result.location << "\n";
            std::cout << "  Residual: " << hp_result.residual << "\n";
            std::cout << "  Error bound: " << hp_result.max_error << "\n";
            std::cout << "  Iterations: " << hp_result.iterations << "\n";

            // Create refined root entry
            RefinedRoot hp_root;
            hp_root.location.push_back(static_cast<double>(hp_result.location));
            hp_root.residual.push_back(static_cast<double>(hp_result.residual));
            hp_root.multiplicity = hp_result.multiplicity;
            hp_root.needs_higher_precision = false;
            hp_root.verified = true;
            hp_refined_roots.push_back(hp_root);
        } else {
            std::cout << "  ❌ HP refinement failed: " << hp_result.error_message << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "Step 3: HP Refinement\n";
    std::cout << "  Total roots verified: " << hp_refined_roots.size() << "\n\n";
#endif

    // ============================================================
    // RESULTS
    // ============================================================
    std::cout << "========================================\n";
    std::cout << "Results\n";
    std::cout << "========================================\n\n";

#ifdef ENABLE_HIGH_PRECISION
    // Display HP-refined roots if available
    const auto& final_roots = hp_refined_roots.empty() ? refined.roots : hp_refined_roots;
#else
    const auto& final_roots = refined.roots;
#endif

    // Display verified roots
    if (!final_roots.empty()) {
        std::cout << "Verified Roots:\n";
        for (std::size_t i = 0; i < final_roots.size(); ++i) {
            const auto& root = final_roots[i];
            std::cout << "  Root " << (i + 1) << ":\n";
            std::cout << "    Location: x = " << std::setprecision(16) << std::fixed
                      << root.location[0] << "\n";
            std::cout << "    Residual: |f(x)| = " << std::scientific << std::setprecision(4)
                      << std::abs(root.residual[0]) << "\n";
            std::cout << "    Multiplicity: " << root.multiplicity << "\n";
            std::cout << "    Condition estimate: " << std::scientific << std::setprecision(2)
                      << root.condition_estimate << "\n";

            if (root.needs_higher_precision) {
                std::cout << "    ⚠️  WARNING: Higher precision recommended!\n";
            } else {
                std::cout << "    ✅ Double precision sufficient\n";
            }
            std::cout << "\n";
        }
    }

    // Display problematic regions
    if (!refined.problematic_regions.empty()) {
        std::cout << "Problematic Regions (need higher precision):\n";
        for (std::size_t i = 0; i < refined.problematic_regions.size(); ++i) {
            const auto& region = refined.problematic_regions[i];
            std::cout << "  Region " << (i + 1) << ":\n";
            std::cout << "    Interval: [" << std::fixed << std::setprecision(10)
                      << region.lower << ", " << region.upper << "]\n";
            std::cout << "    Number of boxes merged: " << region.box_indices.size() << "\n";

            if (region.refinement_attempted) {
                std::cout << "    Refined root: x = " << std::setprecision(16)
                          << region.refined_root << "\n";
                std::cout << "    Residual: " << std::scientific << std::setprecision(4)
                          << std::abs(region.residual) << "\n";
                std::cout << "    Estimated multiplicity: " << region.multiplicity << "\n";
                std::cout << "    ⚠️  Refinement " << (region.refinement_succeeded ? "succeeded but" : "failed,")
                          << " needs higher precision\n";
            }
            std::cout << "\n";
        }
    }

    std::cout << "========================================\n";
    std::cout << "Summary\n";
    std::cout << "========================================\n\n";
    std::cout << "High-multiplicity roots are challenging in double precision.\n";
    std::cout << "The multiplicity-6 root at x=0.6 requires higher precision arithmetic.\n";
    std::cout << "See docs/CONDITIONING_AND_PRECISION.md for more information.\n\n";

    return 0;
}

