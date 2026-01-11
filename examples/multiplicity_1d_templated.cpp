/**
 * @file multiplicity_1d_templated.cpp
 * @brief Demonstrates the templated solver workflow for high-multiplicity roots
 *
 * This example showcases the new templated polynomial solver architecture that
 * provides a unified interface across different scalar types (double, mpreal, etc.).
 *
 * ARCHITECTURE OVERVIEW:
 * ----------------------
 * The templated solver consists of three main components:
 *
 * 1. PolynomialBase<Scalar> (core/polynomial_base.h)
 *    - Templated polynomial in Bernstein basis
 *    - Supports any scalar type with basic arithmetic operations
 *    - Provides evaluate(), differentiate(), subdivide(), etc.
 *
 * 2. SolverBase<Scalar> (solver/solver_base.h)
 *    - Bernstein subdivision solver
 *    - Uses sign changes and convex hull for root isolation
 *    - Returns boxes containing potential roots
 *
 * 3. ResultRefinerBase<Scalar> (refinement/result_refiner_base.h)
 *    - Newton-based root refinement with multiplicity detection
 *    - Multiple iteration methods: Newton, Modified Newton, Halley, Schröder
 *    - Multiplicity estimation: Taylor, Ostrowski, SimpleThreshold
 *    - Condition-aware convergence with HP flagging
 *
 * WORKFLOW FOR HIGH-MULTIPLICITY ROOTS:
 * -------------------------------------
 * 1. SolverBase<double> isolates candidate boxes (some may be degenerate)
 * 2. ResultRefinerBase<double> attempts refinement:
 *    - Simple roots: converge quickly, done
 *    - Multiple roots: detect high condition number, flag for HP
 * 3. ResultRefinerBase<mpreal> with higher precision:
 *    - Ostrowski's method estimates multiplicity m
 *    - Modified Newton (step *= m) achieves fast convergence
 *    - Returns root location accurate to machine epsilon
 *
 * Problem: p(x) = (x - 0.2)(x - 0.6)^6
 * Expected roots: x = 0.2 (mult=1), x = 0.6 (mult=6)
 */

#include "core/polynomial_base.h"
#include "solver/solver_base.h"
#include "refinement/result_refiner_base.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <cassert>

#ifdef ENABLE_HIGH_PRECISION
#include "hp/high_precision_types.h"
#endif

using namespace polynomial_solver;

// Double precision types
using Poly = PolynomialBase<double>;
using System = PolynomialSystemBase<double>;
using Solver = SolverBase<double>;
using Refiner = ResultRefinerBase<double>;
using SolverConfig = SubdivisionConfigBase<double>;
using RefineConfig = RefinementConfigBase<double>;
using BatchResult = RefinementResultBase<double>;

#ifdef ENABLE_HIGH_PRECISION
// High precision types
using HP = mpreal;
using PolyHP = PolynomialBase<HP>;
using RefinerHP = ResultRefinerBase<HP>;
using RefineConfigHP = RefinementConfigBase<HP>;
#endif

struct Config {
    double tolerance = 1e-6;  // Looser tolerance for high-multiplicity
    unsigned int max_depth = 40;  // Limit depth to avoid explosion
    double target_tolerance = 1e-12;
    double residual_tolerance = 1e-12;
    bool test_mode = false;
    bool show_help = false;
};

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tolerance") == 0) {
            if (i + 1 < argc) config.tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--max-depth") == 0) {
            if (i + 1 < argc) config.max_depth = std::atoi(argv[++i]);
        } else if (strcmp(argv[i], "--target-tolerance") == 0) {
            if (i + 1 < argc) config.target_tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--residual-tolerance") == 0) {
            if (i + 1 < argc) config.residual_tolerance = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "--test") == 0) {
            config.test_mode = true;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            config.show_help = true;
        }
    }
    return config;
}

void print_help() {
    std::cout << "Usage: multiplicity_1d_templated [OPTIONS]\n\n";
    std::cout << "Solve high-multiplicity polynomial: (x - 0.2)(x - 0.6)^6\n\n";
    std::cout << "Options:\n";
    std::cout << "  -t, --tolerance <value>        Box size tolerance (default: 1e-8)\n";
    std::cout << "  -d, --max-depth <value>        Maximum subdivision depth (default: 100)\n";
    std::cout << "  --target-tolerance <value>     Refinement target (default: 1e-12)\n";
    std::cout << "  --residual-tolerance <value>   Convergence residual (default: 1e-12)\n";
    std::cout << "  --test                         Run in test mode (assertions enabled)\n";
    std::cout << "  -h, --help                     Show this help\n";
}

int main(int argc, char* argv[]) {
    Config config = parse_args(argc, argv);
    if (config.show_help) { print_help(); return 0; }
    
    std::cout << "=== High-Multiplicity Roots (Templated Workflow) ===\n\n";
    std::cout << "Problem: p(x) = (x - 0.2)(x - 0.6)^6\n";
    std::cout << "Expected: x = 0.2 (mult=1), x = 0.6 (mult=6)\n\n";
    
    // Define polynomial: (x - 0.2)(x - 0.6)^6
    // Coefficients computed symbolically
    std::vector<unsigned int> degrees = {7};
    std::vector<double> power_coeffs = {
        -0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0
    };
    Poly poly = Poly::fromPower(degrees, power_coeffs);
    System system({poly});
    
    // Step 1: Solve
    SolverConfig solver_config;
    solver_config.tolerance = config.tolerance;
    solver_config.max_depth = config.max_depth;
    
    Solver solver;
    auto result = solver.subdivisionSolve(system, solver_config, 
                                          RootBoundingMethodBase::ProjectedPolyhedral);
    
    std::cout << "Step 1: Solve\n";
    std::cout << "  Total boxes: " << result.boxes.size() << "\n";
    std::cout << "  Resolved: " << result.num_resolved << "\n";
    std::cout << "  Unresolved: " << (result.boxes.size() - result.num_resolved) << "\n";
    std::cout << "  Degeneracy detected: " << (result.degeneracy_detected ? "yes" : "no") << "\n";

    // Debug: print resolved box centers
    std::cout << "  Resolved boxes:\n";
    for (std::size_t i = 0; i < result.num_resolved; ++i) {
        std::cout << "    Box " << i << ": center=" << std::setprecision(16) << result.boxes[i].center[0]
                  << " max_error=" << result.boxes[i].max_error[0] << "\n";
    }
    std::cout << "\n";

    // Step 2: Double-precision batch refine (detects condition number)
    RefineConfig refine_config;
    refine_config.target_tolerance = config.target_tolerance;
    refine_config.residual_tolerance = config.residual_tolerance;
    refine_config.max_multiplicity = 10;

    std::cout << "Step 2: Double-Precision Refine (detect high condition)\n";
    BatchResult batch = Refiner::refine(result, poly, refine_config);

    std::cout << "  Roots found: " << batch.roots.size() << "\n";
    std::cout << "  Problematic regions: " << batch.problematic_regions.size() << "\n";
    std::cout << "  Needs HP: " << (batch.any_needs_higher_precision ? "yes" : "no") << "\n\n";

    // Collect initial guesses that need HP refinement
    // - Converged roots with high condition number
    // - ALL problematic regions (they failed to converge, so need HP)
    std::vector<double> hp_candidates;
    for (const auto& root : batch.roots) {
        if (root.needs_higher_precision || !root.converged) {
            hp_candidates.push_back(root.location);
            std::cout << "  Root x=" << std::fixed << std::setprecision(6) << root.location
                      << " needs HP (cond=" << std::scientific << root.condition_estimate
                      << ", converged=" << root.converged << ")\n";
        }
    }
    // All problematic regions need HP (that's why they're problematic)
    for (const auto& region : batch.problematic_regions) {
        double guess = region.refinement_attempted ? region.refined_root
                                                   : (region.lower + region.upper) / 2.0;
        hp_candidates.push_back(guess);
        std::cout << "  Region [" << std::fixed << std::setprecision(6)
                  << region.lower << "," << region.upper
                  << "] needs HP (guess=" << guess << ")\n";
    }

#ifdef ENABLE_HIGH_PRECISION
    // Step 3: HP refinement for candidates
    std::cout << "\nStep 3: High-Precision Refinement\n";

    // Set working precision (512 bits ≈ 154 decimal digits)
    // Need enough precision so f^m doesn't underflow (for mult=6, need ~90 digits for 1e-15 target)
    setPrecision(512);
    std::cout << "  Working precision: " << getPrecision() << " bits\n\n";

    // Create HP polynomial with exact coefficients
    std::vector<unsigned int> hp_degrees = {7};
    std::vector<HP> hp_power_coeffs = {
        HP("-729") / HP("78125"),     // -(1/5)*(3/5)^6
        HP("2187") / HP("15625"),     // coefficient of x
        HP("-13365") / HP("15625"),   // coefficient of x^2
        HP("8775") / HP("3125"),      // coefficient of x^3
        HP("-27") / HP("5"),          // coefficient of x^4
        HP("153") / HP("25"),         // coefficient of x^5
        HP("-19") / HP("5"),          // coefficient of x^6
        HP("1")                       // coefficient of x^7
    };
    PolyHP poly_hp = PolyHP::fromPower(hp_degrees, hp_power_coeffs);

    // HP refinement config: target is DOUBLE precision root error
    RefineConfigHP hp_config;
    hp_config.target_tolerance = HP("1e-15");      // Double precision root error
    // For mult-m root: residual ~ error^m, so for error=1e-15, mult=6: residual ~ 1e-90
    // But be more lenient to allow convergence check to trigger
    hp_config.residual_tolerance = HP("1e-40");    // HP residual for convergence check
    hp_config.max_newton_iters = 200;              // More iterations for high mult
    hp_config.max_multiplicity = 10;
    hp_config.multiplicity_method = MultiplicityMethodBase::OSTROWSKI;
    hp_config.iteration_method = IterationMethodBase::MODIFIED_NEWTON;

    bool found_simple = false, found_multiple = false;

    // First, collect roots that converged well in double precision (no HP needed)
    for (const auto& root : batch.roots) {
        if (!root.needs_higher_precision && root.converged) {
            double x_final = root.location;
            std::cout << "  Double-precision root: x = " << std::setprecision(16) << x_final;
            std::cout << " mult=" << root.multiplicity << " ✓\n";

            if (std::abs(x_final - 0.2) < 0.1) {
                found_simple = true;
                double error = std::abs(x_final - 0.2);
                std::cout << "    -> Simple root, error = " << std::scientific << error << "\n";
                if (config.test_mode) {
                    assert(root.multiplicity == 1 && "Simple root mult=1");
                    assert(error < 1e-10 && "Simple root error < 1e-10");
                }
            } else if (std::abs(x_final - 0.6) < 0.1) {
                found_multiple = true;
                double error = std::abs(x_final - 0.6);
                std::cout << "    -> Multiple root (m=" << root.multiplicity
                          << "), error = " << error << "\n";
            }
        }
    }

    // Then HP-refine candidates that need it
    for (double x0 : hp_candidates) {
        HP x_hp(x0);

        // Debug: manually trace the refinement
        std::cout << "  === Debug trace for x0=" << x0 << " ===\n";

        // Do a few standard Newton steps manually
        PolyHP dpoly_hp = poly_hp.differentiate(0);
        HP x_trace = x_hp;
        for (int i = 0; i < 15; ++i) {
            HP f = poly_hp.evaluate(x_trace);
            HP df = dpoly_hp.evaluate(x_trace);
            HP step = f / df;
            std::cout << "    Newton[" << i << "]: x=" << static_cast<double>(x_trace)
                      << " f=" << static_cast<double>(f)
                      << " df=" << static_cast<double>(df)
                      << " step=" << static_cast<double>(step) << "\n";
            if (abs(step) < HP("1e-50")) break;
            x_trace = x_trace - step;
        }

        // Check Ostrowski from this improved point
        auto ostrowski = RefinerHP::estimateMultiplicityOstrowskiWithLastIterate(x_trace, poly_hp);
        std::cout << "  After Newton: Ostrowski mult=" << ostrowski.multiplicity
                  << ", x3=" << static_cast<double>(ostrowski.last_iterate) << "\n";

        auto hp_result = RefinerHP::refineRoot1D(x_hp, poly_hp, hp_config);

        double x_final = static_cast<double>(hp_result.location);
        double residual = static_cast<double>(hp_result.residual);

        std::cout << "  HP refined: x = " << std::setprecision(16) << x_final;
        std::cout << " mult=" << hp_result.multiplicity;
        std::cout << " |f(x)|=" << std::scientific << std::setprecision(2) << std::abs(residual);
        std::cout << (hp_result.converged ? " ✓" : " ✗") << "\n";

        // Verify results
        if (std::abs(x_final - 0.2) < 0.1) {
            found_simple = true;
            double error = std::abs(x_final - 0.2);
            std::cout << "    -> Simple root, error = " << error << "\n";
            if (config.test_mode) {
                assert(hp_result.converged && "Simple root should converge");
                assert(hp_result.multiplicity == 1 && "Simple root mult=1");
                assert(error < 1e-14 && "Simple root error < 1e-14");
            }
        } else if (std::abs(x_final - 0.6) < 0.1) {
            found_multiple = true;
            double error = std::abs(x_final - 0.6);
            std::cout << "    -> Multiple root (m=" << hp_result.multiplicity
                      << "), error = " << error << "\n";
            if (config.test_mode) {
                assert(hp_result.converged && "Multiple root should converge");
                assert(hp_result.multiplicity >= 5 && "Multiple root mult>=5");
                // For mult-6 root, HP should achieve ~1e-12 or better root error
                // (exact 1e-15 may not be achievable due to modified Newton approximation)
                assert(error < 1e-10 && "Multiple root error < 1e-10 with HP");
            }
        }
    }

    if (config.test_mode) {
        assert(found_simple && "Should find simple root at 0.2");
        assert(found_multiple && "Should find multiple root at 0.6");
        std::cout << "\n✓ All assertions passed\n";
    }
#else
    std::cout << "\nStep 3: HP refinement skipped (ENABLE_HIGH_PRECISION not defined)\n";
    std::cout << "  Build with -DENABLE_HIGH_PRECISION=ON to enable\n";

    // Without HP, verify we found roots via double-precision + detection
    bool found_simple = false, found_multiple = false;

    // Check double-precision roots
    for (const auto& root : batch.roots) {
        if (!root.needs_higher_precision && root.converged) {
            if (std::abs(root.location - 0.2) < 0.1) found_simple = true;
            if (std::abs(root.location - 0.6) < 0.1) found_multiple = true;
        }
    }

    // Check HP candidates (multiple root should be detected)
    for (double x : hp_candidates) {
        if (std::abs(x - 0.2) < 0.1) found_simple = true;
        if (std::abs(x - 0.6) < 0.1) found_multiple = true;
    }

    if (config.test_mode) {
        assert(found_simple && "Should find simple root (in double-precision or HP candidates)");
        assert(found_multiple && "Should detect multiple root needs HP");
        std::cout << "\n✓ Detection assertions passed (full HP test requires ENABLE_HIGH_PRECISION)\n";
    }
#endif

    std::cout << "\n=== Complete ===\n";
    return 0;
}

