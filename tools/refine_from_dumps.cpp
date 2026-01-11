/**
 * @file refine_from_dumps.cpp
 * @brief Refine roots from example result dumps
 * 
 * This tool loads solver results from dumps/ directory and performs
 * systematic 1D root refinement:
 * 1. Refine resolved boxes to 1e-15 precision
 * 2. Verify roots and determine multiplicity
 * 3. Compute exclusion radius based on derivative
 * 4. Merge unresolved boxes into problematic regions
 */

#include "refinement/result_refiner.h"
#include "core/polynomial.h"
#include "solver/solver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace polynomial_solver;

/**
 * @brief Parse a result dump file
 */
SubdivisionSolverResult parseDumpFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    SubdivisionSolverResult result;
    std::string line;
    int total_boxes = 0;
    int resolved_count = 0;
    int current_box = -1;
    bool in_resolved = false;
    
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            if (line.find("# Total boxes:") != std::string::npos) {
                sscanf(line.c_str(), "# Total boxes: %d", &total_boxes);
            } else if (line.find("# Resolved:") != std::string::npos) {
                sscanf(line.c_str(), "# Resolved: %d", &resolved_count);
            }
            continue;
        }
        
        // Section headers
        if (line.find("## Resolved Boxes") != std::string::npos) {
            in_resolved = true;
            continue;
        } else if (line.find("## Unresolved Boxes") != std::string::npos) {
            in_resolved = false;
            continue;
        }
        
        // Box header
        if (line.find("Box ") == 0) {
            result.boxes.push_back(SubdivisionBoxResult());
            current_box = result.boxes.size() - 1;
            continue;
        }

        // Parse box properties
        if (current_box >= 0 && line.find("  ") == 0) {
            SubdivisionBoxResult& box = result.boxes[current_box];
            
            if (line.find("Lower:") != std::string::npos) {
                double val;
                sscanf(line.c_str(), "  Lower: %lf", &val);
                box.lower.push_back(val);
            } else if (line.find("Upper:") != std::string::npos) {
                double val;
                sscanf(line.c_str(), "  Upper: %lf", &val);
                box.upper.push_back(val);
            } else if (line.find("Center:") != std::string::npos) {
                double val;
                sscanf(line.c_str(), "  Center: %lf", &val);
                box.center.push_back(val);
            } else if (line.find("Depth:") != std::string::npos) {
                sscanf(line.c_str(), "  Depth: %u", &box.depth);
            }
        }
    }
    
    // Compute max_error for each box
    for (auto& box : result.boxes) {
        box.max_error.resize(box.lower.size());
        for (std::size_t i = 0; i < box.lower.size(); ++i) {
            box.max_error[i] = (box.upper[i] - box.lower[i]) / 2.0;
        }
    }
    
    result.num_resolved = resolved_count;
    return result;
}

/**
 * @brief Merge nearby unresolved boxes into problematic regions
 *
 * Note: This function is kept for backward compatibility with the tool.
 * The main ResultRefiner class now has mergeUnverifiedBoxes1D() which
 * provides similar functionality.
 */
std::vector<polynomial_solver::ProblematicRegion> mergeUnresolvedBoxes(
    const SubdivisionSolverResult& result,
    double merge_threshold = 1e-6) {

    std::vector<polynomial_solver::ProblematicRegion> regions;
    
    // Get unresolved boxes
    std::vector<std::size_t> unresolved_indices;
    for (std::size_t i = result.num_resolved; i < result.boxes.size(); ++i) {
        unresolved_indices.push_back(i);
    }
    
    if (unresolved_indices.empty()) {
        return regions;
    }
    
    // Sort by center position
    std::sort(unresolved_indices.begin(), unresolved_indices.end(),
        [&result](std::size_t a, std::size_t b) {
            return result.boxes[a].center[0] < result.boxes[b].center[0];
        });
    
    // Merge nearby boxes
    polynomial_solver::ProblematicRegion current;
    current.lower = result.boxes[unresolved_indices[0]].lower[0];
    current.upper = result.boxes[unresolved_indices[0]].upper[0];
    current.box_indices.push_back(unresolved_indices[0]);
    current.multiplicity = 0;
    current.refined_root = 0.0;
    current.residual = 0.0;
    current.refinement_succeeded = false;
    
    for (std::size_t i = 1; i < unresolved_indices.size(); ++i) {
        std::size_t idx = unresolved_indices[i];
        const auto& box = result.boxes[idx];
        
        // Check if this box is close to current region
        if (box.lower[0] - current.upper <= merge_threshold) {
            // Merge into current region
            current.upper = std::max(current.upper, box.upper[0]);
            current.box_indices.push_back(idx);
        } else {
            // Start new region
            regions.push_back(current);
            current = polynomial_solver::ProblematicRegion();
            current.lower = box.lower[0];
            current.upper = box.upper[0];
            current.box_indices.push_back(idx);
        }
    }
    regions.push_back(current);
    
    return regions;
}

/**
 * @brief Estimate multiplicity for a problematic region and refine using modified Newton
 *
 * Note: This function is kept for backward compatibility. The main ResultRefiner
 * class now has refineProblematicRegion1D() which provides similar functionality.
 */
void refineProblematicRegion(
    polynomial_solver::ProblematicRegion& region,
    const Polynomial& poly,
    const ResultRefiner& refiner,
    const RefinementConfig& config)
{
    double center = (region.lower + region.upper) / 2.0;

    // Estimate multiplicity by checking derivatives at center
    PolynomialSystem system({poly});
    double first_nonzero_deriv = 0.0;
    unsigned int mult = refiner.estimateMultiplicity(
        std::vector<double>{center}, system, config.max_multiplicity, 1e-10, first_nonzero_deriv);

    region.multiplicity = mult;

    // Try to refine using modified Newton method with estimated multiplicity
    double refined_location;
    double residual;

    bool success = refiner.refineRoot1DWithMultiplicity(
        center, region.lower, region.upper, poly, mult, config, refined_location, residual);

    if (success) {
        region.refined_root = refined_location;
        region.residual = residual;
        region.refinement_succeeded = true;
    } else {
        // If failed with estimated multiplicity, try different multiplicities
        for (unsigned int m = 1; m <= config.max_multiplicity; ++m) {
            if (m == mult) continue;  // Already tried this

            success = refiner.refineRoot1DWithMultiplicity(
                center, region.lower, region.upper, poly, m, config, refined_location, residual);

            if (success) {
                region.refined_root = refined_location;
                region.residual = residual;
                region.refinement_succeeded = true;
                region.multiplicity = m;
                break;
            }
        }
    }
}

/**
 * @brief Test refinement on a specific example
 */
void testExample(const std::string& name,
                 const std::string& dump_file,
                 const Polynomial& poly) {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "Example: " << name << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    // Load result from dump
    SubdivisionSolverResult result;
    try {
        result = parseDumpFile(dump_file);
    } catch (const std::exception& e) {
        std::cerr << "Error loading dump: " << e.what() << std::endl;
        return;
    }

    std::cout << "Loaded from dump: " << dump_file << std::endl;
    std::cout << "  Total boxes: " << result.boxes.size() << std::endl;
    std::cout << "  Resolved: " << result.num_resolved << std::endl;
    std::cout << "  Unresolved: " << (result.boxes.size() - result.num_resolved) << std::endl;

    // Create polynomial system
    PolynomialSystem system({poly});

    // Refine resolved boxes
    ResultRefiner refiner;
    RefinementConfig config;
    config.target_tolerance = 1e-15;
    config.residual_tolerance = 1e-12;
    config.max_multiplicity = 10;
    config.exclusion_multiplier = 3.0;

    RefinementResult refined = refiner.refine(result, system, config);

    std::cout << "\n--- Refinement Results ---" << std::endl;
    std::cout << "Verified roots: " << refined.roots.size() << std::endl;
    std::cout << "Cancelled boxes: " << refined.cancelled_boxes.size() << std::endl;
    std::cout << "Unverified boxes: " << refined.unverified_boxes.size() << std::endl;

    // Print unverified boxes details
    if (!refined.unverified_boxes.empty()) {
        std::cout << "\n⚠️  Unverified boxes (failed refinement):" << std::endl;
        for (std::size_t idx : refined.unverified_boxes) {
            const auto& box = result.boxes[idx];
            std::cout << "  Box " << idx << ": center = " << std::setprecision(16)
                      << box.center[0] << ", max_error = " << std::scientific
                      << box.max_error[0] << std::endl;

            // Try to evaluate at center
            double val = poly.evaluate(box.center);
            std::cout << "    f(center) = " << std::scientific << std::setprecision(4)
                      << val << std::endl;
        }
    }

    // Print verified roots
    for (std::size_t i = 0; i < refined.roots.size(); ++i) {
        const auto& root = refined.roots[i];
        std::cout << "\nRoot " << i << ":" << std::endl;
        std::cout << "  Location: x = " << std::setprecision(16) << root.location[0] << std::endl;
        std::cout << "  Residual: |f(x)| = " << std::scientific << std::setprecision(4)
                  << std::abs(root.residual[0]) << std::endl;
        std::cout << "  Multiplicity: " << root.multiplicity << std::endl;
        std::cout << "  First non-zero derivative: " << std::scientific << std::setprecision(6)
                  << root.first_nonzero_derivative << std::endl;
        std::cout << "  Exclusion radius: " << std::scientific << std::setprecision(6)
                  << root.exclusion_radius << std::endl;
        std::cout << "  Source boxes: " << root.source_boxes.size() << std::endl;
        std::cout << "  Max error: " << std::scientific << root.max_error[0] << std::endl;
        std::cout << "  Condition estimate: " << std::scientific << std::setprecision(4)
                  << root.condition_estimate << std::endl;

        // Estimate actual error from condition number
        double estimated_error = root.condition_estimate * std::abs(root.residual[0]) /
                                 std::max(std::abs(root.first_nonzero_derivative), 1e-14);
        std::cout << "  Estimated error: " << std::scientific << std::setprecision(4)
                  << estimated_error << std::endl;

        // Check if refinement succeeded
        if (std::abs(root.residual[0]) > 1e-12) {
            std::cout << "  ⚠️  WARNING: Residual exceeds tolerance!" << std::endl;
        }
        if (root.max_error[0] > 1e-14) {
            std::cout << "  ⚠️  WARNING: Error exceeds target precision!" << std::endl;
        }
        if (root.needs_higher_precision) {
            std::cout << "  ⚠️  WARNING: Condition number suggests higher precision needed!" << std::endl;
            std::cout << "      Small residual does NOT guarantee accurate root for this problem." << std::endl;
            std::cout << "      Consider using quad precision or arbitrary precision arithmetic." << std::endl;
        }
    }

    // Handle unresolved boxes
    if (result.boxes.size() > result.num_resolved) {
        std::cout << "\n--- Unresolved Boxes Analysis ---" << std::endl;
        auto regions = mergeUnresolvedBoxes(result, 1e-6);
        std::cout << "Merged into " << regions.size() << " problematic region(s):" << std::endl;

        // Refine each problematic region
        ResultRefiner refiner;
        RefinementConfig region_config;
        region_config.target_tolerance = 1e-15;
        region_config.residual_tolerance = 1e-12;
        region_config.max_newton_iters = 100;  // More iterations for difficult cases
        region_config.max_multiplicity = 10;

        std::size_t successful_refinements = 0;

        for (std::size_t i = 0; i < regions.size(); ++i) {
            auto& region = regions[i];

            // Try to refine this region
            refineProblematicRegion(region, poly, refiner, region_config);

            if (region.refinement_succeeded) {
                successful_refinements++;
            }
        }

        std::cout << "\nRefined " << successful_refinements << " out of "
                  << regions.size() << " problematic regions." << std::endl;

        // Print detailed results
        for (std::size_t i = 0; i < regions.size(); ++i) {
            const auto& region = regions[i];
            std::cout << "\nProblematic Region " << i << ":" << std::endl;
            std::cout << "  Interval: [" << std::setprecision(16) << region.lower
                      << ", " << region.upper << "]" << std::endl;
            std::cout << "  Width: " << std::scientific << std::setprecision(4)
                      << (region.upper - region.lower) << std::endl;
            std::cout << "  Number of boxes: " << region.box_indices.size() << std::endl;

            if (region.refinement_succeeded) {
                std::cout << "  ✅ Refinement: SUCCESS" << std::endl;
                std::cout << "  Root location: x = " << std::fixed << std::setprecision(16)
                          << region.refined_root << std::endl;
                std::cout << "  Residual: |f(x)| = " << std::scientific << std::setprecision(4)
                          << std::abs(region.residual) << std::endl;
                std::cout << "  Estimated multiplicity: " << region.multiplicity << std::endl;
            } else {
                std::cout << "  ❌ Refinement: FAILED" << std::endl;
                std::cout << "  Center: " << std::fixed << std::setprecision(16)
                          << (region.lower + region.upper) / 2.0 << std::endl;
                double center = (region.lower + region.upper) / 2.0;
                double val = poly.evaluate({center});
                std::cout << "  f(center): " << std::scientific << std::setprecision(4) << val << std::endl;
                std::cout << "  Estimated multiplicity: " << region.multiplicity << std::endl;
            }
        }
    }
}

int main() {
    std::cout << "1D Root Refinement from Example Dumps" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    // Test 1: Cubic polynomial
    {
        std::vector<unsigned int> degrees{3u};
        std::vector<double> power_coeffs{-0.08, 0.66, -1.5, 1.0};
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        testExample("Cubic: (x-0.2)(x-0.5)(x-0.8)",
                    "dumps/cubic_1d_result.txt", p);
    }

    // Test 2: Multiplicity polynomial
    {
        std::vector<unsigned int> degrees{7u};
        std::vector<double> power_coeffs{
            -0.0093312, 0.139968, -0.85536, 2.808, -5.4, 6.12, -3.8, 1.0
        };
        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        testExample("Multiplicity: (x-0.2)(x-0.6)^6",
                    "dumps/multiplicity_1d_result.txt", p);
    }

    // Test 3: Wilkinson-19 polynomial
    {
        // Build (x-1/20)(x-2/20)...(x-19/20) = (x-0.05)(x-0.10)...(x-0.95)
        // Start with (x - 1/20)
        double scale = 1.0 / 20.0;
        std::vector<double> coeffs = {-scale, 1.0};  // -1/20 + x

        // Multiply by (x - k/20) for k = 2, 3, ..., 19
        for (int k = 2; k <= 19; ++k) {
            std::vector<double> new_coeffs(coeffs.size() + 1, 0.0);
            double root = k * scale;

            // Multiply: coeffs * (x - k/20) = coeffs * x - (k/20) * coeffs
            for (std::size_t i = 0; i < coeffs.size(); ++i) {
                new_coeffs[i] -= root * coeffs[i];      // -(k/20) * coeffs
                new_coeffs[i + 1] += coeffs[i];         // coeffs * x
            }

            coeffs = new_coeffs;
        }

        std::vector<unsigned int> degrees{19u};
        Polynomial p = Polynomial::fromPower(degrees, coeffs);
        testExample("Wilkinson-19: roots at 0.05, 0.10, ..., 0.95",
                    "dumps/wilkinson_1d_result.txt", p);
    }

    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "All examples completed!" << std::endl;

    return 0;
}

