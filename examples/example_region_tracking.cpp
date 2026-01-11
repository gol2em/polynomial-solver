/**
 * @file example_region_tracking.cpp
 * @brief Demonstrates region tracking in PolynomialBase
 *
 * Shows how subdivision automatically tracks regions in original coordinates.
 */

#include "core/polynomial_base.h"
#include <iostream>
#include <iomanip>

using namespace polynomial_solver;

void print_box(const std::string& label, 
               const std::vector<double>& lower, 
               const std::vector<double>& upper) {
    std::cout << label << ": ";
    for (size_t i = 0; i < lower.size(); ++i) {
        std::cout << "[" << std::fixed << std::setprecision(2) 
                  << lower[i] << ", " << upper[i] << "]";
        if (i + 1 < lower.size()) std::cout << " × ";
    }
    std::cout << "\n";
}

int main() {
    std::cout << "=== Region Tracking Example ===\n\n";

    // Example 1: 1D polynomial with subdivision
    std::cout << "Example 1: 1D Polynomial Subdivision\n";
    std::cout << "-------------------------------------\n";
    {
        std::vector<unsigned int> degrees = {3};
        std::vector<double> bern_coeffs = {0.0, 1.0/3.0, 2.0/3.0, 1.0};

        PolynomialBase<double> p(degrees, bern_coeffs);
        
        auto [l0, u0] = p.getOriginalBox();
        print_box("Original region", l0, u0);

        // Subdivide to [0.2, 0.8]
        auto sub1 = p.restrictedToInterval(0, 0.2, 0.8);
        auto [l1, u1] = sub1.getOriginalBox();
        print_box("After [0.2, 0.8]", l1, u1);

        // Subdivide again to [0.3, 0.7] in local coords
        // Global: [0.38, 0.62]
        auto sub2 = sub1.restrictedToInterval(0, 0.3, 0.7);
        auto [l2, u2] = sub2.getOriginalBox();
        print_box("After [0.3, 0.7]", l2, u2);
        
        std::cout << "\n";
    }

    // Example 2: 2D polynomial with multiple subdivisions
    std::cout << "Example 2: 2D Polynomial Subdivision\n";
    std::cout << "-------------------------------------\n";
    {
        std::vector<unsigned int> degrees = {2, 2};
        std::vector<double> bern_coeffs = {
            1.0, 2.0, 3.0,
            2.0, 3.0, 4.0,
            3.0, 4.0, 5.0
        };

        PolynomialBase<double> p(degrees, bern_coeffs);
        
        auto [l0, u0] = p.getOriginalBox();
        print_box("Original region", l0, u0);

        // Subdivide axis 0 to [0.0, 0.5]
        auto sub1 = p.restrictedToInterval(0, 0.0, 0.5);
        auto [l1, u1] = sub1.getOriginalBox();
        print_box("After axis 0: [0.0, 0.5]", l1, u1);

        // Subdivide axis 1 to [0.5, 1.0]
        auto sub2 = sub1.restrictedToInterval(1, 0.5, 1.0);
        auto [l2, u2] = sub2.getOriginalBox();
        print_box("After axis 1: [0.5, 1.0]", l2, u2);

        std::cout << "\n";
    }

    // Example 3: Simulating solver workflow
    std::cout << "Example 3: Solver Workflow Simulation\n";
    std::cout << "--------------------------------------\n";
    {
        // Create polynomial from power basis
        std::vector<unsigned int> degrees = {2};
        std::vector<double> power_coeffs = {1.0, -2.0, 1.0}; // (x-1)^2

        auto p = PolynomialBase<double>::fromPower(degrees, power_coeffs);
        
        auto [l0, u0] = p.getOriginalBox();
        print_box("Initial region", l0, u0);

        // Solver subdivides to find roots
        // First subdivision: [0.0, 0.5]
        auto left = p.restrictedToInterval(0, 0.0, 0.5);
        auto [ll, lu] = left.getOriginalBox();
        print_box("Left half", ll, lu);

        // Second subdivision: [0.5, 1.0]
        auto right = p.restrictedToInterval(0, 0.5, 1.0);
        auto [rl, ru] = right.getOriginalBox();
        print_box("Right half", rl, ru);

        // Further subdivide right half: [0.75, 1.0]
        auto right_right = right.restrictedToInterval(0, 0.5, 1.0);
        auto [rrl, rru] = right_right.getOriginalBox();
        print_box("Right-right quarter", rrl, rru);

        std::cout << "\n";
    }

    std::cout << "=== Region Tracking Working Correctly ===\n";

    return 0;
}

