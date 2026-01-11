#include "core/polynomial.h"
#include "solver/solver.h"
#include "core/geometry.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

using namespace polynomial_solver;

namespace {

bool approx_equal(double a, double b,
                  double eps = std::numeric_limits<double>::epsilon() * 1000.0)
{
    return std::fabs(a - b) <= eps;
}

int test_subdivision_solver_uniform_grid()
{
    // 1D polynomial p(t) = t on [0,1]. The current subdivision solver does not
    // yet use function values in its decisions, but we attach a simple
    // normalized system so that the workflow is exercised end-to-end.
    std::vector<unsigned int> degrees{1u};
    std::vector<double> power_coeffs{0.0, 1.0}; // p(t) = t
    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    PolynomialSystem system({p});
    if (system.dimension() != 1u || system.equationCount() != 1u) {
        std::cerr << "Subdivision solver: unexpected system shape" << '\n';
        return 1;
    }

    Solver solver;

    SubdivisionConfig config;
    config.tolerance = 0.25;       // Stop once boxes reach width 1/4.
    config.max_depth = 10u;        // Large enough to not be active here.

    const SubdivisionSolverResult result =
        solver.subdivisionSolve(system, config, RootBoundingMethod::None);

    if (result.boxes.size() != 4u) {
        std::cerr << "Subdivision solver: expected 4 final boxes, got "
                  << result.boxes.size() << '\n';
        return 1;
    }

    if (result.num_resolved != 4u) {
        std::cerr << "Subdivision solver: expected 4 resolved boxes, got "
                  << result.num_resolved << '\n';
        return 1;
    }

    const double expected_width = 0.25;
    const double eps = std::numeric_limits<double>::epsilon() * 1000.0;

    for (std::size_t i = 0; i < result.boxes.size(); ++i) {
        const SubdivisionBoxResult& box = result.boxes[i];

        if (box.lower.size() != 1u || box.upper.size() != 1u) {
            std::cerr << "Subdivision solver: box dimension mismatch" << '\n';
            return 1;
        }

        const double width = box.upper[0] - box.lower[0];
        if (!approx_equal(width, expected_width, eps)) {
            std::cerr << "Subdivision solver: unexpected box width " << width
                      << ", expected " << expected_width << '\n';
            return 1;
        }

        if (!box.converged) {
            std::cerr << "Subdivision solver: box not marked as converged" << '\n';
            return 1;
        }

        if (box.depth != 2u) {
            std::cerr << "Subdivision solver: expected depth 2, got "
                      << box.depth << '\n';
            return 1;
        }

        if (box.lower[0] < -eps || box.upper[0] > 1.0 + eps) {
            std::cerr << "Subdivision solver: box outside [0,1]" << '\n';
            return 1;
        }
    }

    return 0;
}

// Removed: test_subdivision_solver_graph_hull_1d() - now covered by test_linear_graph_hull.cpp

// Removed: test_subdivision_solver_graph_hull_1d_debug() - debug CSV dump not needed

// Removed: test_subdivision_solver_graph_hull_2d() - now covered by test_linear_graph_hull.cpp

// Removed: test_subdivision_solver_graph_hull_2d_debug() - debug CSV dump not needed

} // namespace

int test_subdivision_solver_graph_hull_2d_quadratic()
{
    // 2D system with one quadratic equation:
    //   f1(x, y) = x^2 - 0.25  (root at x = 0.5)
    //   f2(x, y) = y - 0.3     (root at y = 0.3)
    // Root at (0.5, 0.3).

    // f1 = x^2 - 0.25 in power basis with degrees (2, 1)
    // Layout (last dimension fastest): [x^0y^0, x^0y^1, x^1y^0, x^1y^1, x^2y^0, x^2y^1]
    // f1(x,y) = -0.25 + 0*y + 0*x + 0*xy + 1*x^2 + 0*x^2y
    std::vector<unsigned int> degrees_f1{2u, 1u};
    std::vector<double> power_coeffs_f1{-0.25, 0.0, 0.0, 0.0, 1.0, 0.0};
    Polynomial p1 = Polynomial::fromPower(degrees_f1, power_coeffs_f1);

    // f2 = y - 0.3 in power basis: [-0.3, 1] with degrees (0, 1)
    std::vector<unsigned int> degrees_f2{0u, 1u};
    std::vector<double> power_coeffs_f2{-0.3, 1.0};
    Polynomial p2 = Polynomial::fromPower(degrees_f2, power_coeffs_f2);

    PolynomialSystem system({p1, p2});

    Solver solver;

    SubdivisionConfig config;
    config.tolerance = 0.05;
    config.max_depth = 10u;

    // First solve with None method to get baseline.
    const SubdivisionSolverResult result_none =
        solver.subdivisionSolve(system, config, RootBoundingMethod::None);

    // Then solve with GraphHull method.
    const SubdivisionSolverResult result_hull =
        solver.subdivisionSolve(system, config, RootBoundingMethod::GraphHull);

    if (result_hull.boxes.empty()) {
        std::cerr << "GraphHull 2D quadratic: expected at least one box" << '\n';
        return 1;
    }

    // With GraphHull, we should get fewer or equal number of boxes.
    if (result_hull.boxes.size() > result_none.boxes.size()) {
        std::cerr << "GraphHull 2D quadratic: expected fewer boxes than None method, got "
                  << result_hull.boxes.size() << " vs " << result_none.boxes.size() << '\n';
        return 1;
    }

    // Check that at least one box contains the root (0.5, 0.3).
    bool found_root = false;
    for (const SubdivisionBoxResult& box : result_hull.boxes) {
        if (box.lower[0] <= 0.5 && box.upper[0] >= 0.5 &&
            box.lower[1] <= 0.3 && box.upper[1] >= 0.3) {
            found_root = true;
            break;
        }
    }

    if (!found_root) {
        std::cerr << "GraphHull 2D quadratic: root (0.5, 0.3) not found in any box" << '\n';
        return 1;
    }

    std::cout << "GraphHull 2D quadratic: None method produced " << result_none.boxes.size()
              << " boxes, GraphHull produced " << result_hull.boxes.size() << " boxes" << std::endl;

    return 0;
}

// Removed: test_subdivision_solver_graph_hull_2d_quadratic_debug() - debug CSV dump not needed

int main(int argc, char* argv[])
{
    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--debug") == 0) {
            getGeometryConfig().debug = true;
            std::cout << "Debug mode enabled" << std::endl;
        }
    }

    int status = 0;
    status |= test_subdivision_solver_uniform_grid();
    status |= test_subdivision_solver_graph_hull_2d_quadratic();

    if (status == 0) {
        std::cout << "Subdivision solver tests passed." << std::endl;
    }

    return status;
}

