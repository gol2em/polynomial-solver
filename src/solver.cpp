#include "solver.h"

/**
 * @file solver.cpp
 * @brief Implementation of the main solver interface
 */

namespace polynomial_solver {

PolynomialSystem::PolynomialSystem()
    : dimension_(0u)
{
}

PolynomialSystem::PolynomialSystem(const std::vector<Polynomial>& equations)
    : dimension_(0u),
      equations_(equations)
{
    if (!equations_.empty()) {
        dimension_ = equations_[0].dimension();
        // TODO: Optionally verify that all equations share this dimension.
    }
}

std::size_t PolynomialSystem::dimension() const {
    return dimension_;
}

std::size_t PolynomialSystem::equationCount() const {
    return equations_.size();
}

const Polynomial& PolynomialSystem::equation(std::size_t i) const {
    return equations_[i];
}

const std::vector<Polynomial>& PolynomialSystem::equations() const {
    return equations_;
}
std::vector<GraphControlNet> PolynomialSystem::graphControlNets() const {
    std::vector<GraphControlNet> nets;
    nets.reserve(equations_.size());

    for (const Polynomial& poly : equations_) {
        GraphControlNet net;
        net.degrees = poly.degrees();
        poly.graphControlPoints(net.control_points);
        nets.push_back(std::move(net));
    }

    return nets;
}




Solver::Solver() {
    // TODO: Implement constructor
}

Solver::~Solver() {
    // TODO: Implement destructor
}

// TODO: Implement solver operations

} // namespace polynomial_solver

