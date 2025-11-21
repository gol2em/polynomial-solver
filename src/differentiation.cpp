#include "differentiation.h"
#include <algorithm>
#include <functional>
#include <stdexcept>

/**
 * @file differentiation.cpp
 * @brief Implementation of differentiation utilities for Bernstein polynomials
 */

namespace polynomial_solver {

namespace {

// Helper: Compute strides for tensor-product layout (last dimension fastest)
void compute_strides(const std::vector<unsigned int>& degrees,
                     std::vector<std::size_t>& strides)
{
    const std::size_t dim = degrees.size();
    strides.assign(dim, 0u);

    std::size_t s = 1u;
    for (std::size_t k = dim; k-- > 0;) {
        strides[k] = s;
        s *= static_cast<std::size_t>(degrees[k] + 1u);
    }
}

// Helper: Flatten multi-index to linear index
std::size_t flatten_index(const std::vector<unsigned int>& multi_index,
                          const std::vector<std::size_t>& strides)
{
    std::size_t idx = 0u;
    for (std::size_t d = 0; d < multi_index.size(); ++d) {
        idx += static_cast<std::size_t>(multi_index[d]) * strides[d];
    }
    return idx;
}

// Helper: Increment multi-index except along specified axis
bool increment_multi_except_axis(std::vector<unsigned int>& multi_index,
                                 const std::vector<unsigned int>& degrees,
                                 std::size_t skip_axis)
{
    const std::size_t dim = degrees.size();
    for (std::size_t d = dim; d-- > 0;) {
        if (d == skip_axis) {
            continue;
        }
        if (multi_index[d] < degrees[d]) {
            ++multi_index[d];
            for (std::size_t e = d + 1; e < dim; ++e) {
                if (e != skip_axis) {
                    multi_index[e] = 0u;
                }
            }
            return true;
        }
    }
    return false;
}

} // anonymous namespace

// ============================================================================
// Differentiation class implementation
// ============================================================================

Polynomial Differentiation::differentiateAxis(const Polynomial& p, std::size_t axis)
{
    const std::size_t dim = p.dimension();
    if (axis >= dim) {
        // Invalid axis: return zero polynomial
        return Polynomial(p.degrees(), std::vector<double>(p.coefficientCount(), 0.0));
    }

    const std::vector<unsigned int>& degrees = p.degrees();
    const std::vector<double>& coeffs = p.bernsteinCoefficients();

    const unsigned int deg_axis = degrees[axis];
    
    // If degree is 0 along this axis, derivative is zero
    if (deg_axis == 0u) {
        return Polynomial(degrees, std::vector<double>(p.coefficientCount(), 0.0));
    }

    // New degrees: reduce by 1 along the differentiation axis
    std::vector<unsigned int> new_degrees = degrees;
    new_degrees[axis] = deg_axis - 1u;

    const std::size_t len_axis = static_cast<std::size_t>(deg_axis + 1u);
    const std::size_t new_len_axis = static_cast<std::size_t>(deg_axis);

    // Compute strides for both old and new layouts
    std::vector<std::size_t> old_strides, new_strides;
    compute_strides(degrees, old_strides);
    compute_strides(new_degrees, new_strides);

    // Compute new coefficient count
    std::size_t new_count = 1u;
    for (std::size_t i = 0; i < dim; ++i) {
        new_count *= static_cast<std::size_t>(new_degrees[i] + 1u);
    }

    std::vector<double> new_coeffs(new_count, 0.0);

    // Iterate over all 1D slices along the given axis
    std::vector<unsigned int> multi_index(dim, 0u);
    std::vector<double> line(len_axis);
    std::vector<double> deriv_line(new_len_axis);

    std::fill(multi_index.begin(), multi_index.end(), 0u);
    bool first = true;
    while (first || increment_multi_except_axis(multi_index, degrees, axis)) {
        first = false;

        // Gather coefficients along the current axis-line
        for (std::size_t k = 0; k < len_axis; ++k) {
            multi_index[axis] = static_cast<unsigned int>(k);
            const std::size_t idx = flatten_index(multi_index, old_strides);
            line[k] = coeffs[idx];
        }

        // Apply Bernstein derivative formula: d_i = n * (b_{i+1} - b_i)
        const double n = static_cast<double>(deg_axis);
        for (std::size_t i = 0; i < new_len_axis; ++i) {
            deriv_line[i] = n * (line[i + 1] - line[i]);
        }

        // Scatter the derivative coefficients back into the tensor
        for (std::size_t k = 0; k < new_len_axis; ++k) {
            multi_index[axis] = static_cast<unsigned int>(k);
            const std::size_t idx = flatten_index(multi_index, new_strides);
            new_coeffs[idx] = deriv_line[k];
        }
    }

    return Polynomial(new_degrees, new_coeffs);
}

Polynomial Differentiation::derivative(const Polynomial& p, std::size_t axis, unsigned int order)
{
    if (order == 0u) {
        return p;
    }

    // Iteratively apply first-order differentiation
    Polynomial result = p;
    for (unsigned int k = 0; k < order; ++k) {
        result = differentiateAxis(result, axis);
    }

    return result;
}

std::vector<Polynomial> Differentiation::gradient(const Polynomial& p)
{
    const std::size_t dim = p.dimension();
    std::vector<Polynomial> grad;
    grad.reserve(dim);

    for (std::size_t i = 0; i < dim; ++i) {
        grad.push_back(differentiateAxis(p, i));
    }

    return grad;
}

std::vector<std::vector<Polynomial>> Differentiation::hessian(const Polynomial& p)
{
    const std::size_t dim = p.dimension();
    std::vector<std::vector<Polynomial>> hess(dim);

    for (std::size_t i = 0; i < dim; ++i) {
        hess[i].resize(dim);
        Polynomial df_dxi = differentiateAxis(p, i);
        
        for (std::size_t j = 0; j < dim; ++j) {
            hess[i][j] = differentiateAxis(df_dxi, j);
        }
    }

    return hess;
}

// ============================================================================
// DerivativeCache class implementation
// ============================================================================

DerivativeCache::DerivativeCache(const Polynomial& p)
    : dimension_(p.dimension())
{
    // Store the original polynomial as derivative order {0, 0, ..., 0}
    std::vector<unsigned int> zero_orders(dimension_, 0u);
    cache_[zero_orders] = p;
}

const Polynomial& DerivativeCache::get(const std::vector<unsigned int>& orders)
{
    if (orders.size() != dimension_) {
        throw std::invalid_argument("DerivativeCache::get: orders size must match dimension");
    }

    return computeIfNeeded(orders);
}

const Polynomial& DerivativeCache::getPartial(std::size_t axis, unsigned int order)
{
    if (axis >= dimension_) {
        throw std::invalid_argument("DerivativeCache::getPartial: axis out of range");
    }

    std::vector<unsigned int> orders(dimension_, 0u);
    orders[axis] = order;
    return get(orders);
}

void DerivativeCache::precomputeUpToOrder(unsigned int maxOrder)
{
    // Generate all multi-indices with sum(orders) <= maxOrder
    // Use recursive enumeration
    std::vector<unsigned int> orders(dimension_, 0u);

    // Helper lambda for recursive enumeration
    std::function<void(std::size_t, unsigned int)> enumerate;
    enumerate = [&](std::size_t axis, unsigned int remaining) {
        if (axis == dimension_) {
            // Base case: compute this derivative
            computeIfNeeded(orders);
            return;
        }

        // Try all possible orders for this axis
        for (unsigned int k = 0; k <= remaining; ++k) {
            orders[axis] = k;
            enumerate(axis + 1, remaining - k);
        }
    };

    enumerate(0, maxOrder);
}

std::size_t DerivativeCache::dimension() const
{
    return dimension_;
}

const Polynomial& DerivativeCache::computeIfNeeded(const std::vector<unsigned int>& orders)
{
    // Check if already cached
    auto it = cache_.find(orders);
    if (it != cache_.end()) {
        return it->second;
    }

    // Find which axis to differentiate
    // Strategy: find the first non-zero order, reduce it by 1, and differentiate
    std::vector<unsigned int> prev_orders = orders;
    std::size_t diff_axis = dimension_;

    for (std::size_t i = 0; i < dimension_; ++i) {
        if (orders[i] > 0u) {
            prev_orders[i] = orders[i] - 1u;
            diff_axis = i;
            break;
        }
    }

    if (diff_axis >= dimension_) {
        // All orders are zero - should be the original polynomial
        // This should have been cached in the constructor
        throw std::logic_error("DerivativeCache: zero-order derivative not found");
    }

    // Recursively get the previous derivative
    const Polynomial& prev = computeIfNeeded(prev_orders);

    // Differentiate once more
    Polynomial result = Differentiation::differentiateAxis(prev, diff_axis);

    // Cache and return
    cache_[orders] = result;
    return cache_[orders];
}

} // namespace polynomial_solver

