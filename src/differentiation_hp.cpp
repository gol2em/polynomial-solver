#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "hp/differentiation_hp.h"
#include "hp/high_precision_types.h"
#include <algorithm>
#include <stdexcept>
#include <iostream>

// Temporary debug flag (must match result_refiner_hp.cpp)
// #define DEBUG_HP_REFINER 1

/**
 * @file differentiation_hp.cpp
 * @brief Implementation of high-precision differentiation utilities
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
// DifferentiationHP class implementation
// ============================================================================

PolynomialHP DifferentiationHP::differentiateAxis(const PolynomialHP& p, std::size_t axis)
{
    const std::size_t dim = p.dimension();
    if (axis >= dim) {
        // Invalid axis: return zero polynomial
        std::vector<mpreal> zero_coeffs(p.coefficientCount(), mpreal(0));
        return PolynomialHP(p.degrees(), zero_coeffs);
    }

    const std::vector<unsigned int>& degrees = p.degrees();
    const std::vector<mpreal>& coeffs = p.bernsteinCoefficients();

    const unsigned int deg_axis = degrees[axis];
    
    // If degree is 0 along this axis, derivative is zero
    if (deg_axis == 0u) {
        std::vector<mpreal> zero_coeffs(p.coefficientCount(), mpreal(0));
        return PolynomialHP(degrees, zero_coeffs);
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

    std::vector<mpreal> new_coeffs(new_count, mpreal(0));

    // Iterate over all 1D slices along the given axis
    std::vector<unsigned int> multi_index(dim, 0u);
    std::vector<mpreal> line(len_axis);
    std::vector<mpreal> deriv_line(new_len_axis);

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
        // Use high-precision arithmetic
        const mpreal n = mpreal(deg_axis);
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

    return PolynomialHP(new_degrees, new_coeffs);
}

PolynomialHP DifferentiationHP::differentiateAxisPower(const PolynomialHP& p, std::size_t axis)
{
    const std::size_t dim = p.dimension();
    if (axis >= dim) {
        // Invalid axis: return zero polynomial
        std::vector<mpreal> zero_coeffs(p.coefficientCount(), mpreal(0));
        return fromPowerHP(p.degrees(), zero_coeffs);
    }

    const std::vector<unsigned int>& degrees = p.degrees();
    const std::vector<mpreal>& coeffs = p.powerCoefficients();  // This will error if power not valid

    const unsigned int deg_axis = degrees[axis];

    // If degree is 0 along this axis, derivative is zero
    if (deg_axis == 0u) {
        std::vector<mpreal> zero_coeffs(p.coefficientCount(), mpreal(0));
        return fromPowerHP(degrees, zero_coeffs);
    }

    // New degrees: reduce degree along differentiation axis
    std::vector<unsigned int> new_degrees = degrees;
    new_degrees[axis] = deg_axis - 1u;

    // Compute strides for tensor indexing
    std::vector<std::size_t> old_strides, new_strides;
    old_strides.resize(dim);
    new_strides.resize(dim);

    old_strides[dim - 1] = 1;
    new_strides[dim - 1] = 1;
    for (std::size_t i = dim - 1; i > 0; --i) {
        old_strides[i - 1] = old_strides[i] * static_cast<std::size_t>(degrees[i] + 1u);
        new_strides[i - 1] = new_strides[i] * static_cast<std::size_t>(new_degrees[i] + 1u);
    }

    // Allocate output coefficients
    std::size_t new_count = 1;
    for (unsigned int d : new_degrees) {
        new_count *= static_cast<std::size_t>(d + 1u);
    }
    std::vector<mpreal> new_coeffs(new_count, mpreal(0));

    // Helper: flatten multi-index to linear index
    auto flatten_index = [](const std::vector<unsigned int>& idx, const std::vector<std::size_t>& strides) -> std::size_t {
        std::size_t linear = 0;
        for (std::size_t i = 0; i < idx.size(); ++i) {
            linear += static_cast<std::size_t>(idx[i]) * strides[i];
        }
        return linear;
    };

    // Helper: increment multi-index except along given axis
    auto increment_multi_except_axis = [](std::vector<unsigned int>& idx, const std::vector<unsigned int>& degs, std::size_t skip_axis) -> bool {
        for (std::size_t i = idx.size(); i > 0; --i) {
            std::size_t d = i - 1;
            if (d == skip_axis) continue;
            if (idx[d] < degs[d]) {
                ++idx[d];
                return true;
            }
            idx[d] = 0;
        }
        return false;
    };

    // Temporary storage for 1D slices
    const std::size_t len_axis = static_cast<std::size_t>(deg_axis + 1u);
    const std::size_t new_len_axis = static_cast<std::size_t>(deg_axis);
    std::vector<mpreal> line(len_axis);
    std::vector<mpreal> deriv_line(new_len_axis);

    // Multi-index for iteration
    std::vector<unsigned int> multi_index(dim, 0u);

    // Iterate over all 1D slices along the given axis
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

        // Apply power basis derivative formula: d/dx(a_i * x^i) = i * a_i * x^(i-1)
        // So coefficient of x^i in derivative is (i+1) * a_{i+1} from original
        for (std::size_t i = 0; i < new_len_axis; ++i) {
            deriv_line[i] = mpreal(i + 1) * line[i + 1];
        }

        // Scatter the derivative coefficients back into the tensor
        for (std::size_t k = 0; k < new_len_axis; ++k) {
            multi_index[axis] = static_cast<unsigned int>(k);
            const std::size_t idx = flatten_index(multi_index, new_strides);
            new_coeffs[idx] = deriv_line[k];
        }
    }

    return fromPowerHP(new_degrees, new_coeffs);
}

PolynomialHP DifferentiationHP::derivative(const PolynomialHP& p, std::size_t axis, unsigned int order)
{
    if (order == 0u) {
        return p;
    }

    // Iteratively apply first-order differentiation
    // Choose method based on current polynomial's representation at each step
    PolynomialHP result = p;
    for (unsigned int k = 0; k < order; ++k) {
        // Check representation at each iteration
        bool use_power = (result.primaryRepresentation() == PolynomialRepresentation::POWER &&
                         result.hasPowerCoefficients());

        #ifdef DEBUG_HP_REFINER
        std::cout << "  [HP DEBUG] DifferentiationHP::derivative order " << k+1 << "/" << order
                  << ", using " << (use_power ? "POWER" : "BERNSTEIN") << " basis\n";
        #endif

        if (use_power) {
            result = differentiateAxisPower(result, axis);
        } else {
            result = differentiateAxis(result, axis);
        }
    }

    return result;
}

std::vector<PolynomialHP> DifferentiationHP::gradient(const PolynomialHP& p)
{
    const std::size_t dim = p.dimension();
    std::vector<PolynomialHP> grad;
    grad.reserve(dim);

    for (std::size_t i = 0; i < dim; ++i) {
        grad.push_back(differentiateAxis(p, i));
    }

    return grad;
}

PolynomialHP DifferentiationHP::derivativeFromDouble(const Polynomial& p, std::size_t axis, unsigned int order)
{
    // Convert to high precision
    PolynomialHP p_hp(p);

    // Differentiate in high precision
    return derivative(p_hp, axis, order);
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

