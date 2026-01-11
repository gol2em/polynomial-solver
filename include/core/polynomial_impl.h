#ifndef POLYNOMIAL_IMPL_H
#define POLYNOMIAL_IMPL_H

/**
 * @file detail/polynomial_impl.h
 * @brief Template implementations for polynomial algorithms
 *
 * This header provides templated implementations of core polynomial algorithms
 * that work with any scalar type (double, mpreal, etc.). It is part of Tier 3
 * (template-based) implementation.
 *
 * Design principle: All algorithms are templated on Scalar type. The public
 * Polynomial class (double) and PolynomialBase<T> use these implementations.
 */

#include <vector>
#include <cstddef>
#include <algorithm>
#include <stdexcept>

namespace polynomial_solver {
namespace detail {

//=============================================================================
// Tensor Layout Helpers (templated on nothing, just utilities)
//=============================================================================

/**
 * @brief Compute strides for tensor-product layout (last dimension fastest)
 */
inline void compute_strides(const std::vector<unsigned int>& degrees,
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

/**
 * @brief Flatten multi-index using precomputed strides
 */
inline std::size_t flatten_index(const std::vector<unsigned int>& multi_index,
                                  const std::vector<std::size_t>& strides)
{
    std::size_t idx = 0u;
    const std::size_t dim = multi_index.size();
    for (std::size_t k = 0; k < dim; ++k) {
        idx += static_cast<std::size_t>(multi_index[k]) * strides[k];
    }
    return idx;
}

/**
 * @brief Increment multi-index over all dimensions except given axis
 * @return false when all combinations exhausted
 */
inline bool increment_multi_except_axis(std::vector<unsigned int>& index,
                                         const std::vector<unsigned int>& degrees,
                                         std::size_t axis)
{
    const std::size_t dim = degrees.size();

    for (std::size_t d = dim; d-- > 0;) {
        if (d == axis) {
            continue;
        }
        if (index[d] < degrees[d]) {
            ++index[d];
            // Reset all dimensions greater than d (except axis) to zero.
            for (std::size_t e = d + 1; e < dim; ++e) {
                if (e == axis) {
                    continue;
                }
                index[e] = 0u;
            }
            return true;
        }
    }

    return false;
}

/**
 * @brief Increment multi-index over all dimensions
 * @return false when all combinations exhausted
 */
inline bool increment_multi_index(std::vector<unsigned int>& index,
                                   const std::vector<unsigned int>& degrees)
{
    const std::size_t dim = degrees.size();

    for (std::size_t d = dim; d-- > 0;) {
        if (index[d] < degrees[d]) {
            ++index[d];
            // Reset all dimensions greater than d to zero.
            for (std::size_t e = d + 1; e < dim; ++e) {
                index[e] = 0u;
            }
            return true;
        }
    }

    return false;
}

//=============================================================================
// Horner Evaluation (Power Basis)
//=============================================================================

/**
 * @brief Horner's method for 1D power basis evaluation
 *
 * p(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n
 * Horner's form: p(x) = a_0 + x*(a_1 + x*(a_2 + ... + x*a_n))
 *
 * @tparam Scalar Coefficient/evaluation type (double, mpreal, etc.)
 * @param power_coeffs_1d Power basis coefficients [a_0, a_1, ..., a_n]
 * @param x Evaluation point
 * @return p(x)
 */
template<typename Scalar>
Scalar horner_eval_1d(const std::vector<Scalar>& power_coeffs_1d, const Scalar& x)
{
    if (power_coeffs_1d.empty()) {
        return Scalar(0);
    }

    // Start from highest degree coefficient
    const int n = static_cast<int>(power_coeffs_1d.size()) - 1;
    Scalar result = power_coeffs_1d[n];

    // Work backwards through coefficients
    for (int i = n - 1; i >= 0; --i) {
        result = result * x + power_coeffs_1d[i];
    }

    return result;
}

/**
 * @brief Horner's method for tensor-product power basis (multivariate)
 *
 * Evaluates dimension by dimension, similar to De Casteljau tensor product.
 *
 * @tparam Scalar Coefficient/evaluation type
 * @param degrees Degree in each dimension
 * @param power_coeffs Tensor-product power coefficients
 * @param parameters Evaluation point (size = dim)
 * @return p(parameters)
 */
template<typename Scalar>
Scalar horner_eval_tensor(const std::vector<unsigned int>& degrees,
                          const std::vector<Scalar>& power_coeffs,
                          const std::vector<Scalar>& parameters)
{
    const std::size_t dim = degrees.size();

    if (dim == 0) {
        return power_coeffs.empty() ? Scalar(0) : power_coeffs[0];
    }

    // Univariate case: directly use Horner's method
    if (dim == 1u && parameters.size() == 1u) {
        return horner_eval_1d(power_coeffs, parameters[0]);
    }

    // Bivariate case: evaluate dimension by dimension
    if (dim == 2u && parameters.size() == 2u) {
        const std::size_t nx = static_cast<std::size_t>(degrees[0] + 1u);
        const std::size_t ny = static_cast<std::size_t>(degrees[1] + 1u);

        const Scalar& x = parameters[0];
        const Scalar& y = parameters[1];

        // First evaluate along the second dimension (y) for each fixed x index
        std::vector<Scalar> row_values(nx, Scalar(0));

        for (std::size_t ix = 0; ix < nx; ++ix) {
            const std::size_t offset = ix * ny;
            std::vector<Scalar> row(power_coeffs.begin() + static_cast<std::ptrdiff_t>(offset),
                                    power_coeffs.begin() + static_cast<std::ptrdiff_t>(offset + ny));
            row_values[ix] = horner_eval_1d(row, y);
        }

        // Then evaluate the resulting 1D polynomial in x
        return horner_eval_1d(row_values, x);
    }

    // Higher dimensions: iterative evaluation
    std::vector<std::size_t> strides;
    compute_strides(degrees, strides);

    std::vector<Scalar> current_coeffs = power_coeffs;

    for (std::size_t d = dim; d > 0; --d) {
        const std::size_t axis = d - 1;
        const unsigned int deg = degrees[axis];
        const std::size_t len = static_cast<std::size_t>(deg + 1);
        const Scalar& x = parameters[axis];

        // Number of 1D slices along this axis
        std::size_t num_slices = 1;
        for (std::size_t i = 0; i < axis; ++i) {
            num_slices *= static_cast<std::size_t>(degrees[i] + 1);
        }

        std::vector<Scalar> next_coeffs(num_slices);

        for (std::size_t slice = 0; slice < num_slices; ++slice) {
            // Extract 1D polynomial along this axis
            std::vector<Scalar> line(len);
            for (std::size_t k = 0; k < len; ++k) {
                line[k] = current_coeffs[slice * len + k];
            }
            // Evaluate using Horner
            next_coeffs[slice] = horner_eval_1d(line, x);
        }

        current_coeffs = next_coeffs;
    }

    return current_coeffs[0];
}

//=============================================================================
// De Casteljau Evaluation (Bernstein Basis)
//=============================================================================

/**
 * @brief De Casteljau algorithm for 1D Bernstein polynomial evaluation
 *
 * @tparam Scalar Coefficient/evaluation type
 * @param bernstein_coeffs Bernstein coefficients
 * @param t Evaluation point (typically in [0,1])
 * @return p(t)
 */
template<typename Scalar>
Scalar de_casteljau_1d(const std::vector<Scalar>& bernstein_coeffs, const Scalar& t)
{
    const std::size_t n = bernstein_coeffs.size();
    if (n == 0) {
        return Scalar(0);
    }

    // Standard De Casteljau algorithm
    std::vector<Scalar> tmp(bernstein_coeffs.begin(), bernstein_coeffs.end());

    Scalar one_minus_t = Scalar(1) - t;
    for (std::size_t r = 1; r < n; ++r) {
        const std::size_t upper = n - r;
        for (std::size_t i = 0; i < upper; ++i) {
            tmp[i] = one_minus_t * tmp[i] + t * tmp[i + 1];
        }
    }

    return tmp[0];
}

/**
 * @brief De Casteljau for tensor-product Bernstein basis (multivariate)
 *
 * @tparam Scalar Coefficient/evaluation type
 * @param degrees Degree in each dimension
 * @param bernstein_coeffs Tensor-product Bernstein coefficients
 * @param parameters Evaluation point
 * @return p(parameters)
 */
template<typename Scalar>
Scalar de_casteljau_tensor(const std::vector<unsigned int>& degrees,
                            const std::vector<Scalar>& bernstein_coeffs,
                            const std::vector<Scalar>& parameters)
{
    const std::size_t dim = degrees.size();

    // Univariate case
    if (dim == 1 && parameters.size() == 1) {
        const unsigned int deg = degrees[0];
        const std::size_t count = static_cast<std::size_t>(deg + 1);
        if (bernstein_coeffs.size() != count) {
            return Scalar(0);
        }
        return de_casteljau_1d(bernstein_coeffs, parameters[0]);
    }

    // Bivariate tensor-product case
    if (dim == 2 && parameters.size() == 2) {
        const unsigned int deg_x = degrees[0];
        const unsigned int deg_y = degrees[1];
        const std::size_t nx = static_cast<std::size_t>(deg_x + 1);
        const std::size_t ny = static_cast<std::size_t>(deg_y + 1);

        if (bernstein_coeffs.size() != nx * ny) {
            return Scalar(0);
        }

        const Scalar& t_x = parameters[0];
        const Scalar& t_y = parameters[1];

        // First, evaluate along y-direction for each x-slice
        std::vector<Scalar> y_vals(nx);
        for (std::size_t i = 0; i < nx; ++i) {
            std::vector<Scalar> y_coeffs(ny);
            for (std::size_t j = 0; j < ny; ++j) {
                y_coeffs[j] = bernstein_coeffs[i * ny + j];
            }
            y_vals[i] = de_casteljau_1d(y_coeffs, t_y);
        }

        // Then evaluate along x-direction
        return de_casteljau_1d(y_vals, t_x);
    }

    // General n-dimensional case: recursive dimension reduction
    // Evaluate dimension by dimension from last to first
    if (dim >= 3 && parameters.size() == dim) {
        // Get sizes for each dimension
        std::vector<std::size_t> sizes(dim);
        std::size_t total = 1;
        for (std::size_t i = 0; i < dim; ++i) {
            sizes[i] = static_cast<std::size_t>(degrees[i] + 1);
            total *= sizes[i];
        }

        if (bernstein_coeffs.size() != total) {
            return Scalar(0);
        }

        // Start with all coefficients
        std::vector<Scalar> current = bernstein_coeffs;
        std::vector<std::size_t> current_sizes = sizes;

        // Reduce dimension by dimension from last to first
        for (std::size_t d = dim; d >= 1; --d) {
            std::size_t axis = d - 1;  // Current axis to reduce
            std::size_t axis_size = current_sizes[axis];
            Scalar t = parameters[axis];

            // Compute size of reduced array
            std::size_t reduced_size = 1;
            for (std::size_t i = 0; i < axis; ++i) {
                reduced_size *= current_sizes[i];
            }
            for (std::size_t i = axis + 1; i < d; ++i) {
                reduced_size *= current_sizes[i];
            }

            if (axis == 0 && d == 1) {
                // Last reduction: just evaluate 1D
                return de_casteljau_1d(current, t);
            }

            // Compute stride for axis
            std::size_t stride_after = 1;
            for (std::size_t i = axis + 1; i < d; ++i) {
                stride_after *= current_sizes[i];
            }

            std::size_t stride_before = 1;
            for (std::size_t i = 0; i < axis; ++i) {
                stride_before *= current_sizes[i];
            }

            // Reduce along this axis
            std::vector<Scalar> reduced(reduced_size);
            std::size_t out_idx = 0;

            for (std::size_t before = 0; before < stride_before; ++before) {
                for (std::size_t after = 0; after < stride_after; ++after) {
                    // Extract 1D slice
                    std::vector<Scalar> slice(axis_size);
                    for (std::size_t k = 0; k < axis_size; ++k) {
                        std::size_t idx = before * (axis_size * stride_after) + k * stride_after + after;
                        slice[k] = current[idx];
                    }
                    reduced[out_idx++] = de_casteljau_1d(slice, t);
                }
            }

            current = reduced;

            // Update sizes: remove the current axis
            std::vector<std::size_t> new_sizes;
            for (std::size_t i = 0; i < d; ++i) {
                if (i != axis) {
                    new_sizes.push_back(current_sizes[i]);
                }
            }
            current_sizes = new_sizes;
        }

        return current.empty() ? Scalar(0) : current[0];
    }

    // Fallback
    throw std::runtime_error("de_casteljau_tensor: unsupported dimension");
}

//=============================================================================
// Power-to-Bernstein Conversion
//=============================================================================

/**
 * @brief Convert 1D power basis to Bernstein basis
 *
 * Uses the formula: b_k = sum_{i=0}^k a_i * C(k,i) / C(n,i)
 * with stable downward recurrence for the binomial ratios.
 *
 * @tparam Scalar Coefficient type
 * @param degree Polynomial degree
 * @param power_coeffs_1d Power coefficients [a_0, a_1, ..., a_n]
 * @param bernstein_coeffs_1d Output Bernstein coefficients [b_0, b_1, ..., b_n]
 */
template<typename Scalar>
void power_to_bernstein_1d(unsigned int degree,
                            const std::vector<Scalar>& power_coeffs_1d,
                            std::vector<Scalar>& bernstein_coeffs_1d)
{
    const std::size_t len = static_cast<std::size_t>(degree + 1u);
    bernstein_coeffs_1d.assign(len, Scalar(0));

    if (len == 0u) {
        return;
    }

    const int n = static_cast<int>(degree);

    // Accumulate contributions from each power-basis coefficient a_i.
    for (int i = 0; i <= n; ++i) {
        const Scalar& a_i = power_coeffs_1d[static_cast<std::size_t>(i)];

        // For fixed i, start from r_{n,i} = C(n,i)/C(n,i) = 1 and
        // move downward in k using
        //   r_{k,i} = r_{k+1,i} * (k + 1 - i) / (k + 1).
        Scalar r = Scalar(1);

        // Contribution to k = n.
        bernstein_coeffs_1d[static_cast<std::size_t>(n)] += a_i * r;

        // Contributions to k = n-1, ..., i.
        for (int k = n - 1; k >= i; --k) {
            r = r * Scalar(k + 1 - i) / Scalar(k + 1);
            bernstein_coeffs_1d[static_cast<std::size_t>(k)] += a_i * r;
        }
    }
}

/**
 * @brief Convert multivariate tensor-product power basis to Bernstein basis
 *
 * Applies 1D conversion along each axis sequentially.
 *
 * @tparam Scalar Coefficient type
 * @param degrees Degrees in each dimension
 * @param power_coeffs Input power coefficients
 * @param bernstein_coeffs Output Bernstein coefficients
 */
template<typename Scalar>
void power_to_bernstein_tensor(const std::vector<unsigned int>& degrees,
                                const std::vector<Scalar>& power_coeffs,
                                std::vector<Scalar>& bernstein_coeffs)
{
    const std::size_t dim = degrees.size();

    if (dim == 0u) {
        bernstein_coeffs = power_coeffs;
        return;
    }

    // Copy power coefficients to working array
    std::vector<Scalar> coeffs = power_coeffs;

    // Precompute strides for the tensor layout
    std::vector<std::size_t> strides;
    compute_strides(degrees, strides);

    // Transform along each dimension separately
    std::vector<unsigned int> multi_index(dim, 0u);

    for (std::size_t axis = 0; axis < dim; ++axis) {
        const unsigned int deg_axis = degrees[axis];
        const std::size_t len_axis = static_cast<std::size_t>(deg_axis + 1u);

        if (len_axis <= 1u) {
            continue; // nothing to do along this dimension
        }

        std::vector<Scalar> line_in(len_axis);
        std::vector<Scalar> line_out;

        std::fill(multi_index.begin(), multi_index.end(), 0u);
        bool first = true;

        while (first || increment_multi_except_axis(multi_index, degrees, axis)) {
            first = false;

            // Gather a 1D slice along the current axis
            for (std::size_t k = 0; k < len_axis; ++k) {
                multi_index[axis] = static_cast<unsigned int>(k);
                const std::size_t idx = flatten_index(multi_index, strides);
                line_in[k] = coeffs[idx];
            }

            // Convert this slice from power basis to Bernstein basis
            power_to_bernstein_1d(deg_axis, line_in, line_out);

            // Scatter converted coefficients back
            for (std::size_t k = 0; k < len_axis; ++k) {
                multi_index[axis] = static_cast<unsigned int>(k);
                const std::size_t idx = flatten_index(multi_index, strides);
                coeffs[idx] = line_out[k];
            }
        }
    }

    bernstein_coeffs = coeffs;
}

//=============================================================================
// De Casteljau Subdivision (for restrictedToInterval)
//=============================================================================

/**
 * @brief Subdivide 1D Bernstein polynomial at parameter t
 *
 * Given Bernstein coefficients on [0,1], compute coefficients for [0,t] and [t,1].
 *
 * @tparam Scalar Coefficient type
 * @param coeffs Input Bernstein coefficients
 * @param t Subdivision parameter
 * @param left Output: coefficients for [0,t]
 * @param right Output: coefficients for [t,1]
 */
template<typename Scalar>
void subdivide_1d(const std::vector<Scalar>& coeffs, const Scalar& t,
                  std::vector<Scalar>& left, std::vector<Scalar>& right)
{
    const std::size_t n = coeffs.size();
    if (n == 0) {
        left.clear();
        right.clear();
        return;
    }

    left.resize(n);
    right.resize(n);

    std::vector<Scalar> tmp = coeffs;
    Scalar one_minus_t = Scalar(1) - t;

    left[0] = tmp[0];
    right[n - 1] = tmp[n - 1];

    for (std::size_t r = 1; r < n; ++r) {
        const std::size_t upper = n - r;
        for (std::size_t i = 0; i < upper; ++i) {
            tmp[i] = one_minus_t * tmp[i] + t * tmp[i + 1];
        }
        left[r] = tmp[0];
        right[n - 1 - r] = tmp[upper - 1];
    }
}

//=============================================================================
// Polynomial Arithmetic Helpers
//=============================================================================

/**
 * @brief Multiply two tensor-product polynomials in power basis
 *
 * @tparam Scalar Coefficient type
 * @param degrees_a Degrees of first polynomial
 * @param coeffs_a Coefficients of first polynomial
 * @param degrees_b Degrees of second polynomial
 * @param coeffs_b Coefficients of second polynomial
 * @param degrees_result Output: degrees of result
 * @param coeffs_result Output: coefficients of result
 */
template<typename Scalar>
void multiply_power(const std::vector<unsigned int>& degrees_a,
                    const std::vector<Scalar>& coeffs_a,
                    const std::vector<unsigned int>& degrees_b,
                    const std::vector<Scalar>& coeffs_b,
                    std::vector<unsigned int>& degrees_result,
                    std::vector<Scalar>& coeffs_result)
{
    const std::size_t dim = degrees_a.size();
    if (dim == 0 || dim != degrees_b.size()) {
        degrees_result.clear();
        coeffs_result.clear();
        return;
    }

    // Result degrees = sum of degrees
    degrees_result.resize(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        degrees_result[i] = degrees_a[i] + degrees_b[i];
    }

    // Compute result size
    std::size_t result_size = 1;
    for (std::size_t i = 0; i < dim; ++i) {
        result_size *= (degrees_result[i] + 1);
    }
    coeffs_result.assign(result_size, Scalar(0));

    // Compute strides
    std::vector<std::size_t> strides_a, strides_b, strides_result;
    compute_strides(degrees_a, strides_a);
    compute_strides(degrees_b, strides_b);
    compute_strides(degrees_result, strides_result);

    // Iterate over all multi-indices of a
    std::vector<unsigned int> idx_a(dim, 0);
    do {
        std::size_t flat_a = flatten_index(idx_a, strides_a);
        const Scalar& coeff_a = coeffs_a[flat_a];

        // Iterate over all multi-indices of b
        std::vector<unsigned int> idx_b(dim, 0);
        do {
            std::size_t flat_b = flatten_index(idx_b, strides_b);
            const Scalar& coeff_b = coeffs_b[flat_b];

            // Result index = idx_a + idx_b
            std::vector<unsigned int> idx_result(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                idx_result[i] = idx_a[i] + idx_b[i];
            }
            std::size_t flat_result = flatten_index(idx_result, strides_result);
            coeffs_result[flat_result] += coeff_a * coeff_b;

        } while (increment_multi_index(idx_b, degrees_b));
    } while (increment_multi_index(idx_a, degrees_a));
}

/**
 * @brief Add two tensor-product polynomials in power basis
 *
 * @tparam Scalar Coefficient type
 */
template<typename Scalar>
void add_power(const std::vector<unsigned int>& degrees_a,
               const std::vector<Scalar>& coeffs_a,
               const std::vector<unsigned int>& degrees_b,
               const std::vector<Scalar>& coeffs_b,
               std::vector<unsigned int>& degrees_result,
               std::vector<Scalar>& coeffs_result)
{
    const std::size_t dim = degrees_a.size();
    if (dim == 0 || dim != degrees_b.size()) {
        degrees_result.clear();
        coeffs_result.clear();
        return;
    }

    // Result degrees = max of degrees
    degrees_result.resize(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        degrees_result[i] = std::max(degrees_a[i], degrees_b[i]);
    }

    // Compute result size
    std::size_t result_size = 1;
    for (std::size_t i = 0; i < dim; ++i) {
        result_size *= (degrees_result[i] + 1);
    }
    coeffs_result.assign(result_size, Scalar(0));

    // Compute strides
    std::vector<std::size_t> strides_a, strides_b, strides_result;
    compute_strides(degrees_a, strides_a);
    compute_strides(degrees_b, strides_b);
    compute_strides(degrees_result, strides_result);

    // Add coefficients from a
    std::vector<unsigned int> idx(dim, 0);
    do {
        std::size_t flat_a = flatten_index(idx, strides_a);
        std::size_t flat_result = flatten_index(idx, strides_result);
        coeffs_result[flat_result] += coeffs_a[flat_a];
    } while (increment_multi_index(idx, degrees_a));

    // Add coefficients from b
    idx.assign(dim, 0);
    do {
        std::size_t flat_b = flatten_index(idx, strides_b);
        std::size_t flat_result = flatten_index(idx, strides_result);
        coeffs_result[flat_result] += coeffs_b[flat_b];
    } while (increment_multi_index(idx, degrees_b));
}

/**
 * @brief Subtract two tensor-product polynomials in power basis (a - b)
 */
template<typename Scalar>
void subtract_power(const std::vector<unsigned int>& degrees_a,
                    const std::vector<Scalar>& coeffs_a,
                    const std::vector<unsigned int>& degrees_b,
                    const std::vector<Scalar>& coeffs_b,
                    std::vector<unsigned int>& degrees_result,
                    std::vector<Scalar>& coeffs_result)
{
    const std::size_t dim = degrees_a.size();
    if (dim == 0 || dim != degrees_b.size()) {
        degrees_result.clear();
        coeffs_result.clear();
        return;
    }

    // Result degrees = max of degrees
    degrees_result.resize(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        degrees_result[i] = std::max(degrees_a[i], degrees_b[i]);
    }

    // Compute result size
    std::size_t result_size = 1;
    for (std::size_t i = 0; i < dim; ++i) {
        result_size *= (degrees_result[i] + 1);
    }
    coeffs_result.assign(result_size, Scalar(0));

    // Compute strides
    std::vector<std::size_t> strides_a, strides_b, strides_result;
    compute_strides(degrees_a, strides_a);
    compute_strides(degrees_b, strides_b);
    compute_strides(degrees_result, strides_result);

    // Add coefficients from a
    std::vector<unsigned int> idx(dim, 0);
    do {
        std::size_t flat_a = flatten_index(idx, strides_a);
        std::size_t flat_result = flatten_index(idx, strides_result);
        coeffs_result[flat_result] += coeffs_a[flat_a];
    } while (increment_multi_index(idx, degrees_a));

    // Subtract coefficients from b
    idx.assign(dim, 0);
    do {
        std::size_t flat_b = flatten_index(idx, strides_b);
        std::size_t flat_result = flatten_index(idx, strides_result);
        coeffs_result[flat_result] -= coeffs_b[flat_b];
    } while (increment_multi_index(idx, degrees_b));
}

/**
 * @brief Scale polynomial coefficients by scalar
 */
template<typename Scalar>
void scale_coefficients(const std::vector<Scalar>& coeffs,
                        const Scalar& scalar,
                        std::vector<Scalar>& result)
{
    result.resize(coeffs.size());
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        result[i] = coeffs[i] * scalar;
    }
}

//=============================================================================
// Differentiation Helpers
//=============================================================================

/**
 * @brief Differentiate 1D Bernstein polynomial
 *
 * Uses the Bernstein derivative formula: d_i = n * (b_{i+1} - b_i)
 * The result has degree n-1 (one less than input).
 *
 * @tparam Scalar Coefficient type
 * @param degree Polynomial degree
 * @param coeffs Input Bernstein coefficients [b_0, ..., b_n]
 * @param deriv_coeffs Output derivative coefficients [d_0, ..., d_{n-1}]
 */
template<typename Scalar>
void differentiate_bernstein_1d(unsigned int degree,
                                 const std::vector<Scalar>& coeffs,
                                 std::vector<Scalar>& deriv_coeffs)
{
    if (degree == 0u) {
        // Derivative of constant is zero (keep same degree for consistency)
        deriv_coeffs.assign(1, Scalar(0));
        return;
    }

    const std::size_t n = static_cast<std::size_t>(degree);
    deriv_coeffs.resize(n);

    // d_i = n * (b_{i+1} - b_i) for i = 0, ..., n-1
    const Scalar deg_scalar = Scalar(static_cast<int>(degree));
    for (std::size_t i = 0; i < n; ++i) {
        deriv_coeffs[i] = deg_scalar * (coeffs[i + 1] - coeffs[i]);
    }
}

/**
 * @brief Differentiate 1D power basis polynomial
 *
 * Uses the power derivative formula: d/dx(a_i * x^i) = i * a_i * x^{i-1}
 * Coefficient of x^i in derivative is (i+1) * a_{i+1} from original.
 * The result has degree n-1 (one less than input).
 *
 * @tparam Scalar Coefficient type
 * @param degree Polynomial degree
 * @param coeffs Input power coefficients [a_0, ..., a_n]
 * @param deriv_coeffs Output derivative coefficients [d_0, ..., d_{n-1}]
 */
template<typename Scalar>
void differentiate_power_1d(unsigned int degree,
                             const std::vector<Scalar>& coeffs,
                             std::vector<Scalar>& deriv_coeffs)
{
    if (degree == 0u) {
        // Derivative of constant is zero (keep same degree for consistency)
        deriv_coeffs.assign(1, Scalar(0));
        return;
    }

    const std::size_t n = static_cast<std::size_t>(degree);
    deriv_coeffs.resize(n);

    // Coefficient of x^i in derivative is (i+1) * a_{i+1}
    for (std::size_t i = 0; i < n; ++i) {
        deriv_coeffs[i] = Scalar(static_cast<int>(i + 1)) * coeffs[i + 1];
    }
}

/**
 * @brief Differentiate multivariate tensor-product Bernstein polynomial along one axis
 *
 * @tparam Scalar Coefficient type
 * @param degrees Input degrees per dimension
 * @param coeffs Input Bernstein coefficients
 * @param axis Differentiation axis
 * @param new_degrees Output degrees (reduced by 1 along axis)
 * @param new_coeffs Output derivative coefficients
 */
template<typename Scalar>
void differentiate_bernstein_tensor(const std::vector<unsigned int>& degrees,
                                     const std::vector<Scalar>& coeffs,
                                     std::size_t axis,
                                     std::vector<unsigned int>& new_degrees,
                                     std::vector<Scalar>& new_coeffs)
{
    const std::size_t dim = degrees.size();

    // Invalid axis: return zero polynomial with same shape
    if (axis >= dim) {
        new_degrees = degrees;
        new_coeffs.assign(coeffs.size(), Scalar(0));
        return;
    }

    const unsigned int deg_axis = degrees[axis];

    // If degree is 0 along this axis, derivative is zero
    if (deg_axis == 0u) {
        new_degrees = degrees;
        new_coeffs.assign(coeffs.size(), Scalar(0));
        return;
    }

    // New degrees: reduce by 1 along the differentiation axis
    new_degrees = degrees;
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

    new_coeffs.assign(new_count, Scalar(0));

    // Iterate over all 1D slices along the given axis
    std::vector<unsigned int> multi_index(dim, 0u);
    std::vector<Scalar> line(len_axis);
    std::vector<Scalar> deriv_line;

    bool first = true;
    while (first || increment_multi_except_axis(multi_index, degrees, axis)) {
        first = false;

        // Gather coefficients along the current axis-line
        for (std::size_t k = 0; k < len_axis; ++k) {
            multi_index[axis] = static_cast<unsigned int>(k);
            const std::size_t idx = flatten_index(multi_index, old_strides);
            line[k] = coeffs[idx];
        }

        // Apply Bernstein derivative formula
        differentiate_bernstein_1d(deg_axis, line, deriv_line);

        // Scatter the derivative coefficients back into the tensor
        for (std::size_t k = 0; k < new_len_axis; ++k) {
            multi_index[axis] = static_cast<unsigned int>(k);
            const std::size_t idx = flatten_index(multi_index, new_strides);
            new_coeffs[idx] = deriv_line[k];
        }
    }
}

/**
 * @brief Differentiate multivariate tensor-product power polynomial along one axis
 *
 * @tparam Scalar Coefficient type
 * @param degrees Input degrees per dimension
 * @param coeffs Input power coefficients
 * @param axis Differentiation axis
 * @param new_degrees Output degrees (reduced by 1 along axis)
 * @param new_coeffs Output derivative coefficients
 */
template<typename Scalar>
void differentiate_power_tensor(const std::vector<unsigned int>& degrees,
                                 const std::vector<Scalar>& coeffs,
                                 std::size_t axis,
                                 std::vector<unsigned int>& new_degrees,
                                 std::vector<Scalar>& new_coeffs)
{
    const std::size_t dim = degrees.size();

    // Invalid axis: return zero polynomial with same shape
    if (axis >= dim) {
        new_degrees = degrees;
        new_coeffs.assign(coeffs.size(), Scalar(0));
        return;
    }

    const unsigned int deg_axis = degrees[axis];

    // If degree is 0 along this axis, derivative is zero
    if (deg_axis == 0u) {
        new_degrees = degrees;
        new_coeffs.assign(coeffs.size(), Scalar(0));
        return;
    }

    // New degrees: reduce by 1 along the differentiation axis
    new_degrees = degrees;
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

    new_coeffs.assign(new_count, Scalar(0));

    // Iterate over all 1D slices along the given axis
    std::vector<unsigned int> multi_index(dim, 0u);
    std::vector<Scalar> line(len_axis);
    std::vector<Scalar> deriv_line;

    bool first = true;
    while (first || increment_multi_except_axis(multi_index, degrees, axis)) {
        first = false;

        // Gather coefficients along the current axis-line
        for (std::size_t k = 0; k < len_axis; ++k) {
            multi_index[axis] = static_cast<unsigned int>(k);
            const std::size_t idx = flatten_index(multi_index, old_strides);
            line[k] = coeffs[idx];
        }

        // Apply power derivative formula
        differentiate_power_1d(deg_axis, line, deriv_line);

        // Scatter the derivative coefficients back into the tensor
        for (std::size_t k = 0; k < new_len_axis; ++k) {
            multi_index[axis] = static_cast<unsigned int>(k);
            const std::size_t idx = flatten_index(multi_index, new_strides);
            new_coeffs[idx] = deriv_line[k];
        }
    }
}

} // namespace detail
} // namespace polynomial_solver

#endif // POLYNOMIAL_IMPL_H
