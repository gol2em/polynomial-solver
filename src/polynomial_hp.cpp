#ifdef ENABLE_HIGH_PRECISION

#include "polynomial_hp.h"
#include "precision_conversion.h"
#include <stdexcept>

namespace {

// Helper functions for power-to-Bernstein conversion in HP

/**
 * @brief Compute strides for tensor-product layout
 */
void compute_strides_hp(const std::vector<unsigned int>& degrees,
                        std::vector<std::size_t>& strides) {
    const std::size_t dim = degrees.size();
    strides.resize(dim);
    if (dim == 0) return;

    strides[dim - 1] = 1;
    for (std::size_t i = dim - 1; i > 0; --i) {
        strides[i - 1] = strides[i] * static_cast<std::size_t>(degrees[i] + 1);
    }
}

/**
 * @brief Flatten multi-index to linear index
 */
std::size_t flatten_index_hp(const std::vector<unsigned int>& multi_index,
                              const std::vector<std::size_t>& strides) {
    std::size_t idx = 0;
    for (std::size_t i = 0; i < multi_index.size(); ++i) {
        idx += static_cast<std::size_t>(multi_index[i]) * strides[i];
    }
    return idx;
}

/**
 * @brief Increment multi-index, skipping one axis
 */
bool increment_multi_except_axis_hp(std::vector<unsigned int>& multi_index,
                                    const std::vector<unsigned int>& degrees,
                                    std::size_t skip_axis) {
    const std::size_t dim = degrees.size();
    for (std::size_t i = dim; i > 0; --i) {
        const std::size_t axis = i - 1;
        if (axis == skip_axis) continue;

        if (multi_index[axis] < degrees[axis]) {
            ++multi_index[axis];
            return true;
        }
        multi_index[axis] = 0;
    }
    return false;
}

/**
 * @brief Convert 1D power basis to Bernstein basis (high-precision)
 *
 * Formula: b_k = sum_{i=0}^k a_i * C(k,i) / C(n,i)
 * Uses stable downward recurrence for ratios.
 */
void power_to_bernstein_1d_hp(unsigned int degree,
                               const std::vector<polynomial_solver::mpreal>& power_coeffs_1d,
                               std::vector<polynomial_solver::mpreal>& bernstein_coeffs_1d) {
    using polynomial_solver::mpreal;

    const std::size_t len = static_cast<std::size_t>(degree + 1u);
    bernstein_coeffs_1d.assign(len, mpreal(0));

    if (len == 0u) return;

    const int n = static_cast<int>(degree);

    // Accumulate contributions from each power-basis coefficient a_i
    for (int i = 0; i <= n; ++i) {
        const mpreal& a_i = power_coeffs_1d[static_cast<std::size_t>(i)];

        // For fixed i, start from r_{n,i} = C(n,i)/C(n,i) = 1
        // Move downward in k using: r_{k,i} = r_{k+1,i} * (k + 1 - i) / (k + 1)
        mpreal r = mpreal(1);

        // Contribution to k = n
        bernstein_coeffs_1d[static_cast<std::size_t>(n)] += a_i * r;

        // Contributions to k = n-1, ..., i
        for (int k = n - 1; k >= i; --k) {
            r *= mpreal(k + 1 - i) / mpreal(k + 1);
            bernstein_coeffs_1d[static_cast<std::size_t>(k)] += a_i * r;
        }
    }
}

} // anonymous namespace

namespace polynomial_solver {

// Helper: De Casteljau evaluation for 1D high-precision
static mpreal evaluate1D_HP(const std::vector<mpreal>& bernstein_coeffs, const mpreal& t) {
    const std::size_t n = bernstein_coeffs.size();
    if (n == 0) {
        return mpreal(0);
    }

    // Standard De Casteljau algorithm (univariate)
    std::vector<mpreal> tmp(bernstein_coeffs.begin(), bernstein_coeffs.end());

    mpreal one_minus_t = mpreal(1) - t;
    for (std::size_t r = 1; r < n; ++r) {
        const std::size_t upper = n - r;
        for (std::size_t i = 0; i < upper; ++i) {
            tmp[i] = one_minus_t * tmp[i] + t * tmp[i + 1];
        }
    }

    return tmp[0];
}

// Helper: De Casteljau evaluation for tensor-product high-precision
static mpreal evaluateTensorProduct_HP(
    const std::vector<unsigned int>& degrees,
    const std::vector<mpreal>& bernstein_coeffs,
    const std::vector<mpreal>& parameters)
{
    const std::size_t dim = degrees.size();

    // Univariate case
    if (dim == 1 && parameters.size() == 1) {
        const unsigned int deg = degrees[0];
        const std::size_t count = static_cast<std::size_t>(deg + 1);
        if (bernstein_coeffs.size() != count) {
            return mpreal(0);
        }
        return evaluate1D_HP(bernstein_coeffs, parameters[0]);
    }

    // Bivariate tensor-product case
    if (dim == 2 && parameters.size() == 2) {
        const unsigned int deg_x = degrees[0];
        const unsigned int deg_y = degrees[1];
        const std::size_t nx = static_cast<std::size_t>(deg_x + 1);
        const std::size_t ny = static_cast<std::size_t>(deg_y + 1);

        if (bernstein_coeffs.size() != nx * ny) {
            return mpreal(0);
        }

        const mpreal& t_x = parameters[0];
        const mpreal& t_y = parameters[1];

        // First, evaluate along y-direction for each x-slice
        std::vector<mpreal> y_vals(nx);
        for (std::size_t i = 0; i < nx; ++i) {
            std::vector<mpreal> y_coeffs(ny);
            for (std::size_t j = 0; j < ny; ++j) {
                y_coeffs[j] = bernstein_coeffs[i * ny + j];
            }
            y_vals[i] = evaluate1D_HP(y_coeffs, t_y);
        }

        // Then evaluate along x-direction
        return evaluate1D_HP(y_vals, t_x);
    }

    // Higher dimensions not implemented
    throw std::runtime_error("PolynomialHP: evaluateTensorProduct only supports 1D and 2D");
}

// PolynomialHP implementation

PolynomialHP::PolynomialHP()
    : degrees_(), bernstein_coeffs_()
{
}

PolynomialHP::PolynomialHP(const std::vector<unsigned int>& degrees,
                           const std::vector<mpreal>& bernstein_coeffs)
    : degrees_(degrees), bernstein_coeffs_(bernstein_coeffs)
{
}

PolynomialHP::PolynomialHP(const Polynomial& poly)
    : degrees_(poly.degrees())
{
    // Convert double coefficients to high-precision
    const std::vector<double>& double_coeffs = poly.bernsteinCoefficients();
    bernstein_coeffs_.reserve(double_coeffs.size());
    for (double coeff : double_coeffs) {
        bernstein_coeffs_.push_back(toHighPrecision(coeff));
    }
}

PolynomialHP::~PolynomialHP() {
}

std::size_t PolynomialHP::dimension() const {
    return degrees_.size();
}

const std::vector<unsigned int>& PolynomialHP::degrees() const {
    return degrees_;
}

std::size_t PolynomialHP::coefficientCount() const {
    return bernstein_coeffs_.size();
}

const std::vector<mpreal>& PolynomialHP::bernsteinCoefficients() const {
    return bernstein_coeffs_;
}

mpreal PolynomialHP::evaluate(const std::vector<mpreal>& parameters) const {
    return evaluateTensorProduct_HP(degrees_, bernstein_coeffs_, parameters);
}

mpreal PolynomialHP::evaluate(const mpreal& t) const {
    std::vector<mpreal> params(1, t);
    return evaluate(params);
}

bool PolynomialHP::empty() const {
    return bernstein_coeffs_.empty();
}

// Conversion utility

PolynomialHP convertToHighPrecision(const Polynomial& poly) {
    return PolynomialHP(poly);
}

PolynomialHP fromPowerHP(const std::vector<unsigned int>& degrees,
                         const std::vector<mpreal>& power_coeffs) {
    const std::size_t dim = degrees.size();

    // Compute expected number of coefficients
    std::size_t count = 1u;
    for (std::size_t i = 0; i < degrees.size(); ++i) {
        count *= static_cast<std::size_t>(degrees[i] + 1u);
    }

    if (count == 0u) {
        return PolynomialHP(degrees, std::vector<mpreal>());
    }

    // Copy and resize input coefficients as needed
    std::vector<mpreal> coeffs = power_coeffs;
    if (coeffs.size() < count) {
        coeffs.resize(count, mpreal(0));
    } else if (coeffs.size() > count) {
        coeffs.resize(count);
    }

    if (dim == 0u) {
        return PolynomialHP(degrees, coeffs);
    }

    // Precompute strides for the tensor layout
    std::vector<std::size_t> strides;
    compute_strides_hp(degrees, strides);

    // Transform along each dimension separately from power basis to Bernstein
    std::vector<unsigned int> multi_index(dim, 0u);

    for (std::size_t axis = 0; axis < dim; ++axis) {
        const unsigned int deg_axis = degrees[axis];
        const std::size_t len_axis = static_cast<std::size_t>(deg_axis + 1u);

        if (len_axis <= 1u) {
            continue; // nothing to do along this dimension
        }

        std::vector<mpreal> line_in(len_axis);
        std::vector<mpreal> line_out;

        std::fill(multi_index.begin(), multi_index.end(), 0u);
        bool first = true;

        while (first || increment_multi_except_axis_hp(multi_index, degrees, axis)) {
            first = false;

            // Gather a 1D slice along the current axis
            for (std::size_t k = 0; k < len_axis; ++k) {
                multi_index[axis] = static_cast<unsigned int>(k);
                const std::size_t idx = flatten_index_hp(multi_index, strides);
                line_in[k] = coeffs[idx];
            }

            // Convert this slice from power basis to Bernstein basis
            power_to_bernstein_1d_hp(deg_axis, line_in, line_out);

            // Scatter converted coefficients back
            for (std::size_t k = 0; k < len_axis; ++k) {
                multi_index[axis] = static_cast<unsigned int>(k);
                const std::size_t idx = flatten_index_hp(multi_index, strides);
                coeffs[idx] = line_out[k];
            }
        }
    }

    // Raise degree to at least 1 in each dimension (in Bernstein basis)
    std::vector<unsigned int> adjusted_degrees = degrees;
    bool needs_degree_raising = false;
    for (std::size_t i = 0; i < dim; ++i) {
        if (adjusted_degrees[i] == 0u) {
            adjusted_degrees[i] = 1u;
            needs_degree_raising = true;
        }
    }

    if (needs_degree_raising) {
        // Compute new coefficient count
        std::size_t new_count = 1u;
        for (std::size_t i = 0; i < dim; ++i) {
            new_count *= static_cast<std::size_t>(adjusted_degrees[i] + 1u);
        }

        std::vector<mpreal> raised_coeffs(new_count, mpreal(0));

        // Copy coefficients, duplicating along dimensions where degree was raised
        std::vector<std::size_t> old_strides, new_strides;
        compute_strides_hp(degrees, old_strides);
        compute_strides_hp(adjusted_degrees, new_strides);

        std::vector<unsigned int> old_idx(dim, 0u);
        std::vector<unsigned int> new_idx(dim, 0u);

        for (std::size_t i = 0; i < count; ++i) {
            // Compute old multi-index
            std::size_t temp = i;
            for (std::size_t d = 0; d < dim; ++d) {
                old_idx[d] = static_cast<unsigned int>(temp / old_strides[d]);
                temp %= old_strides[d];
            }

            // Map to new multi-index (same values, but duplicated for degree-0 dimensions)
            for (std::size_t d = 0; d < dim; ++d) {
                new_idx[d] = old_idx[d];
            }

            // Write to all positions in raised dimensions
            std::vector<unsigned int> write_idx = new_idx;
            do {
                const std::size_t new_i = flatten_index_hp(write_idx, new_strides);
                raised_coeffs[new_i] = coeffs[i];

                // Increment along raised dimensions
                bool carry = true;
                for (std::size_t d = 0; d < dim && carry; ++d) {
                    if (degrees[d] == 0u && adjusted_degrees[d] == 1u) {
                        if (write_idx[d] == 0u) {
                            write_idx[d] = 1u;
                            carry = false;
                        } else {
                            write_idx[d] = new_idx[d];
                        }
                    }
                }
                if (carry) break;
            } while (true);
        }

        return PolynomialHP(adjusted_degrees, raised_coeffs);
    }

    return PolynomialHP(degrees, coeffs);
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

