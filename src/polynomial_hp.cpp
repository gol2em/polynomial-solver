#include "config.h"

#ifdef ENABLE_HIGH_PRECISION

#include "hp/polynomial_hp.h"
#include "hp/precision_conversion.h"
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

/**
 * @brief Horner's method for 1D power basis evaluation (high-precision)
 */
polynomial_solver::mpreal horner_eval_1d_hp(const std::vector<polynomial_solver::mpreal>& power_coeffs_1d,
                                            const polynomial_solver::mpreal& x) {
    using polynomial_solver::mpreal;

    if (power_coeffs_1d.empty()) {
        return mpreal(0);
    }
    const int n = static_cast<int>(power_coeffs_1d.size()) - 1;
    mpreal result = power_coeffs_1d[n];
    for (int i = n - 1; i >= 0; --i) {
        result = result * x + power_coeffs_1d[i];
    }
    return result;
}

/**
 * @brief Horner's method for tensor-product power basis (high-precision)
 */
polynomial_solver::mpreal horner_eval_tensor_hp(const std::vector<unsigned int>& degrees,
                                                const std::vector<polynomial_solver::mpreal>& power_coeffs,
                                                const std::vector<polynomial_solver::mpreal>& parameters) {
    using polynomial_solver::mpreal;

    const std::size_t dim = degrees.size();

    if (dim == 0) {
        return power_coeffs.empty() ? mpreal(0) : power_coeffs[0];
    }

    // Univariate case
    if (dim == 1) {
        return horner_eval_1d_hp(power_coeffs, parameters[0]);
    }

    // Multivariate: evaluate dimension by dimension
    // Start from the last dimension and work backwards
    std::vector<std::size_t> strides;
    compute_strides_hp(degrees, strides);

    std::vector<mpreal> current_coeffs = power_coeffs;

    for (std::size_t d = dim; d > 0; --d) {
        const std::size_t axis = d - 1;
        const unsigned int deg = degrees[axis];
        const std::size_t len = static_cast<std::size_t>(deg + 1);
        const mpreal& x = parameters[axis];

        // Number of 1D slices along this axis
        std::size_t num_slices = 1;
        for (std::size_t i = 0; i < axis; ++i) {
            num_slices *= static_cast<std::size_t>(degrees[i] + 1);
        }

        std::vector<mpreal> next_coeffs(num_slices);

        for (std::size_t slice = 0; slice < num_slices; ++slice) {
            // Extract 1D polynomial along this axis
            std::vector<mpreal> line(len);
            for (std::size_t k = 0; k < len; ++k) {
                line[k] = current_coeffs[slice * len + k];
            }
            // Evaluate using Horner
            next_coeffs[slice] = horner_eval_1d_hp(line, x);
        }

        current_coeffs = next_coeffs;
    }

    return current_coeffs[0];
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
    : degrees_(),
      bernstein_coeffs_(),
      power_coeffs_(),
      primary_rep_(PolynomialRepresentation::BERNSTEIN),
      bernstein_valid_(false),
      power_valid_(false)
{
}

PolynomialHP::PolynomialHP(const std::vector<unsigned int>& degrees,
                           const std::vector<mpreal>& bernstein_coeffs)
    : degrees_(degrees),
      bernstein_coeffs_(bernstein_coeffs),
      power_coeffs_(),
      primary_rep_(PolynomialRepresentation::BERNSTEIN),
      bernstein_valid_(true),
      power_valid_(false)
{
}

PolynomialHP::PolynomialHP(const Polynomial& poly)
    : degrees_(poly.degrees()),
      primary_rep_(poly.primaryRepresentation()),
      bernstein_valid_(false),
      power_valid_(false)
{
    // Convert double coefficients to high-precision
    // Use high-precision conversion for power-to-Bernstein if polynomial has power as primary

    if (poly.hasPowerCoefficients()) {
        // Convert power coefficients to HP
        const std::vector<double>& double_power = poly.powerCoefficients();
        power_coeffs_.reserve(double_power.size());
        for (double coeff : double_power) {
            power_coeffs_.push_back(toHighPrecision(coeff));
        }
        power_valid_ = true;

        // If power is primary, use HP conversion to Bernstein (more accurate)
        if (primary_rep_ == PolynomialRepresentation::POWER) {
            // Bernstein will be computed lazily using HP conversion
            bernstein_valid_ = false;
        } else {
            // Also convert Bernstein coefficients
            const std::vector<double>& double_bern = poly.bernsteinCoefficients();
            bernstein_coeffs_.reserve(double_bern.size());
            for (double coeff : double_bern) {
                bernstein_coeffs_.push_back(toHighPrecision(coeff));
            }
            bernstein_valid_ = true;
        }
    } else {
        // Only Bernstein available
        const std::vector<double>& double_bern = poly.bernsteinCoefficients();
        bernstein_coeffs_.reserve(double_bern.size());
        for (double coeff : double_bern) {
            bernstein_coeffs_.push_back(toHighPrecision(coeff));
        }
        bernstein_valid_ = true;
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
    std::size_t count = 1;
    for (unsigned int deg : degrees_) {
        count *= static_cast<std::size_t>(deg + 1);
    }
    return count;
}

const std::vector<mpreal>& PolynomialHP::bernsteinCoefficients() const {
    if (!bernstein_valid_) {
        std::cerr << "\n========================================\n";
        std::cerr << "ERROR: Implicit Power->Bernstein conversion detected!\n";
        std::cerr << "========================================\n";
        std::cerr << "PolynomialHP::bernsteinCoefficients() called on polynomial with invalid Bernstein representation.\n";
        std::cerr << "This triggers implicit conversion which may introduce errors.\n";
        std::cerr << "\nTo fix this:\n";
        std::cerr << "1. If you need Bernstein coefficients, call convertPowerToBernstein() explicitly first\n";
        std::cerr << "2. If you're differentiating, use power-basis differentiation instead\n";
        std::cerr << "3. If you're evaluating, ensure primary representation matches your needs\n";
        std::cerr << "\nPrimary representation: " << (primary_rep_ == PolynomialRepresentation::POWER ? "POWER" : "BERNSTEIN") << "\n";
        std::cerr << "Power valid: " << (power_valid_ ? "YES" : "NO") << "\n";
        std::cerr << "Bernstein valid: " << (bernstein_valid_ ? "YES" : "NO") << "\n";
        std::cerr << "\nStack trace hint: Set a breakpoint at polynomial_hp.cpp:330 to see call stack\n";
        std::cerr << "========================================\n";
        std::abort();  // Use abort() instead of exit() to get better stack trace
    }
    return bernstein_coeffs_;
}

const std::vector<mpreal>& PolynomialHP::powerCoefficients() const {
    if (!power_valid_) {
        std::cerr << "\n========================================\n";
        std::cerr << "ERROR: Implicit Bernstein->Power conversion detected!\n";
        std::cerr << "========================================\n";
        std::cerr << "PolynomialHP::powerCoefficients() called on polynomial with invalid Power representation.\n";
        std::cerr << "This triggers implicit conversion which may introduce errors.\n";
        std::cerr << "\nTo fix this:\n";
        std::cerr << "1. If you need power coefficients, call convertBernsteinToPower() explicitly first\n";
        std::cerr << "2. Ensure the polynomial was created with the correct primary representation\n";
        std::cerr << "\nPrimary representation: " << (primary_rep_ == PolynomialRepresentation::POWER ? "POWER" : "BERNSTEIN") << "\n";
        std::cerr << "Power valid: " << (power_valid_ ? "YES" : "NO") << "\n";
        std::cerr << "Bernstein valid: " << (bernstein_valid_ ? "YES" : "NO") << "\n";
        std::cerr << "========================================\n";
        std::exit(EXIT_FAILURE);
    }
    return power_coeffs_;
}

PolynomialRepresentation PolynomialHP::primaryRepresentation() const {
    return primary_rep_;
}

bool PolynomialHP::hasPowerCoefficients() const {
    return power_valid_;
}

bool PolynomialHP::hasBernsteinCoefficients() const {
    return bernstein_valid_;
}

void PolynomialHP::ensureBernsteinPrimary() {
    if (primary_rep_ == PolynomialRepresentation::BERNSTEIN) {
        return;
    }
    if (!bernstein_valid_) {
        convertPowerToBernstein();
    }
    primary_rep_ = PolynomialRepresentation::BERNSTEIN;
}

void PolynomialHP::ensurePowerPrimary() {
    if (primary_rep_ == PolynomialRepresentation::POWER) {
        return;
    }
    if (!power_valid_) {
        convertBernsteinToPower();
    }
    primary_rep_ = PolynomialRepresentation::POWER;
}

mpreal PolynomialHP::evaluate(const std::vector<mpreal>& parameters) const {
    // Choose evaluation method based on PRIMARY representation
    if (primary_rep_ == PolynomialRepresentation::POWER) {
        if (!power_valid_) {
            convertBernsteinToPower();
        }
        return horner_eval_tensor_hp(degrees_, power_coeffs_, parameters);
    } else {
        if (!bernstein_valid_) {
            convertPowerToBernstein();
        }
        return evaluateTensorProduct_HP(degrees_, bernstein_coeffs_, parameters);
    }
}

mpreal PolynomialHP::evaluate(const mpreal& t) const {
    std::vector<mpreal> params(1, t);
    return evaluate(params);
}

bool PolynomialHP::empty() const {
    return !bernstein_valid_ && !power_valid_;
}

void PolynomialHP::convertPowerToBernstein() const {
    // Convert power to Bernstein in high precision
    const std::size_t dim = degrees_.size();

    if (dim == 0u) {
        bernstein_coeffs_ = power_coeffs_;
        bernstein_valid_ = true;
        return;
    }

    // Copy power coefficients to working array
    std::vector<mpreal> coeffs = power_coeffs_;

    // Precompute strides for the tensor layout
    std::vector<std::size_t> strides;
    compute_strides_hp(degrees_, strides);

    // Transform along each dimension separately from power basis to Bernstein
    std::vector<unsigned int> multi_index(dim, 0u);

    for (std::size_t axis = 0; axis < dim; ++axis) {
        const unsigned int deg_axis = degrees_[axis];
        const std::size_t len_axis = static_cast<std::size_t>(deg_axis + 1u);

        if (len_axis <= 1u) {
            continue; // nothing to do along this dimension
        }

        std::vector<mpreal> line_in(len_axis);
        std::vector<mpreal> line_out;

        std::fill(multi_index.begin(), multi_index.end(), 0u);
        bool first = true;

        while (first || increment_multi_except_axis_hp(multi_index, degrees_, axis)) {
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

    bernstein_coeffs_ = coeffs;
    bernstein_valid_ = true;
}

void PolynomialHP::convertBernsteinToPower() const {
    // TODO: Implement Bernstein-to-power conversion in HP
    // For now, just copy (this is a placeholder)
    // This conversion is rarely needed in practice
    power_coeffs_ = bernstein_coeffs_;
    power_valid_ = true;
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
        PolynomialHP result;
        result.degrees_ = degrees;
        result.primary_rep_ = PolynomialRepresentation::POWER;
        result.power_valid_ = true;
        result.bernstein_valid_ = false;
        return result;
    }

    // Copy and resize input coefficients as needed
    std::vector<mpreal> coeffs = power_coeffs;
    if (coeffs.size() < count) {
        coeffs.resize(count, mpreal(0));
    } else if (coeffs.size() > count) {
        coeffs.resize(count);
    }

    // NEW DESIGN: Store power coefficients as primary, Bernstein computed lazily
    // This avoids double precision limitations and uses HP conversion when needed
    PolynomialHP result;
    result.degrees_ = degrees;
    result.power_coeffs_ = coeffs;
    result.primary_rep_ = PolynomialRepresentation::POWER;
    result.power_valid_ = true;
    result.bernstein_valid_ = false;  // Will be computed lazily when needed

    return result;
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

