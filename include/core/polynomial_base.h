#ifndef POLYNOMIAL_BASE_H
#define POLYNOMIAL_BASE_H

/**
 * @file polynomial_base.h
 * @brief Templated polynomial class for Tier 3 (flexible precision)
 *
 * This header provides PolynomialBase<Scalar>, a fully templated polynomial class
 * that works with any scalar type (double, mpreal, etc.). It unifies the
 * functionality of both Polynomial and PolynomialHP into a single template.
 *
 * Design (NEW):
 * - Single representation (power OR Bernstein, not both)
 * - Region tracking for both bases
 * - Explicit conversion methods (returns new polynomial)
 * - All polynomials evaluate on [0,1]^n
 *
 * Tier 3 usage:
 * - PolynomialBase<double> for standard precision
 * - PolynomialBase<mpreal> for high precision (when ENABLE_HIGH_PRECISION)
 */

#include "core/polynomial_impl.h"
#include <vector>
#include <map>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <stdexcept>

// Forward declare or include the enum from polynomial.h
// polynomial_base.h should be used WITH polynomial.h, not standalone
// This avoids duplicate enum definitions
#include "core/polynomial.h"

namespace polynomial_solver {

/**
 * @enum Basis
 * @brief Which basis the polynomial coefficients are stored in.
 */
enum class Basis {
    POWER,      ///< Power basis (monomial): sum c_i * x^i
    BERNSTEIN   ///< Bernstein basis: sum c_i * B_{i,n}(x)
};

/**
 * @class PolynomialBase
 * @brief Templated multivariate polynomial with single coefficient storage.
 *
 * NEW DESIGN: This class stores coefficients in a SINGLE representation
 * (either power or Bernstein basis). Conversion between bases returns
 * a NEW polynomial rather than modifying in-place.
 *
 * Key features:
 * - Single representation: Either power OR Bernstein, not both
 * - Region tracking: Both bases track original region for reconstruction
 * - Unified API: All operations work on any polynomial
 * - Automatic normalization: Algebraic ops normalize transparently
 *
 * @tparam Scalar The coefficient/evaluation type (double, mpreal, etc.)
 */
template<typename Scalar>
class PolynomialBase {
public:
    /**
     * @brief Default constructor creating an empty polynomial.
     */
    PolynomialBase()
        : dimension_(0u)
        , basis_(Basis::BERNSTEIN)
    {
        // Empty polynomial has no region
    }

    /**
     * @brief Construct from Bernstein coefficients.
     *
     * @param degrees  Degrees per variable (size = dimension).
     * @param bernstein_coeffs Coefficients in tensor-product Bernstein basis.
     */
    PolynomialBase(const std::vector<unsigned int>& degrees,
                   const std::vector<Scalar>& bernstein_coeffs)
        : dimension_(degrees.size())
        , degrees_(degrees)
        , coeffs_(bernstein_coeffs)
        , basis_(Basis::BERNSTEIN)
    {
        // Initialize region to [0,1]^n (default for backward compatibility)
        region_lower_.assign(dimension_, Scalar(0));
        region_upper_.assign(dimension_, Scalar(1));
    }

    /**
     * @brief Destructor.
     */
    ~PolynomialBase() = default;

    /**
     * @brief Copy constructor.
     */
    PolynomialBase(const PolynomialBase&) = default;

    /**
     * @brief Move constructor.
     */
    PolynomialBase(PolynomialBase&&) = default;

    /**
     * @brief Copy assignment.
     */
    PolynomialBase& operator=(const PolynomialBase&) = default;

    /**
     * @brief Move assignment.
     */
    PolynomialBase& operator=(PolynomialBase&&) = default;

    //=========================================================================
    // Factory Methods
    //=========================================================================

    /**
     * @brief Factory: construct from Bernstein coefficients.
     *
     * Creates polynomial in Bernstein basis on region [0,1]^n.
     */
    static PolynomialBase fromBernstein(const std::vector<unsigned int>& degrees,
                                        const std::vector<Scalar>& bernstein_coeffs)
    {
        return PolynomialBase(degrees, bernstein_coeffs);
    }

    /**
     * @brief Factory: construct from power basis coefficients.
     *
     * Creates polynomial in power basis on region [0,1]^n.
     */
    static PolynomialBase fromPower(const std::vector<unsigned int>& degrees,
                                    const std::vector<Scalar>& power_coeffs)
    {
        const std::size_t dim = degrees.size();

        // Compute expected number of coefficients
        std::size_t count = 1u;
        for (std::size_t i = 0; i < degrees.size(); ++i) {
            count *= static_cast<std::size_t>(degrees[i] + 1u);
        }

        // Copy and resize input coefficients as needed
        std::vector<Scalar> coeffs = power_coeffs;
        if (coeffs.size() < count) {
            coeffs.resize(count, Scalar(0));
        } else if (coeffs.size() > count) {
            coeffs.resize(count);
        }

        PolynomialBase result;
        result.dimension_ = dim;
        result.degrees_ = degrees;
        result.coeffs_ = coeffs;
        result.basis_ = Basis::POWER;
        // Initialize region to [0,1]^n (default for backward compatibility)
        result.region_lower_.assign(dim, Scalar(0));
        result.region_upper_.assign(dim, Scalar(1));

        return result;
    }

    //=========================================================================
    // Basic Accessors
    //=========================================================================

    /**
     * @brief Number of variables (polynomial dimension).
     */
    std::size_t dimension() const { return dimension_; }

    /**
     * @brief Degrees per variable.
     */
    const std::vector<unsigned int>& degrees() const { return degrees_; }

    /**
     * @brief Total number of coefficients.
     */
    std::size_t coefficientCount() const {
        std::size_t count = 1u;
        for (std::size_t i = 0; i < degrees_.size(); ++i) {
            count *= static_cast<std::size_t>(degrees_[i] + 1u);
        }
        return count;
    }

    /**
     * @brief Check if polynomial is empty.
     */
    bool empty() const {
        return coeffs_.empty() && dimension_ == 0;
    }

    /**
     * @brief Get current basis.
     */
    Basis basis() const { return basis_; }

    /**
     * @brief Check if polynomial is in Bernstein basis.
     */
    bool isBernstein() const { return basis_ == Basis::BERNSTEIN; }

    /**
     * @brief Check if polynomial is in power basis.
     */
    bool isPower() const { return basis_ == Basis::POWER; }

    //=========================================================================
    // Coefficient Access
    //=========================================================================

    /**
     * @brief Access underlying coefficients (regardless of basis).
     */
    const std::vector<Scalar>& coefficients() const {
        return coeffs_;
    }

    /**
     * @brief Access Bernstein coefficients.
     *
     * Requires polynomial to be in Bernstein basis.
     * Use convertToBernstein() first if needed.
     */
    const std::vector<Scalar>& bernsteinCoefficients() const {
        if (basis_ != Basis::BERNSTEIN) {
            std::cerr << "\n========================================\n";
            std::cerr << "ERROR: PolynomialBase::bernsteinCoefficients() called but basis is POWER!\n";
            std::cerr << "Call convertToBernstein() first.\n";
            std::cerr << "========================================\n";
            std::abort();
        }
        return coeffs_;
    }

    /**
     * @brief Access power coefficients.
     *
     * Requires polynomial to be in power basis.
     * Use convertToPower() first if needed.
     */
    const std::vector<Scalar>& powerCoefficients() const {
        if (basis_ != Basis::POWER) {
            std::cerr << "\n========================================\n";
            std::cerr << "ERROR: PolynomialBase::powerCoefficients() called but basis is BERNSTEIN!\n";
            std::cerr << "Call convertToPower() first.\n";
            std::cerr << "========================================\n";
            std::abort();
        }
        return coeffs_;
    }

    /**
     * @brief Query which representation is primary (for backward compatibility).
     *
     * @deprecated Use basis() instead.
     */
    PolynomialRepresentation primaryRepresentation() const {
        return (basis_ == Basis::POWER) ? PolynomialRepresentation::POWER
                                        : PolynomialRepresentation::BERNSTEIN;
    }

    /**
     * @brief Get the region bounds.
     *
     * Returns the region in original coordinates that this polynomial represents.
     * Polynomial evaluates on [0,1]^n, but represents function on [lower, upper].
     *
     * @return Pair of (lower_bounds, upper_bounds)
     */
    std::pair<std::vector<Scalar>, std::vector<Scalar>> getOriginalBox() const {
        return {region_lower_, region_upper_};
    }

    /**
     * @brief Check if power coefficients are currently available (backward compat).
     */
    bool hasPowerCoefficients() const { return basis_ == Basis::POWER; }

    /**
     * @brief Check if Bernstein coefficients are currently available (backward compat).
     */
    bool hasBernsteinCoefficients() const { return basis_ == Basis::BERNSTEIN; }

    //=========================================================================
    // Basis Conversion (NEW: returns new polynomial)
    //=========================================================================

    /**
     * @brief Convert to Bernstein basis.
     *
     * Returns a NEW polynomial in Bernstein basis.
     * Region bounds are preserved.
     */
    PolynomialBase convertToBernstein() const {
        if (basis_ == Basis::BERNSTEIN) {
            return *this;  // Already Bernstein
        }

        // Power -> Bernstein conversion
        PolynomialBase result;
        result.dimension_ = dimension_;
        result.degrees_ = degrees_;
        result.basis_ = Basis::BERNSTEIN;
        result.region_lower_ = region_lower_;
        result.region_upper_ = region_upper_;

        if (dimension_ == 0u || coeffs_.empty()) {
            result.coeffs_ = coeffs_;
            return result;
        }

        detail::power_to_bernstein_tensor(degrees_, coeffs_, result.coeffs_);
        return result;
    }

    /**
     * @brief Convert to power basis.
     *
     * Returns a NEW polynomial in power basis.
     * Region bounds are preserved.
     *
     * Note: Bernstein->Power conversion is numerically less stable.
     */
    PolynomialBase convertToPower() const {
        if (basis_ == Basis::POWER) {
            return *this;  // Already power
        }

        // Bernstein -> Power conversion
        PolynomialBase result;
        result.dimension_ = dimension_;
        result.degrees_ = degrees_;
        result.basis_ = Basis::POWER;
        result.region_lower_ = region_lower_;
        result.region_upper_ = region_upper_;

        if (dimension_ == 0u || coeffs_.empty()) {
            result.coeffs_ = coeffs_;
            return result;
        }

        // TODO: Implement proper Bernstein-to-power conversion
        // For now, placeholder
        result.coeffs_ = coeffs_;
        return result;
    }

    //=========================================================================
    // Backward compatibility: in-place conversion (deprecated)
    //=========================================================================

    /**
     * @deprecated Use convertToBernstein() instead.
     */
    void ensureBernsteinPrimary() {
        if (basis_ == Basis::BERNSTEIN) return;
        *this = convertToBernstein();
    }

    /**
     * @deprecated Use convertToPower() instead.
     */
    void ensurePowerPrimary() {
        if (basis_ == Basis::POWER) return;
        *this = convertToPower();
    }

    /**
     * @deprecated Use convertToBernstein() instead.
     */
    void convertPowerToBernstein() {
        if (basis_ == Basis::POWER) {
            *this = convertToBernstein();
        }
    }

    /**
     * @deprecated Use convertToPower() instead.
     */
    void convertBernsteinToPower() {
        if (basis_ == Basis::BERNSTEIN) {
            *this = convertToPower();
        }
    }

    //=========================================================================
    // Evaluation
    //=========================================================================

    /**
     * @brief Evaluate the polynomial at a parameter point u ∈ [0,1]^n.
     *
     * Uses Horner's method for power basis, De Casteljau for Bernstein.
     * Parameter u is in [0,1]^n; for power basis with non-standard region,
     * coordinates are transformed automatically.
     *
     * @param parameters Parameter values for each variable (size = dimension()), in [0,1]^n.
     * @return Value of the polynomial at the given point.
     */
    Scalar evaluate(const std::vector<Scalar>& parameters) const {
        if (basis_ == Basis::POWER) {
            // For power basis, transform coordinates from [0,1] to [region_lower, region_upper]
            std::vector<Scalar> transformed = parameters;
            for (std::size_t i = 0; i < dimension_; ++i) {
                if (i < region_lower_.size() && i < region_upper_.size()) {
                    transformed[i] = region_lower_[i] + parameters[i] * (region_upper_[i] - region_lower_[i]);
                }
            }
            return detail::horner_eval_tensor(degrees_, coeffs_, transformed);
        } else {
            // Bernstein: evaluate directly on [0,1]
            return detail::de_casteljau_tensor(degrees_, coeffs_, parameters);
        }
    }

    /**
     * @brief Evaluate univariate polynomial at t.
     *
     * Convenience overload for dimension() == 1.
     */
    Scalar evaluate(const Scalar& t) const {
        std::vector<Scalar> params(1u, t);
        return evaluate(params);
    }

    //=========================================================================
    // Subdivision / Restriction
    //=========================================================================

    /**
     * @brief Restrict this polynomial to [a,b] along the given axis.
     *
     * For the specified axis, constructs a new polynomial q such that
     * q(u) = p(a + (b - a) * u), i.e., restricted to [a,b] and reparameterized
     * back to [0,1].
     *
     * - For Bernstein basis: Uses De Casteljau to modify coefficients
     * - For Power basis: Updates region bounds (O(1) operation)
     *
     * @param axis The axis to restrict (0 <= axis < dimension())
     * @param a Lower bound of interval (in [0,1])
     * @param b Upper bound of interval (in [0,1])
     * @return New polynomial restricted to [a,b] along axis
     */
    PolynomialBase restrictedToInterval(std::size_t axis, const Scalar& a, const Scalar& b) const {
        const std::size_t dim = degrees_.size();
        if (axis >= dim) {
            return *this;
        }

        // Clamp a and b to [0,1]
        Scalar a_clamped = a;
        Scalar b_clamped = b;
        if (a_clamped < Scalar(0)) a_clamped = Scalar(0);
        if (b_clamped > Scalar(1)) b_clamped = Scalar(1);

        if (!(a_clamped < b_clamped)) {
            return *this;
        }

        // For Power basis: Just update region bounds (O(1) - FAST!)
        if (basis_ == Basis::POWER) {
            PolynomialBase result = *this;
            if (axis < region_lower_.size()) {
                Scalar old_lower = region_lower_[axis];
                Scalar old_upper = region_upper_[axis];
                Scalar range = old_upper - old_lower;

                result.region_lower_[axis] = old_lower + a_clamped * range;
                result.region_upper_[axis] = old_lower + b_clamped * range;
            }
            return result;
        }

        // For Bernstein basis: Use De Casteljau subdivision
        const unsigned int deg_axis = degrees_[axis];
        const std::size_t len_axis = static_cast<std::size_t>(deg_axis + 1u);

        if (len_axis <= 1u) {
            return *this;
        }

        // Precompute strides
        std::vector<std::size_t> strides;
        detail::compute_strides(degrees_, strides);

        std::vector<Scalar> new_coeffs(coeffs_.size(), Scalar(0));

        std::vector<unsigned int> multi_index(dim, 0u);
        std::vector<Scalar> line(len_axis);
        std::vector<Scalar> left, right;

        // Iterate over all 1D slices along the given axis
        bool first = true;
        while (first || detail::increment_multi_except_axis(multi_index, degrees_, axis)) {
            first = false;

            // Gather coefficients along the current axis-line
            for (std::size_t k = 0; k < len_axis; ++k) {
                multi_index[axis] = static_cast<unsigned int>(k);
                const std::size_t idx = detail::flatten_index(multi_index, strides);
                line[k] = coeffs_[idx];
            }

            // Restrict this 1D Bernstein polynomial to [a,b]
            std::vector<Scalar> segment = line;

            if (a_clamped > Scalar(0)) {
                detail::subdivide_1d(segment, a_clamped, left, right);
                segment = right; // [a,1]
            }

            Scalar denom = Scalar(1) - a_clamped;
            Scalar t_rel = (denom > Scalar(0)) ? (b_clamped - a_clamped) / denom : Scalar(0);
            detail::subdivide_1d(segment, t_rel, left, right);

            const std::vector<Scalar>& ab_coeffs = left; // represents [a,b]

            // Scatter the restricted coefficients back into the tensor
            for (std::size_t k = 0; k < len_axis; ++k) {
                multi_index[axis] = static_cast<unsigned int>(k);
                const std::size_t idx = detail::flatten_index(multi_index, strides);
                new_coeffs[idx] = ab_coeffs[k];
            }
        }

        // Create result polynomial with updated region
        PolynomialBase result(degrees_, new_coeffs);
        result.region_lower_ = region_lower_;
        result.region_upper_ = region_upper_;

        if (axis < region_lower_.size()) {
            Scalar old_lower = region_lower_[axis];
            Scalar old_upper = region_upper_[axis];
            Scalar range = old_upper - old_lower;

            result.region_lower_[axis] = old_lower + a_clamped * range;
            result.region_upper_[axis] = old_lower + b_clamped * range;
        }

        return result;
    }

    //=========================================================================
    // Graph / Visualization
    //=========================================================================

    /**
     * @brief Build control points for the graph hypersurface.
     *
     * For each Bernstein coefficient with multi-index (i_0, ..., i_{n-1}),
     * the corresponding control point is in NORMALIZED [0,1]^n coordinates:
     * ( i_0/d_0, ..., i_{n-1}/d_{n-1}, c_{i_0,...,i_{n-1}} )
     *
     * NOTE: Coordinates are in [0,1]^n (normalized), NOT mapped to original region.
     * This avoids extra multiplications and rounding error in the solver.
     * Use getOriginalBox() to map results back to original coordinates when needed.
     *
     * Requires polynomial to be in Bernstein basis.
     * Output is a flat array of length coefficientCount() * (dimension() + 1).
     */
    void graphControlPoints(std::vector<Scalar>& control_points) const {
        // Convert to Bernstein if needed
        if (basis_ != Basis::BERNSTEIN) {
            PolynomialBase bern = convertToBernstein();
            bern.graphControlPoints(control_points);
            return;
        }

        const std::size_t dim = dimension_;

        // Handle degree-0 dimensions by duplicating for visualization
        std::vector<unsigned int> visual_degrees = degrees_;
        for (std::size_t i = 0; i < dim; ++i) {
            if (visual_degrees[i] == 0u) {
                visual_degrees[i] = 1u;
            }
        }

        std::size_t visual_count = 1u;
        for (std::size_t i = 0; i < dim; ++i) {
            visual_count *= static_cast<std::size_t>(visual_degrees[i] + 1u);
        }

        const std::size_t stride = dim + 1u;
        control_points.assign(visual_count * stride, Scalar(0));

        if (visual_count == 0u) {
            return;
        }

        std::vector<unsigned int> visual_multi(dim, 0u);
        std::vector<unsigned int> actual_multi(dim, 0u);

        for (std::size_t visual_idx = 0; visual_idx < visual_count; ++visual_idx) {
            const std::size_t base = visual_idx * stride;

            // Map visual multi-index to actual multi-index
            for (std::size_t j = 0; j < dim; ++j) {
                if (degrees_[j] == 0u) {
                    actual_multi[j] = 0u;
                } else {
                    actual_multi[j] = visual_multi[j];
                }
            }

            // Compute actual linear index
            std::size_t actual_idx = 0;
            std::size_t actual_stride = 1;
            for (std::size_t d = dim; d-- > 0;) {
                actual_idx += actual_multi[d] * actual_stride;
                actual_stride *= (degrees_[d] + 1u);
            }

            // Parameter coordinates in NORMALIZED [0,1]^n space
            for (std::size_t j = 0; j < dim; ++j) {
                Scalar coord = Scalar(0);
                if (visual_degrees[j] > 0u) {
                    coord = Scalar(visual_multi[j]) / Scalar(visual_degrees[j]);
                }
                control_points[base + j] = coord;
            }

            // Function value
            control_points[base + dim] = coeffs_[actual_idx];

            // Increment visual multi-index
            if (dim > 0u) {
                for (std::size_t d = dim; d-- > 0;) {
                    if (visual_multi[d] < visual_degrees[d]) {
                        ++visual_multi[d];
                        for (std::size_t e = d + 1; e < dim; ++e) {
                            visual_multi[e] = 0u;
                        }
                        break;
                    }
                }
            }
        }
    }

    /**
     * @brief Compute graph corners for geometric operations.
     *
     * Evaluates f at all 2^n corners of the unit hypercube and returns
     * points in R^{n+1} of the form (x_1, ..., x_n, f(x_1,...,x_n)).
     */
    void graphCorners(std::vector<std::vector<Scalar>>& corners) const {
        const std::size_t dim = dimension_;
        const std::size_t num_corners = (1u << dim);
        corners.resize(num_corners);

        for (std::size_t i = 0; i < num_corners; ++i) {
            std::vector<Scalar> point(dim + 1);
            std::vector<Scalar> params(dim);

            for (std::size_t d = 0; d < dim; ++d) {
                params[d] = ((i >> d) & 1u) ? Scalar(1) : Scalar(0);
                point[d] = params[d];
            }

            point[dim] = evaluate(params);
            corners[i] = point;
        }
    }

    //=========================================================================
    // Differentiation
    //=========================================================================

    /**
     * @brief Compute first partial derivative with respect to one variable.
     *
     * Returns a NEW polynomial representing ∂f/∂x_axis.
     * The degree is reduced by 1 along the specified axis.
     * Works for both Bernstein and Power bases.
     *
     * @param axis The variable index (0 <= axis < dimension)
     * @return The first partial derivative polynomial
     */
    PolynomialBase differentiate(std::size_t axis) const {
        if (dimension_ == 0) {
            return PolynomialBase();
        }

        std::vector<unsigned int> new_degrees;
        std::vector<Scalar> new_coeffs;

        if (basis_ == Basis::BERNSTEIN) {
            detail::differentiate_bernstein_tensor(degrees_, coeffs_, axis,
                                                    new_degrees, new_coeffs);
            PolynomialBase result(new_degrees, new_coeffs);
            // Copy region from this polynomial
            result.region_lower_ = region_lower_;
            result.region_upper_ = region_upper_;
            return result;
        } else {
            detail::differentiate_power_tensor(degrees_, coeffs_, axis,
                                                new_degrees, new_coeffs);
            PolynomialBase result = PolynomialBase::fromPower(new_degrees, new_coeffs);
            // Copy region from this polynomial
            result.region_lower_ = region_lower_;
            result.region_upper_ = region_upper_;
            return result;
        }
    }

    /**
     * @brief Compute higher-order partial derivative with respect to one variable.
     *
     * Returns a NEW polynomial representing ∂^order f / ∂x_axis^order.
     *
     * @param axis The variable index (0 <= axis < dimension)
     * @param order The derivative order (0 = identity, 1 = first derivative, etc.)
     * @return The derivative polynomial
     */
    PolynomialBase derivative(std::size_t axis, unsigned int order = 1u) const {
        if (order == 0u) {
            return *this;
        }

        // Iteratively apply first-order differentiation
        PolynomialBase result = *this;
        for (unsigned int k = 0; k < order; ++k) {
            result = result.differentiate(axis);
        }

        return result;
    }

    /**
     * @brief Compute the gradient (vector of first partial derivatives).
     *
     * Returns [∂f/∂x_0, ∂f/∂x_1, ..., ∂f/∂x_{n-1}].
     *
     * @return Vector of partial derivative polynomials, one per dimension
     */
    std::vector<PolynomialBase> gradient() const {
        std::vector<PolynomialBase> grad;
        grad.reserve(dimension_);

        for (std::size_t i = 0; i < dimension_; ++i) {
            grad.push_back(differentiate(i));
        }

        return grad;
    }

    /**
     * @brief Compute the Hessian matrix (second partial derivatives).
     *
     * Returns matrix H where H[i][j] = ∂²f/∂x_i∂x_j.
     * Note: H[i][j] = H[j][i] by symmetry of mixed partials.
     *
     * @return 2D vector representing the Hessian (dimension × dimension)
     */
    std::vector<std::vector<PolynomialBase>> hessian() const {
        std::vector<std::vector<PolynomialBase>> hess(dimension_);

        for (std::size_t i = 0; i < dimension_; ++i) {
            hess[i].resize(dimension_);
            PolynomialBase df_dxi = differentiate(i);

            for (std::size_t j = 0; j < dimension_; ++j) {
                hess[i][j] = df_dxi.differentiate(j);
            }
        }

        return hess;
    }

    //=========================================================================
    // Polynomial Arithmetic
    //=========================================================================

    /**
     * @brief Multiply two polynomials.
     *
     * Requires both polynomials to be in power basis.
     * Convert with convertToPower() first if needed.
     */
    PolynomialBase operator*(const PolynomialBase& other) const {
        if (dimension_ != other.dimension_) {
            throw std::invalid_argument("Polynomial multiplication requires same dimension");
        }

        if (dimension_ == 0) {
            return PolynomialBase();
        }

        // Require power basis
        const std::vector<Scalar>& a = powerCoefficients();
        const std::vector<Scalar>& b = other.powerCoefficients();

        std::vector<unsigned int> result_degrees;
        std::vector<Scalar> result_coeffs;
        detail::multiply_power(degrees_, a,
                               other.degrees_, b,
                               result_degrees, result_coeffs);

        return PolynomialBase::fromPower(result_degrees, result_coeffs);
    }

    /**
     * @brief Subtract two polynomials.
     *
     * Requires both polynomials to be in power basis.
     */
    PolynomialBase operator-(const PolynomialBase& other) const {
        if (dimension_ != other.dimension_) {
            throw std::invalid_argument("Polynomial subtraction requires same dimension");
        }

        if (dimension_ == 0) {
            return PolynomialBase();
        }

        const std::vector<Scalar>& a = powerCoefficients();
        const std::vector<Scalar>& b = other.powerCoefficients();

        std::vector<unsigned int> result_degrees;
        std::vector<Scalar> result_coeffs;
        detail::subtract_power(degrees_, a, other.degrees_, b,
                               result_degrees, result_coeffs);

        return PolynomialBase::fromPower(result_degrees, result_coeffs);
    }

    /**
     * @brief Add two polynomials.
     *
     * Requires both polynomials to be in power basis.
     */
    PolynomialBase operator+(const PolynomialBase& other) const {
        if (dimension_ != other.dimension_) {
            throw std::invalid_argument("Polynomial addition requires same dimension");
        }

        if (dimension_ == 0) {
            return PolynomialBase();
        }

        const std::vector<Scalar>& a = powerCoefficients();
        const std::vector<Scalar>& b = other.powerCoefficients();

        std::vector<unsigned int> result_degrees;
        std::vector<Scalar> result_coeffs;
        detail::add_power(degrees_, a, other.degrees_, b,
                          result_degrees, result_coeffs);

        return PolynomialBase::fromPower(result_degrees, result_coeffs);
    }

    /**
     * @brief Scalar multiplication.
     *
     * Requires polynomial to be in power basis.
     */
    PolynomialBase operator*(const Scalar& scalar) const {
        const std::vector<Scalar>& a = powerCoefficients();

        std::vector<Scalar> result_coeffs;
        detail::scale_coefficients(a, scalar, result_coeffs);
        return PolynomialBase::fromPower(degrees_, result_coeffs);
    }

    /**
     * @brief Unary negation.
     */
    PolynomialBase operator-() const {
        return (*this) * Scalar(-1);
    }

private:
    /// Number of variables
    std::size_t dimension_;

    /// Degree per variable
    std::vector<unsigned int> degrees_;

    /// Coefficients (single representation - power OR Bernstein)
    std::vector<Scalar> coeffs_;

    /// Which basis the coefficients are stored in
    Basis basis_;

    /// Region bounds - the polynomial evaluates on [0,1]^n but represents
    /// the function on [region_lower, region_upper]
    std::vector<Scalar> region_lower_;
    std::vector<Scalar> region_upper_;
};

//=============================================================================
// Conversion Utilities
//=============================================================================

/**
 * @brief Convert PolynomialBase between scalar types.
 *
 * @tparam ToScalar Target scalar type
 * @tparam FromScalar Source scalar type
 * @param poly Source polynomial
 * @return Polynomial with converted coefficients
 */
template<typename ToScalar, typename FromScalar>
PolynomialBase<ToScalar> convertPolynomial(const PolynomialBase<FromScalar>& poly) {
    const std::vector<FromScalar>& src_coeffs = poly.coefficients();
    std::vector<ToScalar> dst_coeffs;
    dst_coeffs.reserve(src_coeffs.size());
    for (const FromScalar& c : src_coeffs) {
        dst_coeffs.push_back(static_cast<ToScalar>(c));
    }

    if (poly.isPower()) {
        return PolynomialBase<ToScalar>::fromPower(poly.degrees(), dst_coeffs);
    } else {
        return PolynomialBase<ToScalar>::fromBernstein(poly.degrees(), dst_coeffs);
    }
}

//=============================================================================
// Derivative Cache
//=============================================================================

/**
 * @class DerivativeCacheBase
 * @brief Caches computed derivatives of a polynomial for efficient reuse.
 *
 * This class stores a polynomial and all its computed derivatives in a map,
 * allowing efficient computation of higher-order and mixed partial derivatives
 * through iterative application of first-order derivatives.
 *
 * Derivatives are indexed by a multi-index (k_0, k_1, ..., k_{n-1}) where
 * k_i is the derivative order with respect to variable i.
 *
 * Example:
 *   - {0, 0, 0} = f (original polynomial)
 *   - {1, 0, 0} = ∂f/∂x_0
 *   - {0, 1, 0} = ∂f/∂x_1
 *   - {2, 0, 0} = ∂²f/∂x_0²
 *   - {1, 1, 0} = ∂²f/∂x_0∂x_1
 *
 * @tparam Scalar The coefficient/evaluation type
 */
template<typename Scalar>
class DerivativeCacheBase {
public:
    /**
     * @brief Construct a cache for the given polynomial.
     *
     * @param p The polynomial to cache derivatives for
     */
    explicit DerivativeCacheBase(const PolynomialBase<Scalar>& p)
        : dimension_(p.dimension())
    {
        // Store the original polynomial as derivative order {0, 0, ..., 0}
        std::vector<unsigned int> zero_orders(dimension_, 0u);
        cache_[zero_orders] = p;
    }

    /**
     * @brief Get a derivative by multi-index.
     *
     * The multi-index specifies the derivative order for each variable.
     * If the derivative hasn't been computed yet, it will be computed
     * and cached automatically.
     *
     * @param orders Multi-index of derivative orders (size = dimension)
     * @return Reference to the cached derivative polynomial
     */
    const PolynomialBase<Scalar>& get(const std::vector<unsigned int>& orders) {
        if (orders.size() != dimension_) {
            throw std::invalid_argument("DerivativeCacheBase::get: orders size must match dimension");
        }
        return computeIfNeeded(orders);
    }

    /**
     * @brief Get a partial derivative with respect to a single variable.
     *
     * Convenience function for ∂^order f / ∂x_axis^order.
     *
     * @param axis The variable index
     * @param order The derivative order (default: 1)
     * @return Reference to the cached derivative polynomial
     */
    const PolynomialBase<Scalar>& getPartial(std::size_t axis, unsigned int order = 1u) {
        if (axis >= dimension_) {
            throw std::invalid_argument("DerivativeCacheBase::getPartial: axis out of range");
        }

        std::vector<unsigned int> orders(dimension_, 0u);
        orders[axis] = order;
        return get(orders);
    }

    /**
     * @brief Precompute all derivatives up to a given total order.
     *
     * This computes all derivatives where sum(orders) <= maxOrder.
     * Useful for precomputing all derivatives needed for multiplicity detection.
     *
     * @param maxOrder Maximum total derivative order
     */
    void precomputeUpToOrder(unsigned int maxOrder) {
        std::vector<unsigned int> orders(dimension_, 0u);
        precomputeRecursive(0, maxOrder, orders);
    }

    /**
     * @brief Get the dimension of the cached polynomial.
     */
    std::size_t dimension() const {
        return dimension_;
    }

private:
    std::size_t dimension_;
    std::map<std::vector<unsigned int>, PolynomialBase<Scalar>> cache_;

    void precomputeRecursive(std::size_t axis, unsigned int remaining,
                              std::vector<unsigned int>& orders) {
        if (axis == dimension_) {
            computeIfNeeded(orders);
            return;
        }

        for (unsigned int k = 0; k <= remaining; ++k) {
            orders[axis] = k;
            precomputeRecursive(axis + 1, remaining - k, orders);
        }
    }

    const PolynomialBase<Scalar>& computeIfNeeded(const std::vector<unsigned int>& orders) {
        auto it = cache_.find(orders);
        if (it != cache_.end()) {
            return it->second;
        }

        // Find which axis to differentiate
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
            throw std::logic_error("DerivativeCacheBase: zero-order derivative not found");
        }

        // Recursively get the previous derivative
        const PolynomialBase<Scalar>& prev = computeIfNeeded(prev_orders);

        // Differentiate once more
        PolynomialBase<Scalar> result = prev.differentiate(diff_axis);

        // Cache and return
        cache_[orders] = result;
        return cache_[orders];
    }
};

} // namespace polynomial_solver

#endif // POLYNOMIAL_BASE_H
