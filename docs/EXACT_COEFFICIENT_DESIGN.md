# Exact Coefficient Pattern - Design Proposal (v3)

## Summary

A simplified design focused on **error tracking** with the ability to **lift precision locally** when accumulated error grows too large. The solver manages coefficient lineage externally.

---

## Core Principles

1. **Track accumulated error** - Every polynomial knows its error bound
2. **No lineage storage** - Solver handles origin tracking externally  
3. **Lift precision on demand** - When error grows too large, promote to higher precision
4. **Error stays near machine epsilon** - Automatic precision lifting keeps error bounded

---

## Error Tracking Model

```cpp
/**
 * Error bound estimate for rounding error accumulation.
 * 
 * For exact arithmetic (GMP rationals): error = 0 always
 * For floating-point: error grows with operations
 * 
 * NOTE: ErrorScalar is typically double, even when Scalar is higher precision.
 * This is fine as long as the bound is valid (conservative upper bound).
 */
template<typename ErrorScalar = double>
struct ErrorBound {
    ErrorScalar absolute_error;   // Absolute error bound on coefficients
    
    bool isExact() const { return absolute_error == ErrorScalar(0); }
    
    static ErrorBound zero() { return {ErrorScalar(0)}; }
};

// Type trait to detect exact arithmetic types
template<typename T>
struct is_exact_arithmetic : std::false_type {};

template<>
struct is_exact_arithmetic<mpq_class> : std::true_type {};  // GMP rational

// Error estimation functions
// ErrorScalar is typically double for efficiency
template<typename Scalar, typename ErrorScalar = double>
ErrorScalar errorFromConversion(
    const std::vector<unsigned int>& degrees,
    const std::vector<Scalar>& coeffs)
{
    if constexpr (is_exact_arithmetic<Scalar>::value) {
        return ErrorScalar(0);  // Exact: no error
    } else {
        // Floating-point: error = eps * degree^2 * |max_coeff|
        double max_coeff = 0;
        for (const auto& c : coeffs) {
            double abs_c = std::abs(static_cast<double>(c));
            if (abs_c > max_coeff) max_coeff = abs_c;
        }
        unsigned int total_degree = 1;
        for (auto d : degrees) total_degree *= d;
        double eps = std::numeric_limits<double>::epsilon();  // Always use double eps as baseline
        return static_cast<ErrorScalar>(eps * double(total_degree * total_degree) * max_coeff);
    }
}

template<typename Scalar, typename ErrorScalar = double>
ErrorScalar errorFromSubdivision(
    const std::vector<unsigned int>& degrees,
    const std::vector<Scalar>& coeffs)
{
    if constexpr (is_exact_arithmetic<Scalar>::value) {
        return ErrorScalar(0);
    } else {
        double max_coeff = 0;
        for (const auto& c : coeffs) {
            double abs_c = std::abs(static_cast<double>(c));
            if (abs_c > max_coeff) max_coeff = abs_c;
        }
        unsigned int max_degree = *std::max_element(degrees.begin(), degrees.end());
        double eps = std::numeric_limits<double>::epsilon();
        return static_cast<ErrorScalar>(eps * double(max_degree) * max_coeff);
    }
}
```

---

## Simplified Data Model

```cpp
enum class Basis {
    POWER,
    BERNSTEIN
};

template<typename Scalar, typename ErrorScalar = double>
class PolynomialBase {
private:
    std::vector<Scalar> coeffs_;
    Basis basis_;
    std::vector<unsigned int> degrees_;
    std::size_t dimension_;
    ErrorBound<ErrorScalar> error_bound_;  // Error stored as double for efficiency
    
public:
    // =========== CONSTRUCTION ===========
    
    static PolynomialBase fromPower(const std::vector<unsigned int>& degrees,
                                    const std::vector<Scalar>& coeffs) {
        PolynomialBase p;
        p.coeffs_ = coeffs;
        p.basis_ = Basis::POWER;
        p.degrees_ = degrees;
        p.dimension_ = degrees.size();
        p.error_bound_ = ErrorBound<ErrorScalar>::zero();  // User input is exact
        return p;
    }
    
    static PolynomialBase fromBernstein(const std::vector<unsigned int>& degrees,
                                        const std::vector<Scalar>& coeffs) {
        PolynomialBase p;
        p.coeffs_ = coeffs;
        p.basis_ = Basis::BERNSTEIN;
        p.degrees_ = degrees;
        p.dimension_ = degrees.size();
        p.error_bound_ = ErrorBound<ErrorScalar>::zero();
        return p;
    }
    
    // =========== ACCESSORS ===========
    
    Basis basis() const { return basis_; }
    const std::vector<unsigned int>& degrees() const { return degrees_; }
    std::size_t dimension() const { return dimension_; }
    const std::vector<Scalar>& coefficients() const { return coeffs_; }
    ErrorScalar errorBound() const { return error_bound_.absolute_error; }
    bool isExact() const { return error_bound_.isExact(); }
    
    // Convenience with basis check
    const std::vector<Scalar>& powerCoefficients() const {
        if (basis_ != Basis::POWER) {
            std::cerr << "ERROR: powerCoefficients() but basis is BERNSTEIN\n";
            std::abort();
        }
        return coeffs_;
    }
    
    const std::vector<Scalar>& bernsteinCoefficients() const {
        if (basis_ != Basis::BERNSTEIN) {
            std::cerr << "ERROR: bernsteinCoefficients() but basis is POWER\n";
            std::abort();
        }
        return coeffs_;
    }
    
    // =========== EXPLICIT CONVERSION ===========
    
    PolynomialBase convertToBernstein() const {
        if (basis_ == Basis::BERNSTEIN) {
            return *this;
        }
        
        PolynomialBase result;
        result.coeffs_ = power_to_bernstein(degrees_, coeffs_);
        result.basis_ = Basis::BERNSTEIN;
        result.degrees_ = degrees_;
        result.dimension_ = dimension_;
        // Error accumulates additively
        ErrorScalar conversion_err = errorFromConversion<Scalar, ErrorScalar>(degrees_, result.coeffs_);
        result.error_bound_.absolute_error = error_bound_.absolute_error + conversion_err;
        return result;
    }
    
    PolynomialBase convertToPower() const {
        if (basis_ == Basis::POWER) {
            return *this;
        }
        
        PolynomialBase result;
        result.coeffs_ = bernstein_to_power(degrees_, coeffs_);
        result.basis_ = Basis::POWER;
        result.degrees_ = degrees_;
        result.dimension_ = dimension_;
        ErrorScalar conversion_err = errorFromConversion<Scalar, ErrorScalar>(degrees_, result.coeffs_);
        result.error_bound_.absolute_error = error_bound_.absolute_error + conversion_err;
        return result;
    }
    
    // =========== SUBDIVISION ===========
    
    PolynomialBase restrictedToInterval(std::size_t axis, Scalar a, Scalar b) const {
        if (basis_ != Basis::BERNSTEIN) {
            std::cerr << "ERROR: restrictedToInterval requires Bernstein basis\n";
            std::abort();
        }
        
        std::vector<Scalar> new_coeffs = de_casteljau_restrict(degrees_, coeffs_, axis, a, b);
        
        PolynomialBase result;
        result.coeffs_ = new_coeffs;
        result.basis_ = Basis::BERNSTEIN;
        result.degrees_ = degrees_;
        result.dimension_ = dimension_;
        ErrorScalar sub_err = errorFromSubdivision<Scalar, ErrorScalar>(degrees_, new_coeffs);
        result.error_bound_.absolute_error = error_bound_.absolute_error + sub_err;
        return result;
    }
    
    // =========== EVALUATION ===========
    
    Scalar evaluate(const std::vector<Scalar>& params) const {
        if (basis_ == Basis::POWER) {
            return horner_eval(degrees_, coeffs_, params);
        } else {
            return de_casteljau_eval(degrees_, coeffs_, params);
        }
    }
    
    // =========== PRECISION LIFTING ===========
    
    /**
     * Lift to higher precision scalar type.
     * Error bound is converted but NOT reset (lifting doesn't reduce error).
     */
    template<typename TargetScalar>
    PolynomialBase<TargetScalar, ErrorScalar> liftPrecision() const {
        std::vector<TargetScalar> lifted_coeffs;
        lifted_coeffs.reserve(coeffs_.size());
        for (const Scalar& c : coeffs_) {
            lifted_coeffs.push_back(static_cast<TargetScalar>(c));
        }
        
        PolynomialBase<TargetScalar, ErrorScalar> result;
        result.coeffs_ = lifted_coeffs;
        result.basis_ = basis_;
        result.degrees_ = degrees_;
        result.dimension_ = dimension_;
        // Keep same error bound (lifting doesn't fix existing error)
        result.error_bound_ = error_bound_;
        return result;
    }
};
```

---

## Usage Patterns

### Pattern 1: Basic Operations with Error Tracking

```cpp
// Create polynomial with exact coefficients
auto poly = PolynomialBase<double>::fromPower(degrees, coeffs);
std::cout << "Initial error: " << poly.errorBound() << "\n";  // 0

// Convert to Bernstein
auto bern = poly.convertToBernstein();
std::cout << "After conversion: " << bern.errorBound() << "\n";  // small

// Subdivide multiple times
auto sub = bern;
for (int i = 0; i < 20; ++i) {
    sub = sub.restrictedToInterval(0, 0.0, 0.5);
}
std::cout << "After 20 subdivisions: " << sub.errorBound() << "\n";
```

### Pattern 2: Solver Decides When to Lift Precision

```cpp
// Solver keeps original polynomial for reconstruction
void solve(const PolynomialBase<double>& original) {
    auto working = original.convertToBernstein();
    std::vector<double> box_lower = {0.0, 0.0};
    std::vector<double> box_upper = {1.0, 1.0};
    
    // ... subdivision loop ...
    
    // Solver decides when to lift based on error bound
    const double ERROR_TOLERANCE = 1e-10;
    if (working.errorBound() > ERROR_TOLERANCE) {
        // Reconstruct from original in HP using box coordinates
        auto hp_original = original.liftPrecision<mpfr_float>();
        auto hp_bern = hp_original.convertToBernstein();
        auto hp_poly = hp_bern.restrictedToBox(box_lower, box_upper);
        // Continue with hp_poly...
    }
}
```

### Pattern 3: Exact Rational Arithmetic (No Error)

```cpp
// With GMP rationals - error stays zero
auto poly_r = PolynomialBase<mpq_class>::fromPower(degrees, rational_coeffs);
auto bern_r = poly_r.convertToBernstein();
auto sub_r = bern_r.restrictedToInterval(0, mpq_class(3,10), mpq_class(7,10));

// Still exact!
assert(sub_r.isExact());
assert(sub_r.errorBound() == 0);
```

---

## Benefits

1. **Simple data model** - Just coefficients + basis + error bound
2. **Error awareness** - Always know how much error has accumulated
3. **Exact arithmetic support** - Error stays zero for rational types
4. **Precision lifting** - Can promote to higher precision when needed
5. **Solver controls lineage** - External management, no shared_ptr overhead

---

## Implementation Notes

1. **ErrorBound is cheap** - Just 1 scalar value (absolute_error)
2. **Type trait for exact types** - Use `is_exact_arithmetic<T>` to specialize  
3. **No caching** - Single representation, explicit conversion
4. **Error preserved on lift** - `liftPrecision()` keeps existing error bound (doesn't reduce error)
5. **Solver decides when to lift** - Polynomial just reports error, doesn't decide tolerance

---

## Migration from Current Design

1. Remove `primary_rep_`, `bernstein_valid_`, `power_valid_`
2. Add `error_bound_` member
3. Remove implicit conversion logic
4. Add `liftPrecision<T>()` method
5. Update operations to compute new error bounds
6. Solver stores original polynomials + box coordinates for HP reconstruction
