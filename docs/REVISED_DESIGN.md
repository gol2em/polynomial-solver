# Revised Design: Both Bases as Solver Tools

## Key Principle: No "Winner" - Both Are Tools

Both bases serve different purposes in the solver:
- **Power basis**: Efficient subdivision via region tracking (O(1))
- **Bernstein basis**: Efficient subdivision via De Casteljau, convex hull property

**Design goal**: Support both bases efficiently, with explicit conversion and region tracking.

---

## Region Tracking Philosophy

### Bernstein Basis
- **Coefficients**: Represent polynomial on **current region** `[a,b]^n`
- **After subdivision**: Coefficients recomputed (De Casteljau), region updated
- **Region tracking**: Essential for recovery and conversion

### Power Basis
- **Coefficients**: Represent polynomial on **original region** (unchanged)
- **After subdivision**: Coefficients unchanged, region updated
- **Evaluation**: Rescale point to original region
- **Region tracking**: Essential for evaluation and conversion

---

## Operation Comparison

| Operation | Power [a,b]^n | Bernstein [a,b]^n | Notes |
|-----------|---------------|-------------------|-------|
| **Subdivision** | O(1) - update region | O(n²) - De Casteljau | Power: coeffs unchanged<br/>Bernstein: coeffs recomputed |
| **Evaluation** | O(n) - Horner + rescale | O(n²) - De Casteljau | Power: rescale to original region<br/>Bernstein: direct evaluation |
| **Differentiation** | O(n) - simple formula | O(n) - simple formula | Both excellent |
| **Division** | O(d·e) - standard | O(d·e) - homogenization | Power: simpler<br/>Bernstein: more complex |
| **Convex hull** | ❌ Not available | ✅ Coefficients bound values | Bernstein only |

**Key insight**: Both bases are useful tools - choose based on solver needs!

---

## Conversion Rules

### Power → Bernstein

**Requirement**: Power coefficients must be normalized to [0,1]^n first

```cpp
// Power polynomial on [2, 5]
auto power = PolynomialBase::fromPower(degrees, coeffs, region=[2,5]);
// power.coeffs = coefficients on [2, 5] (original)
// power.region = [2, 5]

// Convert to Bernstein
auto bern = power.convertToBernstein();
// Step 1: Normalize power coeffs from [2,5] to [0,1]
// Step 2: Convert normalized coeffs to Bernstein
// bern.coeffs = Bernstein coefficients on [2, 5]
// bern.region = [2, 5] (for recovery)
```

### Bernstein → Power

**Result**: Power coefficients represent polynomial on [0,1]^n

```cpp
// Bernstein polynomial on [3, 4]
auto bern = /* ... */;
// bern.coeffs = Bernstein coefficients on [3, 4]
// bern.region = [3, 4]

// Convert to power
auto power = bern.convertToPower();
// power.coeffs = power coefficients on [0, 1]
// power.region = [3, 4] (for recovery)
```

---

## Design: Explicit Conversion with Region Tracking

### Region Tracking for Both Bases

```cpp
// User provides data
auto power = PolynomialBase::fromPower(degrees, coeffs, region=[2,5]);
auto bern = PolynomialBase::fromBernstein(degrees, coeffs);  // region defaults to [0,1]^n
```

### API Design

```cpp
template<typename Scalar, typename ErrorScalar = double>
class PolynomialBase {
private:
    std::vector<Scalar> coeffs_;
    Basis basis_;
    std::vector<unsigned int> degrees_;
    std::size_t dimension_;

    // Region tracking (for POWER basis only)
    std::vector<Scalar> current_region_lower_;
    std::vector<Scalar> current_region_upper_;
    std::vector<Scalar> original_region_lower_;
    std::vector<Scalar> original_region_upper_;

    // Original region for reconstruction (for BERNSTEIN basis only)
    std::vector<Scalar> original_box_lower_;
    std::vector<Scalar> original_box_upper_;

    // Interpretation:
    // - Bernstein:
    //   * coeffs always represent polynomial on [0,1]^n
    //   * original_box tracks mapping to original coordinates for reconstruction
    //   * After subdivision [a,b], coeffs represent polynomial on [0,1]^n
    //     and original_box is updated to reflect the subdivided region
    // - Power:
    //   * coeffs represent polynomial on original_region (unchanged after subdivision)
    //   * current_region tracks CURRENT region for evaluation rescaling
    //   * original_region tracks what region coefficients represent

    ErrorBound<ErrorScalar> error_bound_;
    
public:
    //=========================================================================
    // CONSTRUCTION
    //=========================================================================

    static PolynomialBase fromPower(
        const std::vector<unsigned int>& degrees,
        const std::vector<Scalar>& coeffs,
        const std::vector<Scalar>& region_lower = {},  // Defaults to [0,1]^n
        const std::vector<Scalar>& region_upper = {}
    ) {
        PolynomialBase result;
        result.coeffs_ = coeffs;
        result.basis_ = Basis::POWER;
        result.degrees_ = degrees;
        result.dimension_ = degrees.size();

        // Set region (default to [0,1]^n if not specified)
        if (region_lower.empty()) {
            result.current_region_lower_.assign(result.dimension_, Scalar(0));
            result.current_region_upper_.assign(result.dimension_, Scalar(1));
            result.original_region_lower_.assign(result.dimension_, Scalar(0));
            result.original_region_upper_.assign(result.dimension_, Scalar(1));
        } else {
            result.current_region_lower_ = region_lower;
            result.current_region_upper_ = region_upper;
            result.original_region_lower_ = region_lower;
            result.original_region_upper_ = region_upper;
        }

        return result;
    }

    static PolynomialBase fromBernstein(
        const std::vector<unsigned int>& degrees,
        const std::vector<Scalar>& coeffs,
        const std::vector<Scalar>& original_box_lower = {},  // Defaults to [0,1]^n
        const std::vector<Scalar>& original_box_upper = {}
    ) {
        PolynomialBase result;
        result.coeffs_ = coeffs;
        result.basis_ = Basis::BERNSTEIN;
        result.degrees_ = degrees;
        result.dimension_ = degrees.size();

        // Bernstein coefficients always represent polynomial on [0,1]^n
        // original_box tracks mapping to original coordinates for reconstruction
        if (original_box_lower.empty()) {
            result.original_box_lower_.assign(result.dimension_, Scalar(0));
            result.original_box_upper_.assign(result.dimension_, Scalar(1));
        } else {
            result.original_box_lower_ = original_box_lower;
            result.original_box_upper_ = original_box_upper;
        }

        return result;
    }

    //=========================================================================
    // SUBDIVISION (different behavior for each basis)
    //=========================================================================

    PolynomialBase restrictedToInterval(size_t axis, Scalar a, Scalar b) const {
        if (basis_ == Basis::POWER) {
            // Power: Update CURRENT region bounds (O(1))
            // Coefficients remain UNCHANGED (still on original_region)
            PolynomialBase result = *this;

            Scalar old_lower = result.current_region_lower_[axis];
            Scalar old_upper = result.current_region_upper_[axis];
            Scalar range = old_upper - old_lower;

            result.current_region_lower_[axis] = old_lower + a * range;
            result.current_region_upper_[axis] = old_lower + b * range;

            // original_region unchanged

            return result;

        } else {
            // Bernstein: Recompute coefficients via De Casteljau (O(n²))
            // Input: coeffs on [0,1]^n, subdivide axis to [a,b]
            // Output: new coeffs on [0,1]^n representing subdivided polynomial
            auto new_coeffs = subdivide_bernstein(
                degrees_, coeffs_, axis, a, b);

            PolynomialBase result;
            result.coeffs_ = new_coeffs;
            result.basis_ = Basis::BERNSTEIN;
            result.degrees_ = degrees_;
            result.dimension_ = dimension_;

            // Update original_box to reflect subdivided region
            Scalar old_lower = original_box_lower_[axis];
            Scalar old_upper = original_box_upper_[axis];
            Scalar range = old_upper - old_lower;

            result.original_box_lower_ = original_box_lower_;
            result.original_box_upper_ = original_box_upper_;
            result.original_box_lower_[axis] = old_lower + a * range;
            result.original_box_upper_[axis] = old_lower + b * range;

            result.error_bound_ = /* updated */;
            return result;
        }
    }
    
    //=========================================================================
    // DIFFERENTIATION (works for both bases)
    //=========================================================================
    
    PolynomialBase differentiate(size_t axis, unsigned int order = 1) const {
        if (basis_ == Basis::BERNSTEIN) {
            // Bernstein: d_i = n * (b_{i+1} - b_i)
            auto deriv_coeffs = differentiate_bernstein(degrees_, coeffs_, axis, order);
            
            PolynomialBase result;
            result.coeffs_ = deriv_coeffs;
            result.basis_ = Basis::BERNSTEIN;
            result.degrees_ = degrees_;  // Degree reduced
            result.dimension_ = dimension_;
            result.error_bound_ = /* updated */;
            return result;
            
        } else {
            // Power: d_i = (i+1) * a_{i+1}
            // Works on ANY region [a,b]^n
            auto deriv_coeffs = differentiate_power(degrees_, coeffs_, axis, order);
            
            PolynomialBase result;
            result.coeffs_ = deriv_coeffs;
            result.basis_ = Basis::POWER;
            result.degrees_ = degrees_;  // Degree reduced
            result.dimension_ = dimension_;
            result.region_lower_ = region_lower_;  // Preserve region
            result.region_upper_ = region_upper_;
            result.error_bound_ = /* updated */;
            return result;
        }
    }
    
    //=========================================================================
    // DIVISION (requires power basis on [0,1]^n)
    //=========================================================================
    
    std::pair<PolynomialBase, PolynomialBase> 
    divideWithRemainder(const PolynomialBase& divisor) const {
        // Check preconditions
        if (basis_ != Basis::POWER) {
            throw std::runtime_error(
                "Division requires power basis. "
                "Call convertToPower() first.");
        }
        
        if (!hasStandardRegion()) {
            throw std::runtime_error(
                "Division requires standard region [0,1]^n. "
                "Call normalizeToStandardRegion() first.");
        }
        
        if (divisor.basis_ != Basis::POWER || !divisor.hasStandardRegion()) {
            throw std::runtime_error(
                "Divisor must be in power basis on [0,1]^n.");
        }
        
        // Divide in power basis
        auto [quot_coeffs, rem_coeffs] = divide_power(
            degrees_, coeffs_,
            divisor.degrees_, divisor.coeffs_);
        
        return {
            PolynomialBase::fromPower(quot_degrees, quot_coeffs),
            PolynomialBase::fromPower(rem_degrees, rem_coeffs)
        };
    }
    
    //=========================================================================
    // EXPLICIT CONVERSION
    //=========================================================================

    PolynomialBase convertToBernstein() const {
        if (basis_ == Basis::BERNSTEIN) return *this;

        // Power coefficients represent polynomial on original_region
        // Need to normalize to [0,1]^n first

        // Step 1: Normalize power coefficients from original_region to [0,1]^n
        auto normalized_coeffs = normalize_power_coefficients(
            degrees_, coeffs_, original_region_lower_, original_region_upper_);

        // Step 2: Convert normalized coeffs to Bernstein
        auto bern_coeffs = power_to_bernstein(degrees_, normalized_coeffs);

        // Step 3: Create Bernstein polynomial
        // - Bernstein coeffs represent polynomial on [0,1]^n
        // - original_box = current_region (for reconstruction)
        return PolynomialBase::fromBernstein(
            degrees_, bern_coeffs, current_region_lower_, current_region_upper_);
    }

    PolynomialBase convertToPower() const {
        if (basis_ == Basis::POWER) return *this;

        // Bernstein coeffs represent polynomial on [0,1]^n
        // Convert to power (also on [0,1]^n)
        auto power_coeffs = bernstein_to_power(degrees_, coeffs_);

        // Create power polynomial
        // - Power coeffs represent polynomial on [0,1]^n
        // - current_region = original_box (for reconstruction)
        // - original_region = [0,1]^n
        PolynomialBase result = PolynomialBase::fromPower(degrees_, power_coeffs);
        result.current_region_lower_ = original_box_lower_;
        result.current_region_upper_ = original_box_upper_;
        return result;
    }

    //=========================================================================
    // RECONSTRUCTION
    //=========================================================================

    std::pair<std::vector<Scalar>, std::vector<Scalar>> getOriginalBox() const {
        if (basis_ == Basis::BERNSTEIN) {
            return {original_box_lower_, original_box_upper_};
        } else {
            return {current_region_lower_, current_region_upper_};
        }
    }

    PolynomialBase normalizeToStandardRegion() const {
        if (basis_ == Basis::BERNSTEIN) {
            // Bernstein coefficients already represent polynomial on current region
            // Just update region to [0,1]^n
            return PolynomialBase::fromBernstein(degrees_, coeffs_);
        } else {
            // Power: transform coefficients to [0,1]^n
            auto normalized_coeffs = normalize_power_coefficients(
                degrees_, coeffs_, region_lower_, region_upper_);

            return PolynomialBase::fromPower(degrees_, normalized_coeffs);
        }
    }
    
    //=========================================================================
    // EVALUATION (works for both bases, handles regions)
    //=========================================================================

    Scalar evaluate(const std::vector<Scalar>& point) const {
        if (basis_ == Basis::POWER) {
            // Power: coefficients are on ORIGINAL region
            // Need to rescale point from CURRENT region to ORIGINAL region

            // For now, assume point is in [0,1]^n and rescale to current region
            std::vector<Scalar> rescaled_point(dimension_);
            for (size_t i = 0; i < dimension_; ++i) {
                Scalar range = region_upper_[i] - region_lower_[i];
                rescaled_point[i] = region_lower_[i] + point[i] * range;
            }

            return evaluate_power(degrees_, coeffs_, rescaled_point);

        } else {
            // Bernstein: coefficients are on CURRENT region
            // Point is in [0,1]^n relative to current region
            return evaluate_bernstein(degrees_, coeffs_, point);
        }
    }
    
    //=========================================================================
    // CONVEX HULL (only for Bernstein)
    //=========================================================================
    
    std::pair<Scalar, Scalar> getConvexHullBounds() const {
        if (basis_ != Basis::BERNSTEIN) {
            throw std::runtime_error(
                "Convex hull property only available in Bernstein basis. "
                "Call convertToBernstein() first.");
        }
        
        auto [min_coeff, max_coeff] = std::minmax_element(
            coeffs_.begin(), coeffs_.end());
        
        return {*min_coeff, *max_coeff};
    }
};
```

---

## Complete Workflow Example

### User-Defined Polynomial → Bernstein Subdivision → Reconstruction

```cpp
// ============================================================================
// Step 1: User defines power polynomial on [2,5] × [3,7]
// ============================================================================
auto power = PolynomialBase<mpreal>::fromPower(
    degrees, coeffs,
    region_lower={2, 3},
    region_upper={5, 7});

// State:
// power.coeffs = coefficients representing P(x,y) on [2,5]×[3,7]
// power.original_region = [2,5]×[3,7]
// power.current_region = [2,5]×[3,7]

// ============================================================================
// Step 2: Convert to Bernstein
// ============================================================================
auto bern = power.convertToBernstein();

// What happens:
// 1. Normalize power coeffs from [2,5]×[3,7] to [0,1]²
// 2. Convert normalized coeffs to Bernstein
// 3. Create Bernstein polynomial

// State:
// bern.coeffs = Bernstein coefficients on [0,1]²
// bern.original_box = [2,5]×[3,7]  (for reconstruction)

// ============================================================================
// Step 3: Solver subdivides (e.g., left half of first axis)
// ============================================================================
auto sub = bern.restrictedToInterval(0, 0.0, 0.5);

// What happens:
// 1. De Casteljau subdivision on axis 0, interval [0, 0.5]
// 2. New coeffs represent polynomial on [0,1]²
// 3. Update original_box to reflect subdivided region

// State:
// sub.coeffs = Bernstein coefficients on [0,1]²
// sub.original_box = [2, 3.5]×[3,7]  ✓ (updated!)

// ============================================================================
// Step 4: Further subdivision (top half of second axis)
// ============================================================================
auto sub2 = sub.restrictedToInterval(1, 0.5, 1.0);

// State:
// sub2.coeffs = Bernstein coefficients on [0,1]²
// sub2.original_box = [2, 3.5]×[5, 7]  ✓ (updated!)

// ============================================================================
// Step 5: Reconstruct box in original region
// ============================================================================
auto [box_lower, box_upper] = sub2.getOriginalBox();

// Result:
// box_lower = [2, 5]
// box_upper = [3.5, 7]
// ✓ Correct box in original [2,5]×[3,7] coordinates!
```

---

## Usage Examples

### Example 1: Subdivision-Based Solver (Power Basis)

```cpp
// Create polynomial in power basis on [2,5]
auto poly = PolynomialBase<mpreal>::fromPower(
    degrees, coeffs, region_lower={2}, region_upper={5});

// Subdivide repeatedly (O(1) each time!)
auto left = poly.restrictedToInterval(0, 0.0, 0.5);   // [2, 3.5]
auto right = poly.restrictedToInterval(0, 0.5, 1.0);  // [3.5, 5]

// Differentiate (stays in power, preserves region)
auto dleft = left.differentiate(0);
auto dright = right.differentiate(0);

// Evaluate (automatic rescaling)
mpreal val_left = dleft.evaluate({0.25});   // Rescales to [2, 3.5]
mpreal val_right = dright.evaluate({0.75}); // Rescales to [3.5, 5]

// Get current box
auto [box_lower, box_upper] = left.getOriginalBox();
// box_lower = [2], box_upper = [3.5]
```

### Example 2: Bernstein Subdivision with Reconstruction

```cpp
// User polynomial on [2,5] × [3,7]
auto power = PolynomialBase<mpreal>::fromPower(
    degrees, coeffs, region_lower={2,3}, region_upper={5,7});

// Convert to Bernstein for subdivision
auto bern = power.convertToBernstein();
// bern.coeffs on [0,1]²
// bern.original_box = [2,5]×[3,7]

// Subdivide to find roots
auto sub1 = bern.restrictedToInterval(0, 0.0, 0.5);  // Left half
auto sub2 = sub1.restrictedToInterval(1, 0.5, 1.0);  // Top half

// Reconstruct box
auto [box_lower, box_upper] = sub2.getOriginalBox();
// box_lower = [2, 5], box_upper = [3.5, 7]
// ✓ Correct box in original coordinates!
```

### Example 3: Convex Hull (Bernstein Only)

```cpp
// Polynomial in Bernstein basis
auto bern = /* ... */;

// Get bounds from convex hull
auto [min_val, max_val] = bern.getConvexHullBounds();

// Use bounds to prune search space
if (min_val > 0 || max_val < 0) {
    // No roots in this region
    return;
}
```

---

## Summary

### Region Tracking Design

**Power Basis:**
- `coeffs`: Represent polynomial on `original_region` (unchanged after subdivision)
- `current_region`: Tracks current region for evaluation rescaling
- `original_region`: Tracks what region coefficients represent
- Subdivision: O(1) - update `current_region` only

**Bernstein Basis:**
- `coeffs`: Always represent polynomial on [0,1]^n
- `original_box`: Tracks mapping to original coordinates for reconstruction
- Subdivision: O(n²) - recompute coefficients, update `original_box`

### Conversion Rules

**Power → Bernstein:**
1. Normalize power coeffs from `original_region` to [0,1]^n
2. Convert to Bernstein
3. Set `original_box = current_region` (for reconstruction)

**Bernstein → Power:**
1. Convert Bernstein coeffs to power (on [0,1]^n)
2. Set `current_region = original_box` (for reconstruction)
3. Set `original_region = [0,1]^n`

### Key Features

✅ **Both bases are tools** (no "winner")
✅ **Explicit conversion only** (no automatic)
✅ **Region tracking** for reconstruction
✅ **Bernstein coeffs always on [0,1]^n** (simplifies implementation)
✅ **Power supports arbitrary regions** (flexible)

**Key insight**: Proper region tracking enables seamless workflow from user input → subdivision → reconstruction!

