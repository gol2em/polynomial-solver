# Region Tracking Design

## Overview

This document explains how region tracking works for both Power and Bernstein bases to enable the complete workflow:
1. User defines polynomial on arbitrary region [a,b]^n
2. Convert to Bernstein for subdivision
3. Solver subdivides to find roots
4. Reconstruct box in original coordinates

---

## Key Principle

**Bernstein coefficients always represent polynomial on [0,1]^n**

This simplifies the implementation because:
- No need to track "coefficient region" for Bernstein
- Subdivision always works on [0,1]^n
- Only need to track mapping to original coordinates

---

## Data Structure

```cpp
template<typename Scalar>
class PolynomialBase {
private:
    std::vector<Scalar> coeffs_;
    Basis basis_;
    
    // For POWER basis only
    std::vector<Scalar> current_region_lower_;
    std::vector<Scalar> current_region_upper_;
    std::vector<Scalar> original_region_lower_;
    std::vector<Scalar> original_region_upper_;
    
    // For BERNSTEIN basis only
    std::vector<Scalar> original_box_lower_;
    std::vector<Scalar> original_box_upper_;
};
```

---

## Power Basis Region Tracking

### Interpretation

- `coeffs`: Represent polynomial on `original_region` (NEVER changes)
- `current_region`: Current region for evaluation
- `original_region`: What region coefficients represent

### Subdivision

```cpp
auto sub = power.restrictedToInterval(axis, a, b);
// coeffs: UNCHANGED (still on original_region)
// current_region: UPDATED to subdivided region
// original_region: UNCHANGED
```

**Cost**: O(1) - just update bounds

### Evaluation

```cpp
Scalar evaluate(const std::vector<Scalar>& point) {
    // point is in [0,1]^n
    // Rescale to current_region
    for (size_t i = 0; i < dimension_; ++i) {
        Scalar range = current_region_upper_[i] - current_region_lower_[i];
        rescaled_point[i] = current_region_lower_[i] + point[i] * range;
    }
    // Evaluate using original coefficients
    return evaluate_power(coeffs_, rescaled_point);
}
```

---

## Bernstein Basis Region Tracking

### Interpretation

- `coeffs`: Always represent polynomial on [0,1]^n
- `original_box`: Mapping to original coordinates for reconstruction

### Subdivision

```cpp
auto sub = bern.restrictedToInterval(axis, a, b);
// coeffs: RECOMPUTED via De Casteljau (still on [0,1]^n)
// original_box: UPDATED to reflect subdivided region in original coordinates
```

**Cost**: O(n²) - De Casteljau + update original_box

### Evaluation

```cpp
Scalar evaluate(const std::vector<Scalar>& point) {
    // point is in [0,1]^n
    // Coefficients are on [0,1]^n
    // Direct evaluation
    return evaluate_bernstein(coeffs_, point);
}
```

---

## Conversion

### Power → Bernstein

```cpp
auto bern = power.convertToBernstein();

// Steps:
// 1. Normalize power coeffs from original_region to [0,1]^n
auto normalized = normalize_power_coefficients(
    coeffs_, original_region_lower_, original_region_upper_);

// 2. Convert to Bernstein
auto bern_coeffs = power_to_bernstein(normalized);

// 3. Set original_box = current_region (for reconstruction)
bern.coeffs = bern_coeffs;  // on [0,1]^n
bern.original_box = power.current_region;
```

### Bernstein → Power

```cpp
auto power = bern.convertToPower();

// Steps:
// 1. Convert Bernstein coeffs to power (both on [0,1]^n)
auto power_coeffs = bernstein_to_power(coeffs_);

// 2. Set regions
power.coeffs = power_coeffs;  // on [0,1]^n
power.original_region = [0,1]^n;
power.current_region = bern.original_box;  // for reconstruction
```

---

## Complete Workflow Example

```cpp
// Step 1: User defines power polynomial on [2,5] × [3,7]
auto power = PolynomialBase::fromPower(degrees, coeffs, {2,3}, {5,7});
// power.coeffs: on [2,5]×[3,7]
// power.original_region: [2,5]×[3,7]
// power.current_region: [2,5]×[3,7]

// Step 2: Convert to Bernstein
auto bern = power.convertToBernstein();
// bern.coeffs: on [0,1]²
// bern.original_box: [2,5]×[3,7]

// Step 3: Subdivide (left half of axis 0)
auto sub1 = bern.restrictedToInterval(0, 0.0, 0.5);
// sub1.coeffs: on [0,1]²
// sub1.original_box: [2,3.5]×[3,7]  ✓

// Step 4: Subdivide (top half of axis 1)
auto sub2 = sub1.restrictedToInterval(1, 0.5, 1.0);
// sub2.coeffs: on [0,1]²
// sub2.original_box: [2,3.5]×[5,7]  ✓

// Step 5: Reconstruct box
auto [box_lower, box_upper] = sub2.getOriginalBox();
// box_lower = [2, 5]
// box_upper = [3.5, 7]
// ✓ Correct box in original [2,5]×[3,7] coordinates!
```

---

## Summary

| Aspect | Power Basis | Bernstein Basis |
|--------|-------------|-----------------|
| **Coefficients** | On `original_region` | Always on [0,1]^n |
| **Region tracking** | `current_region` + `original_region` | `original_box` only |
| **Subdivision** | O(1) - update `current_region` | O(n²) - recompute coeffs, update `original_box` |
| **Evaluation** | Rescale point to `current_region` | Direct (on [0,1]^n) |
| **Reconstruction** | Return `current_region` | Return `original_box` |

**Key insight**: Bernstein coefficients always on [0,1]^n simplifies implementation while `original_box` enables reconstruction!

