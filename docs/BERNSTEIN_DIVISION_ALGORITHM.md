# Bernstein Division Algorithm

## Overview

This document describes the **division algorithm for Bernstein polynomials** based on the paper by Busé et al.

**Key insight**: Division in Bernstein basis uses homogenization and operations on powers of `t` and `(1-t)`.

---

## Mathematical Background

### Homogenization

For a Bernstein polynomial `f(t)` of degree `d`, the **homogenization** is:

```
f^h(t,s) = s^d · f(t/s) ∈ ℝ[t,s]
```

**Properties**:
- `f^h(1,0) ≠ 0` if `f(1) ≠ 0` (polynomial doesn't vanish at infinity)
- Dehomogenization: `f(t) = f^h(t,1)`

### Division Formula (Proposition 4.1)

Given two Bernstein polynomials:
```
f(t) = Σ(i=0 to d) a_i·B^d_i(t)    with f(1) ≠ 0
g(t) = Σ(i=0 to e) b_i·B^e_i(t)    with g(1) ≠ 0
```

There exist unique polynomials `q(t)` and `r(t)` such that:

```
g(t) = q(t)·f(t) + (1-t)^(e-d+1)·r(t)
```

where:
- `q(t) = Σ(i=0 to e-d) q_i·B^(e-d)_i(t)` (quotient)
- `r(t) = Σ(i=0 to d-1) r_i·B^(d-1)_i(t)` (remainder)
- `deg(r) < deg(f)` and `r(1) ≠ 0` (possibly)

**Note**: The remainder is multiplied by `(1-t)^(e-d+1)`, not just `r(t)` alone!

---

## Algorithm

### Input
- `f(t)` with degree `d` and `f(1) ≠ 0`
- `g(t)` with degree `e ≥ d` and `g(1) ≠ 0`

### Output
- Quotient `q(t)` with degree `e-d`
- Remainder `r(t)` with degree `< d`

### Steps

1. **Initialize**:
   ```
   r(t) ← g(t)     // degree e
   q(t) ← 0        // degree e-d
   ```

2. **While** `deg(r) ≥ deg(f)`:

   a. **Extract leading coefficient**:
      ```
      Let r(t) = b·B^deg(r)_deg(r)(t) + ... (leading coefficient is b)
      Let a = f(1) (last Bernstein coefficient of f)
      
      Set coefficient of B^(e-d)_(deg(r)-d)(t) in q(t) to:
          b / (binomial(e-d, deg(r)-d) · a)
      ```

   b. **Update remainder**:
      ```
      r(t) ← r(t) - (b/a)·t^(deg(r)-d)·f(t)
      ```
      
      This uses the multiplication formula (2.1) to compute `t^(deg(r)-d)·f(t)`.

   c. **Remove common powers of (1-t)**:
      ```
      Factor out (1-t) from r(t) to reduce degree
      ```
      
      This uses the division formula (2.4) to compute `r(t)/(1-t)`.

3. **Return** `q(t)` and `r(t)`

---

## Key Formulas

### Multiplication by `t^d`

To compute `t^d · P(t)` where `P(t) = Σ c^n_k·B^n_k(t)`:

```
(t^d · P)(t) = Σ c^(n+d)_k·B^(n+d)_k(t)
```

where the new coefficients are:

```
c^(n+d)_k = ∏(i=0 to d-1) (k-i)/(n+d-i) · c^n_(k-d)    (equation 2.1)
```

**Implementation**:
```cpp
std::vector<Scalar> multiply_by_t_power(
    const std::vector<Scalar>& coeffs, 
    unsigned int n,  // original degree
    unsigned int d   // power of t
) {
    std::vector<Scalar> result(coeffs.size() + d, Scalar(0));
    
    for (size_t k = d; k < result.size(); ++k) {
        Scalar factor = Scalar(1);
        for (unsigned int i = 0; i < d; ++i) {
            factor *= Scalar(k - i) / Scalar(n + d - i);
        }
        result[k] = factor * coeffs[k - d];
    }
    
    return result;
}
```

### Multiplication by `(1-t)^d`

To compute `(1-t)^d · P(t)` where `P(t) = Σ c^n_k·B^n_k(t)`:

```
((1-t)^d · P)(t) = Σ c^(n+d)_k·B^(n+d)_k(t)
```

where:

```
c^(n+d)_k = ∏(i=0 to d-1) (n+d-k-i)/(n+d-i) · c^n_k    (equation 2.2)
```

### Division by `t^j`

To compute `P(t)/t^j` where `P(t) = Σ(k=j to n) c^n_k·B^n_k(t)`:

```
(P/t^j)(t) = Σ c^(n-j)_k·B^(n-j)_k(t)
```

where:

```
c^(n-j)_k = ∏(i=0 to j-1) (n-i)/(k+j-i) · c^n_(k+j)    (equation 2.3)
```

### Division by `(1-t)^j`

To compute `P(t)/(1-t)^j` where `P(t) = Σ(k=0 to n-j) c^n_k·B^n_k(t)`:

```
(P/(1-t)^j)(t) = Σ c^(n-j)_k·B^(n-j)_k(t)
```

where:

```
c^(n-j)_k = ∏(i=0 to j-1) (n-i)/(n-k-i) · c^n_k    (equation 2.4)
```

---

## Example

**Divide** `g(t) = 3B^3_0 + 2B^3_1 + 4B^3_2 + 5B^3_3` by `f(t) = 1B^2_0 + 2B^2_1 + 3B^2_2`:

1. **Initialize**:
   ```
   r = [3, 2, 4, 5]  (degree 3)
   q = [0]           (degree 3-2 = 1)
   ```

2. **Iteration 1** (`deg(r) = 3 ≥ deg(f) = 2`):
   
   a. Leading coefficient of `r`: `b = 5`
      Leading coefficient of `f`: `a = f(1) = 3`
      
      Set `q_1 = 5 / (binomial(1,1) · 3) = 5/3`
   
   b. Update: `r ← r - (5/3)·t^1·f`
   
   c. Remove common `(1-t)` from `r`

3. **Iteration 2** (if `deg(r) ≥ 2`):
   
   Continue until `deg(r) < 2`

4. **Result**: `q(t)` and `r(t)` such that `g(t) = q(t)·f(t) + (1-t)^2·r(t)`

---

## Implementation Notes

### Numerical Stability

- Use **high precision** for binomial coefficients
- Accumulate products incrementally to avoid overflow
- Check for division by zero

### Degree Reduction

The "remove common powers of (1-t)" step is crucial:
- After subtraction, `r(t)` may have leading zeros
- Factor out `(1-t)` to reduce degree
- This ensures termination

### Complexity

- **Time**: O(d·e) where d = deg(f), e = deg(g)
- **Space**: O(e) for storing coefficients

---

## Summary

✅ **Division works in Bernstein basis!**
✅ **Uses homogenization** (not direct coefficient manipulation)
✅ **Key operations**: Multiply/divide by powers of `t` and `(1-t)`
✅ **Numerically stable** (uses binomial coefficient ratios)
✅ **Complexity**: Similar to power basis division

**Recommendation**: Implement this algorithm for Bernstein division instead of converting to power basis!

