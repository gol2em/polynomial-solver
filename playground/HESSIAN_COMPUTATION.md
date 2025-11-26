# Hessian Determinant Computation

## Overview

This document explains exactly how the Hessian determinant is computed for finding the zero set that separates convex/concave regions from saddle point regions.

## Mathematical Background

### The Hessian Matrix

For a function f(x, y), the **Hessian matrix** is the matrix of second partial derivatives:

```
H = [[f_xx,  f_xy],
     [f_xy,  f_yy]]
```

Where:
- `f_xx = ∂²f/∂x²` (second derivative with respect to x)
- `f_yy = ∂²f/∂y²` (second derivative with respect to y)
- `f_xy = ∂²f/∂x∂y` (mixed partial derivative, also called cross derivative)

### Determinant of the Hessian

The determinant is:

```
det(H) = f_xx * f_yy - f_xy²
```

### Interpretation

The sign of det(H) tells us about the local curvature:

- **det(H) > 0**: The function is locally convex or concave
  - If f_xx > 0: local minimum (convex)
  - If f_xx < 0: local maximum (concave)
  
- **det(H) < 0**: The function has saddle point behavior
  - Curves in different directions have opposite curvatures
  
- **det(H) = 0**: Boundary between these regions
  - This is the **zero set** we're finding!

## Numerical Computation

Since we don't have analytical derivatives, we use **finite differences** with step size h = 1e-6.

### Second Derivative f_xx

Using the central difference formula:

```
f_xx = (f(x+h, y) - 2*f(x, y) + f(x-h, y)) / h²
```

**Derivation**: This comes from Taylor expansion:
- f(x+h, y) ≈ f(x,y) + h*f_x + (h²/2)*f_xx + O(h³)
- f(x-h, y) ≈ f(x,y) - h*f_x + (h²/2)*f_xx + O(h³)
- Adding: f(x+h,y) + f(x-h,y) ≈ 2*f(x,y) + h²*f_xx
- Solving: f_xx ≈ (f(x+h,y) - 2*f(x,y) + f(x-h,y)) / h²

### Second Derivative f_yy

Similarly for y:

```
f_yy = (f(x, y+h) - 2*f(x, y) + f(x, y-h)) / h²
```

### Mixed Derivative f_xy

Using the 4-point formula:

```
f_xy = (f(x+h, y+h) - f(x+h, y-h) - f(x-h, y+h) + f(x-h, y-h)) / (4*h²)
```

**Derivation**: This comes from applying the difference operator twice:
- First compute ∂f/∂x at y+h and y-h
- Then compute ∂(∂f/∂x)/∂y using central difference
- Result: f_xy ≈ [f(x+h,y+h) - f(x-h,y+h) - f(x+h,y-h) + f(x-h,y-h)] / (4h²)

### Final Determinant

```
det(H) = f_xx * f_yy - f_xy²
```

## Implementation

### C++ Implementation (test_hessian_determinant.cpp)

```cpp
double second_derivative(double (*func)(double, double), 
                        double u, double v, int var1, int var2, 
                        double h = 1e-6) {
    if (var1 == 0 && var2 == 0) {
        // d²f/du² (second derivative w.r.t. u)
        return (func(u + h, v) - 2.0 * func(u, v) + func(u - h, v)) / (h * h);
    } else if (var1 == 1 && var2 == 1) {
        // d²f/dv² (second derivative w.r.t. v)
        return (func(u, v + h) - 2.0 * func(u, v) + func(u, v - h)) / (h * h);
    } else {
        // d²f/dudv (mixed derivative)
        return (func(u + h, v + h) - func(u + h, v - h) - 
                func(u - h, v + h) + func(u - h, v - h)) / (4.0 * h * h);
    }
}

double hessian_determinant(double u, double v) {
    double fuu = second_derivative(f_transformed, u, v, 0, 0);
    double fvv = second_derivative(f_transformed, u, v, 1, 1);
    double fuv = second_derivative(f_transformed, u, v, 0, 1);
    
    return fuu * fvv - fuv * fuv;
}
```

### Python Implementation (visualize_hessian_curve.py)

```python
def compute_hessian_determinant_numerical(x, y, h=1e-6):
    # Second derivative with respect to x
    f_xx = (f_original(x+h, y) - 2*f_original(x, y) + f_original(x-h, y)) / (h*h)
    
    # Second derivative with respect to y
    f_yy = (f_original(x, y+h) - 2*f_original(x, y) + f_original(x, y-h)) / (h*h)
    
    # Mixed derivative
    f_xy = (f_original(x+h, y+h) - f_original(x+h, y-h) - 
            f_original(x-h, y+h) + f_original(x-h, y-h)) / (4*h*h)
    
    # Determinant
    det_H = f_xx * f_yy - f_xy * f_xy
    
    return det_H, f_xx, f_yy, f_xy
```

## Verification

The computation has been verified in multiple ways:

1. **Mathematical verification**: The formula det(H) = f_xx * f_yy - f_xy² is printed at test points
2. **Contour comparison**: Python's contour method produces the same zero curve
3. **Residual checking**: Points on the zero set have |det(H)| ≤ threshold

## Example Computation

At point (x=0, y=0):
```
f(0, 0) = 3.141593
f_xx = 0.000000
f_yy = 0.000000
f_xy = 0.000000
det(H) = 0.000000 * 0.000000 - 0.000000² = 0.000000
```

This point is on the zero set!

## Accuracy Considerations

- **Step size h = 1e-6**: Small enough for accuracy, large enough to avoid numerical cancellation
- **Central differences**: O(h²) accuracy (better than forward/backward differences)
- **4-point mixed derivative**: More accurate than 2-point approximations
- **Residual filtering**: Ensures |det(H)| ≤ threshold for all kept points

