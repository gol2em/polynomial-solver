#!/usr/bin/env python3
"""
Simple test to verify modified Newton implementation.
"""

import numpy as np
from numpy.polynomial import Polynomial

def evaluate_poly(coeffs, x):
    """Evaluate polynomial with power basis coefficients at x."""
    return sum(c * x**i for i, c in enumerate(coeffs))

def evaluate_derivative(coeffs, x, order=1):
    """Evaluate derivative of polynomial at x."""
    deriv_coeffs = list(coeffs)
    for _ in range(order):
        deriv_coeffs = [(i+1) * deriv_coeffs[i+1] for i in range(len(deriv_coeffs)-1)]
        if len(deriv_coeffs) == 0:
            return 0.0
    return evaluate_poly(deriv_coeffs, x)

# Test polynomial: (x - 0.6)^6
p = Polynomial([-0.6, 1.0]) ** 6
coeffs = p.coef.tolist()

print("Testing polynomial: (x - 0.6)^6")
print(f"Coefficients: {coeffs}\n")

# Test at x = 0.59
x = 0.59
print(f"At x = {x}:")
print(f"  f(x)   = {evaluate_poly(coeffs, x):.6e}")
print(f"  f'(x)  = {evaluate_derivative(coeffs, x, 1):.6e}")
print(f"  f''(x) = {evaluate_derivative(coeffs, x, 2):.6e}")
print(f"  f'''(x) = {evaluate_derivative(coeffs, x, 3):.6e}")
print(f"  f^(4)(x) = {evaluate_derivative(coeffs, x, 4):.6e}")
print(f"  f^(5)(x) = {evaluate_derivative(coeffs, x, 5):.6e}")
print(f"  f^(6)(x) = {evaluate_derivative(coeffs, x, 6):.6e}")

# Test at exact root x = 0.6
x = 0.6
print(f"\nAt x = {x} (exact root):")
print(f"  f(x)   = {evaluate_poly(coeffs, x):.6e}")
print(f"  f'(x)  = {evaluate_derivative(coeffs, x, 1):.6e}")
print(f"  f''(x) = {evaluate_derivative(coeffs, x, 2):.6e}")
print(f"  f'''(x) = {evaluate_derivative(coeffs, x, 3):.6e}")
print(f"  f^(4)(x) = {evaluate_derivative(coeffs, x, 4):.6e}")
print(f"  f^(5)(x) = {evaluate_derivative(coeffs, x, 5):.6e}")
print(f"  f^(6)(x) = {evaluate_derivative(coeffs, x, 6):.6e}")

# Verify using numpy
print("\n=== Verification with numpy ===")
x = 0.59
print(f"At x = {x}:")
for i in range(7):
    p_deriv = p.deriv(i)
    val = p_deriv(x)
    print(f"  f^({i})(x) = {val:.6e}")

# Test modified Newton step
print("\n=== Modified Newton Step ===")
x = 0.59
f_val = evaluate_poly(coeffs, x)
df_val = evaluate_derivative(coeffs, x, 1)

print(f"x = {x}")
print(f"f(x) = {f_val:.6e}")
print(f"f'(x) = {df_val:.6e}")

for m in [1, 2, 3, 4, 5, 6]:
    step = m * f_val / df_val
    x_new = x - step
    print(f"  m={m}: step = {step:.6e}, x_new = {x_new:.6f}, error = {abs(x_new - 0.6):.6e}")

# Test a few iterations with m=6
print("\n=== Modified Newton with m=6 (multiple iterations) ===")
x = 0.59
for i in range(10):
    f_val = evaluate_poly(coeffs, x)
    df_val = evaluate_derivative(coeffs, x, 1)
    step = 6 * f_val / df_val
    x = x - step
    error = abs(x - 0.6)
    print(f"  Iter {i}: x = {x:.17f}, f(x) = {f_val:.6e}, step = {step:.6e}, error = {error:.6e}")
    if abs(step) < 1e-15:
        break

