#!/usr/bin/env python3
"""
Test CORRECT Ojika/Bertini deflation method using AUGMENTED SYSTEM.

Ojika's deflation for 1D polynomial f(x) (CORRECT VERSION):

Stage 0 (original):
  System: F(x) = [f(x)] = 0
  Variables: x
  Jacobian: J = [f'(x)]
  If |f'(x)| < tol → SINGULAR, need deflation

Stage 1 (first deflation):
  System: F(x, v₁) = [f(x)    ] = 0
                     [f'(x)·v₁]
                     [v₁² - 1 ]
  Variables: (x, v₁)
  Jacobian: J = [f'(x)      0     ]
                [f''(x)·v₁  f'(x) ]
                [0          2v₁   ]
  Check rank using SVD. If rank < 3 → need another deflation

Stage 2 (second deflation):
  System: F(x, v₁, v₂) = [f(x)                    ] = 0
                         [f'(x)·v₁                ]
                         [f''(x)·v₁ + f'(x)·v₂   ]  ← USES BOTH v₁ AND v₂!
                         [v₁² + v₂² - 1          ]
  Variables: (x, v₁, v₂)

Stage 3 (third deflation):
  System: F(x, v₁, v₂, v₃) = [f(x)                              ] = 0
                             [f'(x)·v₁                          ]
                             [f''(x)·v₁ + f'(x)·v₂              ]
                             [f'''(x)·v₁ + f''(x)·v₂ + f'(x)·v₃]  ← USES ALL THREE!
                             [v₁² + v₂² + v₃² - 1               ]

General pattern for stage k:
  F₀ = f(x)
  F₁ = f'(x)·v₁
  F₂ = f''(x)·v₁ + f'(x)·v₂
  F₃ = f'''(x)·v₁ + f''(x)·v₂ + f'(x)·v₃
  ...
  Fₖ = f^(k)(x)·v₁ + f^(k-1)(x)·v₂ + ... + f'(x)·vₖ
  Fₖ₊₁ = v₁² + v₂² + ... + vₖ² - 1

Key: Each equation Fᵢ uses ALL variables v₁, ..., vᵢ (dual space approach)!
"""

import numpy as np

def evaluate_poly(coeffs, x):
    """Evaluate polynomial with power basis coefficients at x."""
    return sum(c * x**i for i, c in enumerate(coeffs))

def evaluate_derivative(coeffs, x, order=1):
    """Evaluate derivative of polynomial at x."""
    deriv_coeffs = coeffs.copy()
    for _ in range(order):
        deriv_coeffs = [(i+1) * deriv_coeffs[i+1] for i in range(len(deriv_coeffs)-1)]
        if len(deriv_coeffs) == 0:
            return 0.0
    return evaluate_poly(deriv_coeffs, x)

def augmented_system_residual(f_coeffs, z, num_deflations):
    """
    Evaluate augmented deflation system (CORRECT Ojika/Bertini version).

    z = [x, v₁, v₂, ..., vₖ] where k = num_deflations

    System:
      F₀ = f(x)
      F₁ = f'(x) · v₁
      F₂ = f''(x) · v₁ + f'(x) · v₂
      F₃ = f'''(x) · v₁ + f''(x) · v₂ + f'(x) · v₃
      ...
      Fᵢ = f^(i)(x) · v₁ + f^(i-1)(x) · v₂ + ... + f'(x) · vᵢ
      ...
      Fₖ = f^(k)(x) · v₁ + f^(k-1)(x) · v₂ + ... + f'(x) · vₖ
      Fₖ₊₁ = v₁² + v₂² + ... + vₖ² - 1
    """
    x = z[0]
    v = z[1:] if num_deflations > 0 else []

    residual = []

    # F₀ = f(x)
    residual.append(evaluate_poly(f_coeffs, x))

    # Fᵢ = sum_{j=1}^{i} f^(i-j+1)(x) · vⱼ for i = 1, ..., k
    for i in range(1, num_deflations + 1):
        Fi = 0.0
        for j in range(1, i + 1):
            # f^(i-j+1)(x) · vⱼ
            deriv_order = i - j + 1
            df_val = evaluate_derivative(f_coeffs, x, order=deriv_order)
            Fi += df_val * v[j-1]  # v is 0-indexed
        residual.append(Fi)

    # Normalization constraint: ||v||² = 1
    if num_deflations > 0:
        v_norm_sq = sum(vi**2 for vi in v)
        residual.append(v_norm_sq - 1.0)

    return np.array(residual)

def augmented_system_jacobian(f_coeffs, z, num_deflations):
    """
    Compute Jacobian of augmented system (CORRECT Ojika/Bertini version).

    Variables: z = [x, v₁, v₂, ..., vₖ]

    Jacobian J[i,j] = ∂Fᵢ/∂zⱼ

    For Fᵢ = sum_{j=1}^{i} f^(i-j+1)(x) · vⱼ:
      ∂Fᵢ/∂x = sum_{j=1}^{i} f^(i-j+2)(x) · vⱼ
      ∂Fᵢ/∂vⱼ = f^(i-j+1)(x) for j ≤ i
      ∂Fᵢ/∂vⱼ = 0 for j > i
    """
    x = z[0]
    v = z[1:] if num_deflations > 0 else []

    n = 1 + num_deflations + (1 if num_deflations > 0 else 0)  # Total equations
    m = 1 + num_deflations  # Total variables

    J = np.zeros((n, m))

    # Row 0: F₀ = f(x)
    # ∂F₀/∂x = f'(x)
    J[0, 0] = evaluate_derivative(f_coeffs, x, order=1)

    # Rows 1 to k: Fᵢ = sum_{j=1}^{i} f^(i-j+1)(x) · vⱼ
    for i in range(1, num_deflations + 1):
        # ∂Fᵢ/∂x = sum_{j=1}^{i} f^(i-j+2)(x) · vⱼ
        dFi_dx = 0.0
        for j in range(1, i + 1):
            deriv_order = i - j + 2
            df_val = evaluate_derivative(f_coeffs, x, order=deriv_order)
            dFi_dx += df_val * v[j-1]
        J[i, 0] = dFi_dx

        # ∂Fᵢ/∂vⱼ = f^(i-j+1)(x) for j = 1, ..., i
        for j in range(1, i + 1):
            deriv_order = i - j + 1
            df_val = evaluate_derivative(f_coeffs, x, order=deriv_order)
            J[i, j] = df_val

    # Last row: ||v||² - 1
    if num_deflations > 0:
        for i in range(num_deflations):
            # ∂(||v||² - 1)/∂vᵢ = 2vᵢ
            J[num_deflations + 1, i+1] = 2.0 * v[i]

    return J

def standard_newton(f_coeffs, x0, max_iter=100, tol=1e-15):
    """Standard Newton iteration."""
    x = x0
    for i in range(max_iter):
        f_val = evaluate_poly(f_coeffs, x)
        df_val = evaluate_derivative(f_coeffs, x, 1)
        
        step = f_val / df_val
        
        if i < 10 or i % 10 == 0:
            print(f"  Iter {i}: x = {x:.17f}, |f(x)| = {abs(f_val):.6e}, step = {step:.6e}")
        
        x = x - step
        
        if abs(step) < tol:
            print(f"  Converged at iteration {i+1}")
            return x, i+1
    
    print(f"  Max iterations reached")
    return x, max_iter

def estimate_multiplicity(f_coeffs, x, max_mult=10, tol=1e-10):
    """
    Estimate multiplicity by finding first non-zero derivative.

    At a root of multiplicity m:
      f(x) ≈ 0, f'(x) ≈ 0, ..., f^(m-1)(x) ≈ 0, f^(m)(x) ≠ 0
    """
    for m in range(1, max_mult + 1):
        deriv_val = evaluate_derivative(f_coeffs, x, order=m)
        if abs(deriv_val) > tol:
            return m
    return max_mult

def modified_newton(f_coeffs, x0, m, max_iter=100, tol=1e-15):
    """Modified Newton with multiplicity m."""
    x = x0
    for i in range(max_iter):
        f_val = evaluate_poly(f_coeffs, x)
        df_val = evaluate_derivative(f_coeffs, x, 1)

        if abs(df_val) < 1e-100:
            print(f"  Derivative too small at iteration {i}, stopping")
            return x, i

        step = m * f_val / df_val

        if i < 10 or i % 10 == 0:
            print(f"  Iter {i}: x = {x:.17f}, |f(x)| = {abs(f_val):.6e}, step = {step:.6e}")

        # Check convergence BEFORE taking the step
        # If f(x) is at machine precision, we're done
        if abs(f_val) < 1e-15:
            print(f"  Converged (residual) at iteration {i}")
            return x, i

        x = x - step

        if abs(step) < tol:
            print(f"  Converged (step size) at iteration {i+1}")
            return x, i+1

    print(f"  Max iterations reached")
    return x, max_iter

def adaptive_modified_newton(f_coeffs, x0, max_iter=100, tol=1e-15):
    """
    Modified Newton with adaptive multiplicity estimation.

    At each iteration, estimate multiplicity using high-order derivatives.
    """
    x = x0
    for i in range(max_iter):
        f_val = evaluate_poly(f_coeffs, x)
        df_val = evaluate_derivative(f_coeffs, x, 1)

        if abs(df_val) < 1e-100:
            print(f"  Derivative too small at iteration {i}, stopping")
            return x, i

        # Estimate multiplicity at current point
        m_est = estimate_multiplicity(f_coeffs, x, max_mult=10, tol=1e-8)

        step = m_est * f_val / df_val

        if i < 10 or i % 10 == 0:
            print(f"  Iter {i}: x = {x:.17f}, |f(x)| = {abs(f_val):.6e}, m_est = {m_est}, step = {step:.6e}")

        x = x - step

        if abs(step) < tol:
            print(f"  Converged at iteration {i+1}")
            m_final = estimate_multiplicity(f_coeffs, x, max_mult=10, tol=1e-8)
            print(f"  Final estimated multiplicity: {m_final}")
            return x, i+1

    print(f"  Max iterations reached")
    m_final = estimate_multiplicity(f_coeffs, x, max_mult=10, tol=1e-8)
    print(f"  Final estimated multiplicity: {m_final}")
    return x, max_iter

def iterative_deflation_newton(f_coeffs, x0, max_deflations=10, max_iter=100, tol=1e-15, singular_tol=1e-8):
    """
    Iterative deflation using augmented system and Newton's method.

    Algorithm:
    1. Start with stage=0 (original system)
    2. Apply Newton to current augmented system
    3. Check if Jacobian is singular using SVD
    4. If singular, increment stage (add new variable and equation)
    5. Repeat until Jacobian is non-singular
    """
    stage = 0
    z = np.array([x0])  # Initial: just x

    print("\n=== ITERATIVE DEFLATION WITH AUGMENTED SYSTEM ===\n")

    for deflation_round in range(max_deflations):
        print(f"--- Deflation Stage {stage} ---")
        print(f"Variables: {len(z)} (x" + (f", v1...v{stage}" if stage > 0 else "") + ")")
        print(f"Equations: {len(z) + (1 if stage > 0 else 0)}")

        # Apply Newton's method to current augmented system
        for iter in range(max_iter):
            F = augmented_system_residual(f_coeffs, z, stage)
            J = augmented_system_jacobian(f_coeffs, z, stage)

            # Check Jacobian rank using SVD
            _, s, _ = np.linalg.svd(J)
            # Use ABSOLUTE tolerance for rank determination
            rank = np.sum(s > singular_tol)
            cond = s[0] / s[-1] if s[-1] > 1e-100 else np.inf

            if iter < 5 or iter % 10 == 0:
                print(f"  Iter {iter}: |F| = {np.linalg.norm(F):.6e}, " +
                      f"rank = {rank}/{len(s)}, cond = {cond:.2e}")

            # Solve J * step = -F using least-squares (handles singular J)
            step, residuals, rank_solve, s_solve = np.linalg.lstsq(J, -F, rcond=None)

            z = z + step

            if np.linalg.norm(step) < tol:
                print(f"  Converged at iteration {iter+1}")
                break

        # Check if Jacobian is non-singular
        F_final = augmented_system_residual(f_coeffs, z, stage)
        J_final = augmented_system_jacobian(f_coeffs, z, stage)
        _, s, _ = np.linalg.svd(J_final)
        # Use ABSOLUTE tolerance for rank determination
        rank = np.sum(s > singular_tol)

        print(f"\nAfter Newton convergence:")
        print(f"  x = {z[0]:.17f}")
        print(f"  |F| = {np.linalg.norm(F_final):.6e}")
        print(f"  Jacobian rank: {rank}/{len(s)}")
        print(f"  Singular values: {s}")

        if rank == len(s):
            print(f"  ✓ Jacobian is FULL RANK - deflation complete!\n")
            break
        else:
            print(f"  ✗ Jacobian is RANK DEFICIENT - need more deflation\n")
            # Add new variable v_{stage+1}
            stage += 1
            z = np.append(z, 1.0 / np.sqrt(stage))  # Initialize new v with normalized value

    print(f"Total deflation stages: {stage}")
    print(f"Detected multiplicity: {stage + 1}\n")

    return z[0], stage

def main():
    import sys

    # Determine which test to run
    if len(sys.argv) > 1 and sys.argv[1] == "simple":
        # Test polynomial: p(x) = (x - 0.6)^6
        print("=== DEFLATION TEST FOR PURE MULTIPLE ROOT ===\n")
        print("Polynomial: p(x) = (x - 0.6)^6")
        print("Expected root: x = 0.6 (multiplicity 6)\n")

        # Expanded form of (x - 0.6)^6
        # Using binomial expansion or numpy
        from numpy.polynomial import Polynomial
        p = Polynomial([-0.6, 1.0]) ** 6
        coeffs = p.coef.tolist()
        print(f"Coefficients (power basis): {coeffs}\n")

        x0 = 0.59
        true_root = 0.6
        expected_mult = 6
    else:
        # Test polynomial: p(x) = (x - 0.2)(x - 0.6)^6
        print("=== DEFLATION TEST FOR MULTIPLE ROOTS ===\n")
        print("Polynomial: p(x) = (x - 0.2)(x - 0.6)^6")
        print("Expected roots:")
        print("  x = 0.2 (multiplicity 1)")
        print("  x = 0.6 (multiplicity 6)\n")

        # Expanded form
        coeffs = [
            -0.0093312,   # x^0
            0.139968,     # x^1
            -0.85536,     # x^2
            2.808,        # x^3
            -5.4,         # x^4
            6.12,         # x^5
            -3.8,         # x^6
            1.0           # x^7
        ]

        x0 = 0.59
        true_root = 0.6
        expected_mult = 6

    # Test 1: Standard Newton (will fail)
    print("=== TEST 1: STANDARD NEWTON (NO DEFLATION) ===")
    x1, iters1 = standard_newton(coeffs, x0)
    error1 = abs(x1 - true_root)
    print(f"Final x: {x1:.17f}")
    print(f"Error: {error1:.6e}\n")

    # Test 2: Modified Newton with m=6 (correct multiplicity)
    print("=== TEST 2: MODIFIED NEWTON WITH m=6 (CORRECT MULTIPLICITY) ===")
    x2, iters2 = modified_newton(coeffs, x0, m=6)
    error2 = abs(x2 - true_root)
    print(f"Final x: {x2:.17f}")
    print(f"Error: {error2:.6e}\n")

    # Test 3: Adaptive Modified Newton (estimates multiplicity)
    print("=== TEST 3: ADAPTIVE MODIFIED NEWTON (AUTO-DETECT MULTIPLICITY) ===")
    x3, iters3 = adaptive_modified_newton(coeffs, x0)
    error3 = abs(x3 - true_root)
    print(f"Final x: {x3:.17f}")
    print(f"Error: {error3:.6e}\n")

    # Test 4: Iterative deflation with augmented system
    print("=== TEST 4: ITERATIVE DEFLATION (AUGMENTED SYSTEM) ===")
    x4, stages4 = iterative_deflation_newton(coeffs, x0)
    error4 = abs(x4 - true_root)
    print(f"Final x: {x4:.17f}")
    print(f"Error: {error4:.6e}\n")

    # Summary
    print("=== SUMMARY ===")
    print(f"Test 1 (Standard Newton):           error = {error1:.6e}")
    print(f"Test 2 (Modified Newton m=6):       error = {error2:.6e}")
    print(f"Test 3 (Adaptive Modified Newton):  error = {error3:.6e}")
    print(f"Test 4 (Iterative Deflation):       error = {error4:.6e}, stages = {stages4}")
    print(f"\nExpected multiplicity: {expected_mult}")
    print(f"Detected multiplicity (deflation): {stages4 + 1}")

if __name__ == "__main__":
    main()

