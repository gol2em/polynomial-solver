# Univariate (1D) Polynomial Equation Solvers

## INTRODUCTION

This document provides a comprehensive reference for algorithms that solve univariate polynomial equations of the form:

**p(x) = aₙxⁿ + aₙ₋₁xⁿ⁻¹ + ... + a₁x + a₀ = 0**

### Algorithm Categories

1. **Closed-Form Solutions** (Degrees 1-4): Exact algebraic formulas
2. **Root Counting Methods**: Determine number of roots in intervals
3. **Root Isolation Methods**: Separate roots into distinct intervals
4. **Root Bounds**: Estimate magnitude of roots
5. **Numerical Methods**: Iterative approximation techniques
6. **Interval Arithmetic Methods**: Guaranteed bounds on roots
7. **Algebraic Methods**: Symbolic manipulation and factorization

### Wide Application Markers

Algorithms marked with **[WIDELY USED]** are implemented in major libraries and have broad practical applications across:
- Scientific computing (NumPy, SciPy, MATLAB)
- Computer algebra systems (Mathematica, Maple, SageMath)
- Engineering software (MATLAB, Octave)
- Graphics and CAD systems
- Numerical analysis libraries (GSL, LAPACK)

### Choosing an Algorithm

**For exact solutions (degrees 1-4)**:
- Use closed-form formulas (quadratic, cubic, quartic)
- Widely implemented in all major software

**For numerical approximation**:
- **Single root**: Newton-Raphson (fast, requires derivative), Bisection (robust, slow)
- **All roots**: Companion Matrix Method (industry standard), Aberth Method (parallel)
- **High precision**: MPSolve (arbitrary precision), Interval Newton (guaranteed bounds)

**For symbolic computation**:
- Sturm's Theorem (exact root counting)
- Square-free factorization (multiple roots)
- Thom Encoding (symbolic root comparison)

**For root isolation**:
- Vincent's Theorem (continued fractions)
- Real Roots Isolation Method (geometric bounds)
- Bernstein Subdivision (convex hull property)

---

## CLOSED-FORM SOLUTIONS (Degrees 1-4)

### 1. Linear Equation (Degree 1)
- **Formula**: ax + b = 0
- **Solution**: x = -b/a
- **Complexity**: O(1)
- **Conditions**: a ≠ 0

### 2. Quadratic Equation (Degree 2) **[WIDELY USED]**
- **Formula**: ax² + bx + c = 0
- **Solutions**:
  ```
  x = (-b ± √(b² - 4ac)) / (2a)
  ```
- **Discriminant**: Δ = b² - 4ac
  - Δ > 0: Two distinct real roots
  - Δ = 0: One repeated real root
  - Δ < 0: Two complex conjugate roots
- **Complexity**: O(1)
- **Applications**:
  - Physics (projectile motion, optics)
  - Engineering (circuit analysis, control systems)
  - Computer graphics (ray-sphere intersection)
  - Finance (option pricing)
- **Historical**: Known since ancient Babylon (~2000 BC)
- **Key References**:
  - Euclid's Elements (300 BC)
  - Al-Khwarizmi (820 AD): "The Compendious Book on Calculation by Completion and Balancing"

### 3. Cubic Equation (Degree 3)
- **General Form**: ax³ + bx² + cx + d = 0
- **Depressed Form**: t³ + pt + q = 0 (via substitution x = t - b/3a)

#### 3a. Cardano's Formula
- **Core Algorithm**:
  1. Reduce to depressed cubic: t³ + pt + q = 0
  2. Compute discriminant: Δ = -4p³ - 27q²
  3. If Δ > 0 (three real roots): Use trigonometric form
  4. If Δ ≤ 0 (one real root): Use Cardano's formula with cube roots
- **Formula**:
  ```
  t = ∛(-q/2 + √(q²/4 + p³/27)) + ∛(-q/2 - √(q²/4 + p³/27))
  ```
- **Discriminant**: Δ = -4p³ - 27q²
  - Δ > 0: Three distinct real roots (casus irreducibilis)
  - Δ = 0: Multiple root
  - Δ < 0: One real root, two complex conjugate roots
- **Historical**: Scipione del Ferro (1515), Tartaglia (1535), Cardano (1545)
- **Key References**:
  - Cardano, G. (1545): "Ars Magna"
  - Bombelli, R. (1572): "L'Algebra" (complex numbers)

#### 3b. Trigonometric Solution (Three Real Roots)
- **When**: Δ > 0 (casus irreducibilis)
- **Formula**:
  ```
  t_k = 2√(-p/3) cos((1/3)arccos(3q/(2p)√(-3/p)) - 2πk/3)  for k = 0,1,2
  ```
- **Advantage**: Avoids complex arithmetic when all roots are real

#### 3c. Vieta's Substitution
- **Method**: Substitute t = w - p/(3w)
- **Result**: Reduces to quadratic in w³

### 4. Quartic Equation (Degree 4)
- **General Form**: ax⁴ + bx³ + cx² + dx + e = 0
- **Depressed Form**: y⁴ + py² + qy + r = 0 (via substitution x = y - b/4a)

#### 4a. Ferrari's Method
- **Core Algorithm**:
  1. Reduce to depressed quartic: y⁴ + py² + qy + r = 0
  2. Rewrite as: (y² + p/2)² = -qy - r + p²/4
  3. Add (y + m)² to both sides to make right side a perfect square
  4. Solve resolvent cubic for m
  5. Factor into two quadratics
  6. Solve each quadratic
- **Resolvent Cubic**: 8m³ + 8pm² + (2p² - 8r)m - q² = 0
- **Historical**: Lodovico Ferrari (1540), published by Cardano (1545)

#### 4b. Descartes' Resolvent
- **Formula**: (p² + c)² - (d/p)² = 4e
- **Result**: Cubic in p²

#### 4c. Biquadratic Case
- **When**: a₃ = a₁ = 0 (form: ax⁴ + cx² + e = 0)
- **Method**: Substitute z = x², solve quadratic, then x = ±√z

- **Key References**:
  - Ferrari, L. (1545): Solution in Cardano's "Ars Magna"
  - Euler, L. (1738): Alternative methods

### 5. Quintic and Higher (Degree ≥ 5)
- **Abel-Ruffini Theorem** (1824): No general algebraic solution exists
- **Galois Theory**: Explains why degrees ≥ 5 are not solvable by radicals
- **Special Cases**: Some quintics are solvable (e.g., x⁵ - x - 1)
- **Key References**:
  - Abel, N.H. (1824): "Mémoire sur les équations algébriques"
  - Galois, É. (1832): Galois theory (published posthumously 1846)

## ROOT COUNTING AND ISOLATION METHODS

### 6. Descartes' Rule of Signs
- **Type**: Root counting (bounds)
- **Description**: Counts positive/negative real roots by examining coefficient sign changes
- **Core Algorithm**:
  1. Count sign changes in coefficient sequence (ignoring zeros)
  2. Number of positive roots = sign changes - 2k (for some k ≥ 0)
  3. For negative roots: apply to p(-x)
- **Formula**: 
  - Positive roots ≤ number of sign changes
  - Difference is even
- **Example**: p(x) = x³ + x² - x - 1 has signs (+,+,-,-) → 1 sign change → 1 positive root
- **Complexity**: O(n) for degree n
- **Historical**: René Descartes (1637) in "La Géométrie"
- **Key References**:
  - Descartes, R. (1637): "La Géométrie"
  - Fourier, J. (1820): Extensions and refinements

### 7. Sturm's Theorem **[WIDELY USED]**
- **Type**: Exact root counting in intervals
- **Description**: Counts real roots in interval [a,b] using Sturm sequence
- **Core Algorithm**:
  1. Compute Sturm sequence: P₀ = p, P₁ = p', Pᵢ₊₁ = -rem(Pᵢ₋₁, Pᵢ)
  2. Evaluate sequence at a and b
  3. Count sign changes: V(a) and V(b)
  4. Number of roots in (a,b] = V(a) - V(b)
- **Sturm Sequence**: p₀, p₁, ..., pₘ where:
  - p₀ = p(x)
  - p₁ = p'(x)
  - pᵢ₊₁ = -remainder(pᵢ₋₁ ÷ pᵢ)
- **Complexity**: O(n²) for degree n polynomial
- **Advantages**:
  - Exact count (not just bounds)
  - Works on any interval
  - Provides root isolation
- **Applications**:
  - Root isolation algorithms
  - Computer algebra systems (SymPy, SageMath, Mathematica)
  - Symbolic computation
  - Polynomial analysis and theorem proving
- **Historical**: Jacques Charles François Sturm (1829)
- **Key References**:
  - Sturm, J.C.F. (1829): "Mémoire sur la résolution des équations numériques"
  - Sylvester, J.J. (1853): Extensions
- **C/Python Implementations**:
  - **SymPy** (Python): sturm() function
  - **SageMath** (Python): polynomial.sturm_sequence()
  - **Macaulay2** (C++): Sturm sequence computation

### 8. Budan's Theorem / Budan-Fourier Theorem
- **Type**: Root counting in intervals
- **Description**: Generalization of Descartes' rule using derivatives
- **Core Algorithm**:
  1. Compute derivatives: p, p', p'', ..., p⁽ⁿ⁾
  2. Count sign variations at a: V(a)
  3. Count sign variations at b: V(b)
  4. Number of roots in (a,b) ≤ V(a) - V(b), difference is even
- **Relation to Sturm**: Similar but uses derivatives instead of remainders
- **Historical**: F.D. Budan (1807), J. Fourier (1820)
- **Key References**:
  - Budan, F.D. (1807): "Nouvelle méthode pour la résolution des équations numériques"
  - Fourier, J. (1820): Extensions
  - Akritas, A.G. (1982): "Reflections on a pair of theorems by Budan and Fourier"

### 9. Vincent's Theorem / Continued Fraction Method
- **Type**: Root isolation
- **Description**: Uses continued fraction expansion for root isolation
- **Core Algorithm**:
  1. Apply linear fractional transformations
  2. Use Descartes' rule to count roots
  3. Subdivide intervals using continued fractions
  4. Isolate each root in a unique interval
- **Complexity**: Efficient for sparse polynomials
- **Historical**: A.J.H. Vincent (1836)
- **Modern Use**: Basis for fastest real root isolation algorithms
- **Key References**:
  - Vincent, A.J.H. (1836): "Sur la résolution des équations numériques"
  - Akritas, A.G., Strzeboński, A.W., Vigklas, P.S. (2008): "Improving the performance of the continued fractions method"
- **C/Python Implementations**:
  - **CGAL** (C++): Algebraic kernel with Vincent's method
  - **mpsolve** (C): Uses Vincent-based isolation

## ROOT BOUNDS

### 10. Cauchy Bound
- **Type**: Root magnitude bound
- **Formula**: All roots satisfy |x| < 1 + max(|aₙ₋₁|, |aₙ₋₂|, ..., |a₀|) / |aₙ|
- **Simpler Form**: |x| < 1 + max|aᵢ/aₙ|
- **Use**: Provides initial interval for root finding
- **Historical**: Augustin-Louis Cauchy (1829)

### 11. Lagrange Bound
- **Type**: Root magnitude bound
- **Formula**: |x| < max(1, Σ|aᵢ/aₙ|) for i = 0 to n-1
- **Tighter**: Often better than Cauchy bound
- **Historical**: Joseph-Louis Lagrange (1798)

### 12. Fujiwara Bound
- **Type**: Root magnitude bound
- **Formula**: |x| < 2 · max(|aₙ₋₁/aₙ|^(1/1), |aₙ₋₂/aₙ|^(1/2), ..., |a₀/aₙ|^(1/n))
- **Advantage**: Often tighter than Cauchy and Lagrange

### 13. Kojima's Bound
- **Type**: Root magnitude bound (tighter than Cauchy)
- **Description**: Geometric property-based bound for polynomial roots
- **Advantage**: More numerically friendly than Fujiwara's bound while still tight
- **Use**: Used in Real Roots Isolation Method
- **Comparison**: Tighter than Cauchy bound by several magnitudes in most cases
- **Reference**: Used in ZhepeiWang/Root-Finder (2019)

## NUMERICAL METHODS

### 13. Newton-Raphson Method **[WIDELY USED]**
- **Type**: Numerical (iterative)
- **Description**: Iterative method using derivative
- **Core Algorithm**:
  1. Start with initial guess x₀
  2. Iterate: xₙ₊₁ = xₙ - p(xₙ)/p'(xₙ)
  3. Stop when |xₙ₊₁ - xₙ| < ε
- **Convergence**: Quadratic (near simple roots)
- **Complexity**: O(n) per iteration for degree n
- **Advantages**: Fast convergence when close to root
- **Disadvantages**: Requires good initial guess, may diverge
- **Applications**:
  - Numerical analysis (most common root-finding method)
  - Optimization algorithms
  - Scientific computing (SciPy, MATLAB, Octave)
  - Engineering simulations
  - Machine learning (gradient descent is a variant)
- **Historical**: Isaac Newton (1669), Joseph Raphson (1690)
- **Key References**:
  - Newton, I. (1669): "De analysi per aequationes numero terminorum infinitas"
  - Raphson, J. (1690): "Analysis aequationum universalis"
- **C/Python Implementations**:
  - **SciPy** (Python): scipy.optimize.newton()
  - **GSL** (C): gsl_root_fsolver_newton
  - **NumPy** (Python): Can be implemented easily
  - **MATLAB**: fzero() uses variant

### 14. Bisection Method
- **Type**: Numerical (bracketing)
- **Description**: Repeatedly halves interval containing root
- **Core Algorithm**:
  1. Start with [a,b] where p(a)·p(b) < 0
  2. Compute midpoint c = (a+b)/2
  3. If p(c)·p(a) < 0, set b = c; else set a = c
  4. Repeat until |b-a| < ε
- **Convergence**: Linear (guaranteed)
- **Complexity**: O(log(1/ε)) iterations
- **Advantages**: Always converges, simple
- **Disadvantages**: Slow convergence
- **C/Python Implementations**:
  - **SciPy** (Python): scipy.optimize.bisect()
  - **GSL** (C): gsl_root_fsolver_bisection

### 15. Secant Method
- **Type**: Numerical (iterative)
- **Description**: Like Newton-Raphson but approximates derivative
- **Core Algorithm**:
  1. Start with two guesses x₀, x₁
  2. Iterate: xₙ₊₁ = xₙ - p(xₙ)·(xₙ - xₙ₋₁)/(p(xₙ) - p(xₙ₋₁))
  3. Stop when |xₙ₊₁ - xₙ| < ε
- **Convergence**: Superlinear (φ ≈ 1.618)
- **Advantages**: No derivative needed
- **C/Python Implementations**:
  - **SciPy** (Python): scipy.optimize.newton() with fprime=None
  - **GSL** (C): gsl_root_fsolver_secant

### 16. Laguerre's Method
- **Type**: Numerical (iterative, polynomial-specific)
- **Description**: Uses first and second derivatives, designed for polynomials
- **Core Algorithm**:
  1. Compute p(x), p'(x), p''(x)
  2. Iterate: xₙ₊₁ = xₙ - n·p(xₙ) / (p'(xₙ) ± √H(xₙ))
     where H(x) = (n-1)[(n-1)(p'(x))² - n·p(x)·p''(x)]
  3. Choose sign to maximize denominator
- **Convergence**: Cubic (very fast)
- **Advantages**:
  - Converges to a root from almost any starting point
  - Works well for multiple roots
  - Finds complex roots
- **Historical**: Edmond Laguerre (1880)
- **Key References**:
  - Laguerre, E. (1880): "Sur une méthode pour obtenir par approximation les racines"
- **C/Python Implementations**:
  - **NumPy** (Python): numpy.polynomial.polynomial.polyroots() uses companion matrix
  - **mpmath** (Python): polyroots() uses Laguerre-like methods

### 17. Jenkins-Traub Algorithm
- **Type**: Numerical (iterative, polynomial-specific)
- **Description**: Three-stage algorithm for finding all polynomial roots
- **Core Algorithm**:
  1. Stage 1: No-shift process (5 iterations)
  2. Stage 2: Fixed-shift process (variable iterations)
  3. Stage 3: Variable-shift process (converges to root)
  4. Deflate polynomial and repeat
- **Convergence**: Globally convergent, cubic near roots
- **Advantages**:
  - Finds all roots (real and complex)
  - Very reliable
  - Industry standard
- **Complexity**: O(n²) for degree n
- **Historical**: M.A. Jenkins, J.F. Traub (1970)
- **Key References**:
  - Jenkins, M.A., Traub, J.F. (1970): "A Three-Stage Algorithm for Real Polynomials Using Quadratic Iteration"
  - Jenkins, M.A., Traub, J.F. (1970): "A Three-Stage Variable-Shift Iteration for Polynomial Zeros"
- **C/Python Implementations**:
  - **NumPy** (Python): numpy.roots() (uses eigenvalues of companion matrix, similar approach)
  - **FORTRAN**: Original RPOLY/CPOLY implementations
  - **Netlib**: Jenkins-Traub FORTRAN code

### 18. Aberth Method (Aberth-Ehrlich Method) **[WIDELY USED]**
- **Type**: Numerical (simultaneous iteration)
- **Description**: Finds all roots simultaneously using parallel iteration
- **Core Algorithm**:
  1. Start with n initial approximations z₁, ..., zₙ
  2. Iterate simultaneously for all i:
     zᵢ⁽ᵏ⁺¹⁾ = zᵢ⁽ᵏ⁾ - 1 / (p'(zᵢ)/p(zᵢ) - Σⱼ≠ᵢ 1/(zᵢ - zⱼ))
  3. Stop when all |zᵢ⁽ᵏ⁺¹⁾ - zᵢ⁽ᵏ⁾| < ε
- **Convergence**: Cubic (very fast)
- **Advantages**:
  - Finds all roots at once
  - Parallel-friendly (GPU implementations exist)
  - Works for complex roots
  - Cubic convergence (faster than Newton-Raphson)
- **Applications**:
  - **MPSolve** (primary algorithm, arbitrary precision)
  - High-performance polynomial solving
  - Parallel computing applications
  - Complex root finding
- **Historical**: Oliver Aberth (1973)
- **Key References**:
  - Aberth, O. (1973): "Iteration methods for finding all zeros of a polynomial simultaneously"
  - Ehrlich, L.W. (1967): "A modified Newton method for polynomials"
  - Bini, D.A. (1996): "Numerical computation of polynomial zeros by means of Aberth's method"
- **C/Python Implementations**:
  - **MPSolve** (C): Uses Aberth method as primary algorithm
  - **mpmath** (Python): Can use Aberth-like methods

### 19. Durand-Kerner Method (Weierstrass Method)
- **Type**: Numerical (simultaneous iteration)
- **Description**: Simpler simultaneous iteration method
- **Core Algorithm**:
  1. Start with n initial approximations z₁, ..., zₙ
  2. Iterate simultaneously for all i:
     zᵢ⁽ᵏ⁺¹⁾ = zᵢ⁽ᵏ⁾ - p(zᵢ⁽ᵏ⁾) / ∏ⱼ≠ᵢ (zᵢ⁽ᵏ⁾ - zⱼ⁽ᵏ⁾)
  3. Stop when all |zᵢ⁽ᵏ⁺¹⁾ - zᵢ⁽ᵏ⁾| < ε
- **Convergence**: Linear to quadratic
- **Advantages**: Simple, finds all roots
- **Historical**: Karl Weierstrass (1891), Durand (1960), Kerner (1966)
- **Key References**:
  - Weierstrass, K. (1891): "Neuer Beweis des Satzes, dass jede ganze rationale Function"
  - Durand, E. (1960): "Solutions Numériques des Équations Algébriques"
- **C/Python Implementations**:
  - **mpmath** (Python): Can implement easily
  - Various educational implementations

### 20. Companion Matrix Method **[WIDELY USED]**
- **Type**: Numerical (eigenvalue-based)
- **Description**: Converts polynomial root finding to eigenvalue problem
- **Core Algorithm**:
  1. Construct companion matrix C for polynomial p(x) = Σ aᵢxⁱ:
     ```
     C = [0   0   ...  0   -a₀/aₙ]
         [1   0   ...  0   -a₁/aₙ]
         [0   1   ...  0   -a₂/aₙ]
         [⋮   ⋮   ⋱   ⋮      ⋮   ]
         [0   0   ...  1   -aₙ₋₁/aₙ]
     ```
  2. Compute eigenvalues of C using QR algorithm or similar
  3. Eigenvalues = roots of polynomial
- **Complexity**: O(n³) using standard eigenvalue algorithms
- **Advantages**:
  - Finds all roots (real and complex)
  - Leverages mature eigenvalue solvers (LAPACK)
  - Numerically stable
  - Industry standard method
- **Disadvantages**: Slower than specialized polynomial methods
- **Applications**:
  - **NumPy/SciPy** (default method for numpy.roots())
  - **MATLAB** (roots() function)
  - Scientific computing libraries
  - Control theory (characteristic polynomial roots)
  - Signal processing (filter design)
- **Key References**:
  - Edelman, A., Murakami, H. (1995): "Polynomial roots from companion matrix eigenvalues"
- **C/Python Implementations**:
  - **NumPy** (Python): numpy.roots() uses this method
  - **LAPACK** (C/Fortran): Eigenvalue routines (DGEEV, ZGEEV)
  - **Eigen** (C++): eigenvalues() method
  - **MATLAB**: roots() function

## INTERVAL ARITHMETIC METHODS

### 21. Interval Newton Method
- **Type**: Numerical (interval arithmetic)
- **Description**: Newton method with interval arithmetic for guaranteed bounds
- **Core Algorithm**:
  1. Start with interval [a,b] containing root
  2. Compute interval extension: N([a,b]) = m - p(m)/p'([a,b])
     where m = midpoint of [a,b]
  3. Intersect N([a,b]) with [a,b]
  4. Repeat until interval width < ε
- **Advantages**:
  - Provides guaranteed bounds
  - Can prove existence/uniqueness
  - Rigorous error bounds
- **C/Python Implementations**:
  - **MPFI** (C): Multiple precision interval arithmetic
  - **pyinterval** (Python): Interval arithmetic
  - **mpmath** (Python): Interval arithmetic support

### 22. Bernstein Basis Subdivision
- **Type**: Interval arithmetic (geometric)
- **Description**: Uses Bernstein polynomial representation and convex hull property
- **Core Algorithm**:
  1. Convert polynomial to Bernstein basis on [0,1]
  2. Check sign variations in Bernstein coefficients (Descartes' rule)
  3. If all coefficients same sign: no roots
  4. If one sign change: exactly one root
  5. Otherwise: subdivide at midpoint using de Casteljau algorithm
  6. Recursively process sub-intervals
- **Advantages**:
  - Geometrically intuitive
  - Numerically stable subdivision
  - Natural for computer graphics applications
- **Complexity**: Adaptive, depends on root distribution
- **Key References**:
  - Farouki, R.T. (2012): "The Bernstein polynomial basis: A centennial retrospective"
  - Spencer, J. (1994): "Bernstein polynomials and learning theory"
- **C/Python Implementations**:
  - **CGAL** (C++): Polynomial root isolation with Bernstein
  - Custom implementations in CAD/graphics libraries

## SPECIAL PROPERTIES AND FORMULAS

### 23. Vieta's Formulas
- **Type**: Algebraic relations
- **Description**: Relates roots to coefficients
- **Formulas** (for polynomial aₙxⁿ + aₙ₋₁xⁿ⁻¹ + ... + a₁x + a₀):
  - Sum of roots: r₁ + r₂ + ... + rₙ = -aₙ₋₁/aₙ
  - Sum of products (pairs): Σᵢ<ⱼ rᵢrⱼ = aₙ₋₂/aₙ
  - Product of all roots: r₁·r₂·...·rₙ = (-1)ⁿ·a₀/aₙ
- **Use**:
  - Verify computed roots
  - Find remaining roots when some are known
  - Polynomial reconstruction from roots
- **Historical**: François Viète (1591)
- **Key References**:
  - Viète, F. (1591): "In artem analyticem isagoge"

### 24. Discriminant
- **Type**: Algebraic invariant
- **Description**: Determines nature of roots without solving
- **Formulas**:
  - **Quadratic** (ax² + bx + c): Δ = b² - 4ac
  - **Cubic** (x³ + px + q): Δ = -4p³ - 27q²
  - **Quartic**: More complex formula involving all coefficients
- **Interpretation**:
  - Δ > 0: Distinct real roots (for cubic/quartic: specific patterns)
  - Δ = 0: Multiple roots
  - Δ < 0: Complex roots (for cubic: one real, two complex)
- **Use**: Determine root type before computation

### 25. Multiple Root Detection
- **Type**: Algebraic (GCD-based)
- **Description**: Detect and remove multiple roots
- **Core Algorithm**:
  1. Compute GCD of p(x) and p'(x)
  2. If GCD = 1: no multiple roots
  3. If GCD = g(x): g(x) contains all multiple roots
  4. Square-free part: p(x) / g(x)
- **Use**:
  - Simplify root finding (work with square-free polynomial)
  - Improve numerical stability
- **C/Python Implementations**:
  - **SymPy** (Python): sqf() for square-free factorization
  - **SageMath** (Python): squarefree_decomposition()
  - **Singular** (C++): GCD and square-free factorization

## SUMMARY TABLE

| Method | Type | Degree | Complexity | Exact/Numerical | Key Advantage |
|--------|------|--------|------------|-----------------|---------------|
| Linear Formula | Closed-form | 1 | O(1) | Exact | Trivial |
| Quadratic Formula | Closed-form | 2 | O(1) | Exact | Simple, well-known |
| Cardano's Formula | Closed-form | 3 | O(1) | Exact | Algebraic solution |
| Ferrari's Method | Closed-form | 4 | O(1) | Exact | Algebraic solution |
| Descartes' Rule | Counting | Any | O(n) | Bounds | Simple, fast bounds |
| Sturm's Theorem | Counting | Any | O(n²) | Exact count | Exact root count |
| Vincent's Method | Isolation | Any | Adaptive | Exact intervals | Fast isolation |
| Real Roots Isolation | Geometric | Any | Adaptive | Exact intervals | Faster than Jenkins-Traub |
| Kojima's Bound | Bounds | Any | O(n) | Tighter bounds | Better than Cauchy |
| Newton-Raphson | Numerical | Any | O(n) per iter | Numerical | Fast convergence |
| Laguerre's Method | Numerical | Any | O(n) per iter | Numerical | Robust for polynomials |
| Jenkins-Traub | Numerical | Any | O(n²) | Numerical | Industry standard |
| Aberth Method | Numerical | Any | O(n²) | Numerical | Parallel, all roots |
| Companion Matrix | Numerical | Any | O(n³) | Numerical | Leverages LAPACK |
| Power Iteration | Numerical | Any | O(n) per iter | Numerical | Dominant eigenvalue |
| Bernstein Subdivision | Interval | Any | Adaptive | Guaranteed bounds | Geometric, stable |
| Thom Encoding | Symbolic | Any | Polynomial | Exact | Symbolic root comparison |
| Kronecker Factorization | Algebraic | Any | Varies | Exact | Integer factorization |

## KEY REFERENCES

### Historical Works
1. **Cardano, G.** (1545): "Ars Magna" - Cubic and quartic solutions
2. **Viète, F.** (1591): "In artem analyticem isagoge" - Vieta's formulas
3. **Descartes, R.** (1637): "La Géométrie" - Rule of signs
4. **Newton, I.** (1669): "De analysi" - Newton's method
5. **Sturm, J.C.F.** (1829): "Mémoire sur la résolution des équations numériques"
6. **Abel, N.H.** (1824): Impossibility of quintic formula
7. **Galois, É.** (1832): Galois theory

### Modern References
1. **Pan, V.Y.** (1997): "Solving a Polynomial Equation: Some History and Recent Progress"
2. **McNamee, J.M.** (2007): "Numerical Methods for Roots of Polynomials, Part I"
3. **McNamee, J.M., Pan, V.Y.** (2013): "Numerical Methods for Roots of Polynomials, Part II"
4. **Bini, D.A., Fiorentino, G.** (2000): "Design, analysis, and implementation of a multiprecision polynomial rootfinder"

### Software Documentation
1. **NumPy**: numpy.roots(), numpy.polynomial
2. **SciPy**: scipy.optimize root-finding functions
3. **SymPy**: Symbolic polynomial solving
4. **MPSolve**: High-performance polynomial solver
5. **CGAL**: Computational geometry algorithms library

