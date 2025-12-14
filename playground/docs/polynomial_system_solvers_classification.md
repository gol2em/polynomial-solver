# Classification of 2-Variable Polynomial System Solvers

## INTRODUCTION

This document provides a comprehensive reference for algorithms that solve systems of two polynomial equations in two variables:

**f(x,y) = 0**
**g(x,y) = 0**

### Algorithm Categories

1. **Symbolic Methods**: Exact algebraic computation (Gröbner bases, resultants, triangular decomposition)
2. **Numerical Methods**: Iterative approximation (homotopy continuation, Newton-Raphson, interval arithmetic)
3. **Hybrid Methods**: Combination of symbolic and numerical techniques (CAD, numerical algebraic geometry)

### Wide Application Markers

Algorithms marked with **[WIDELY USED]** are implemented in major computer algebra systems and have broad practical applications across:
- Computer algebra systems (Mathematica, Maple, SageMath, Macaulay2, SINGULAR)
- Numerical libraries (PHCpack, Bertini, HomotopyContinuation.jl)
- Scientific computing (SciPy, NumPy)
- Engineering software (MATLAB)
- Computational geometry (CGAL)

### Choosing an Algorithm

**For exact symbolic solutions**:
- **Small systems, low degree**: Resultant Method (fast, direct)
- **General systems**: Gröbner Basis (F4 algorithm in msolve, Macaulay2, SINGULAR)
- **Parametric solutions**: RUR (Rational Univariate Representation)
- **Geometric problems**: Wu's Method, Regular Chains

**For numerical approximation**:
- **All isolated solutions**: Homotopy Continuation (PHCpack, Bertini, HomotopyContinuation.jl)
- **Local solutions**: Newton-Raphson (fast with good initial guess)
- **Guaranteed bounds**: Interval Arithmetic (subdivision methods)

**For hybrid approaches**:
- **Quantifier elimination**: CAD (Cylindrical Algebraic Decomposition)
- **Positive-dimensional sets**: Numerical Algebraic Geometry

**For special structures**:
- **Bivariate systems**: Hidden Variable Resultant (spoly), Bernstein Subdivision
- **Sparse systems**: Polyhedral homotopy (PHCpack, Bertini)
- **Overdetermined systems**: Least squares, optimization methods

**For subdivision-based methods**:
- **Bernstein basis with convex hull**: ProjectedPolyhedral method (this project), GraphHull method (this project)
- **Interval arithmetic**: MPFI, Arb

---

## POLYNOMIAL BASIS REPRESENTATIONS

Most algorithms work with polynomials in **power basis** (monomial basis): p(x) = a₀ + a₁x + a₂x² + ...

However, other basis representations exist and can be advantageous for certain applications:

### Alternative Polynomial Bases:
1. **Power/Monomial Basis** (standard): {1, x, x², x³, ...}
   - Most common representation
   - Used by most symbolic solvers

2. **Bernstein Basis**: Used in Bézier curves and computer graphics
   - Better numerical stability for certain operations
   - Basis: bᵥ,ₙ(x) = (n choose v) xᵛ(1-x)ⁿ⁻ᵛ for x ∈ [0,1]
   - Applications: CAD/CAM, computer graphics, approximation theory
   - Conversion to/from power basis is possible but can be numerically unstable
   - **Key Properties for Root Finding**:
     * **Convex Hull Property**: Polynomial graph lies within convex hull of control points
     * **Variation Diminishing**: Number of roots ≤ number of sign changes in coefficients
     * **Subdivision**: de Casteljau algorithm provides stable subdivision at any parameter value
     * **Endpoint Interpolation**: b₀ = p(0), bₙ = p(1)
   - **Root Isolation Algorithm**: Convert to Bernstein basis → check sign variations → subdivide if needed
   - Used in specialized subdivision solvers for guaranteed root isolation

3. **Chebyshev Basis**: Optimal for approximation and interpolation
   - First kind: Tₙ(x) defined by Tₙ(cos θ) = cos(nθ)
   - Second kind: Uₙ(x) defined by Uₙ(cos θ) sin θ = sin((n+1)θ)
   - Minimizes interpolation error (Runge's phenomenon)
   - Used in Clenshaw-Curtis quadrature
   - Better conditioning than power basis

4. **Legendre Basis**: Orthogonal polynomials on [-1,1]
   - Used in spectral methods
   - Good for numerical integration

5. **Lagrange Basis**: Interpolation polynomials
   - Each basis polynomial is 1 at one node, 0 at others
   - Used in finite element methods

**Note**: Most polynomial system solvers work exclusively in power basis. Conversion between bases is possible but:
- Can introduce numerical errors
- May lose sparsity structure
- Typically done as preprocessing/postprocessing step

For 2-variable polynomial systems, basis choice mainly affects:
- Numerical stability of coefficient representation
- Efficiency of evaluation
- Conditioning of linear systems in Gröbner basis computation

### Multivariate Bernstein Basis for 2-Variable Systems:

The Bernstein basis extends to bivariate polynomials over rectangular domains [a,b] × [c,d]:

**Tensor Product Bernstein Basis**:
- Basis functions: Bᵢ,ⱼ(x,y) = bᵢ,ₘ(x) · bⱼ,ₙ(y)
- Polynomial: p(x,y) = Σᵢ₌₀ᵐ Σⱼ₌₀ⁿ cᵢⱼ Bᵢ,ⱼ(x,y)
- Control net: {cᵢⱼ} forms a grid of control points

**Properties for 2-Variable Root Finding**:
1. **Convex Hull Property**: Surface p(x,y) lies within convex hull of control net
2. **Subdivision**: Tensor product de Casteljau algorithm for domain subdivision
3. **Root Exclusion**: If all cᵢⱼ > 0 or all cᵢⱼ < 0, no roots in domain
4. **Bounding Box**: Can compute tight axis-aligned bounding box from control points

**Application to 2-Variable Systems**:
- Convert f(x,y) and g(x,y) to Bernstein form on rectangular domain
- Use convex hull property to exclude regions without roots
- Subdivide domain recursively (quadtree subdivision)
- Combine with interval arithmetic or linear programming for tighter bounds
- Particularly effective for polynomial systems arising from geometric problems

**Challenges**:
- Conversion from power basis to Bernstein basis can be expensive for high degrees
- Number of coefficients grows as (m+1)(n+1) for degrees m,n
- Subdivision creates many sub-problems (exponential growth)
- Most effective when combined with other techniques (resultants, Gröbner bases)

## SYMBOLIC METHODS

### 1. Resultant Method (Sylvester/Bezout) **[WIDELY USED]**
- **Type**: Symbolic/Exact
- **Description**: Eliminates one variable by computing determinant of Sylvester matrix
- **Core Algorithm**:
  1. Given f(x,y) and g(x,y), treat as polynomials in y with coefficients in ℚ[x]
  2. Construct Sylvester matrix S (size (deg_y(f) + deg_y(g)) × (deg_y(f) + deg_y(g)))
  3. Compute resultant R(x) = det(S), a univariate polynomial in x
  4. Find roots of R(x) to get x-coordinates
  5. Substitute back to find corresponding y-coordinates
- **Complexity**: O(d³) for degree d polynomials
- **Applications**:
  - Computer algebra systems (SymPy, SageMath, Mathematica, Maple)
  - Computational geometry (curve/surface intersection)
  - Computer graphics (ray-surface intersection, spoly for specular paths)
  - Robotics (kinematic equations)
  - Algebraic geometry (elimination theory)
- **Variants**:
  - **Sylvester Resultant**: Standard matrix-based method
  - **Bezout Resultant**: Alternative formulation
  - **Hidden Variable Resultant**: For bivariate systems (used in spoly)
- **Key References**:
  - Sylvester, J.J. (1840): "A Method of Determining by Mere Inspection the Derivatives from Two Equations of Any Degree"
  - Cox, Little, O'Shea (2007): "Ideals, Varieties, and Algorithms", Chapter 3
  - Fan, Z., et al. (2024): "Specular Polynomials" SIGGRAPH 2024 (hidden variable resultant)
  - Wikipedia: https://en.wikipedia.org/wiki/Resultant
- **C/Python Implementations**:
  - **SymPy** (Python): Built-in resultant computation
  - **SageMath** (Python): Polynomial resultant methods
  - **Macaulay2** (C++ core): Resultant computations
  - **spoly** (C++): Hidden variable resultant for bivariate systems

### 2. Gröbner Basis Methods **[WIDELY USED]**
#### 2a. Buchberger's Algorithm
- **Type**: Symbolic/Exact
- **Description**: Original algorithm for computing Gröbner bases
- **Core Algorithm**:
  1. Start with polynomial set F = {f₁, f₂, ...}
  2. Compute S-polynomials: S(fᵢ, fⱼ) = (lcm/LT(fᵢ))·fᵢ - (lcm/LT(fⱼ))·fⱼ
  3. Reduce S-polynomial by current basis using polynomial division
  4. If remainder ≠ 0, add to basis
  5. Repeat until all S-polynomials reduce to 0
  6. Result is Gröbner basis G
  7. For zero-dimensional systems, convert to lexicographic ordering to get triangular form
- **Complexity**: Doubly exponential in worst case
- **Applications**:
  - Computer algebra systems (Mathematica, Maple, SageMath, Macaulay2, SINGULAR)
  - Algebraic geometry (ideal theory)
  - Robotics (kinematic constraints)
  - Cryptography (algebraic attacks)
  - Theorem proving
- **Key References**:
  - Buchberger, B. (1965): "Ein Algorithmus zum Auffinden der Basiselemente des Restklassenringes nach einem nulldimensionalen Polynomideal" (PhD Thesis)
  - Buchberger, B. (2006): "Bruno Buchberger's PhD thesis 1965: An algorithm for finding the basis elements of the residue class ring of a zero dimensional polynomial ideal", Journal of Symbolic Computation, 41(3-4), 475-511
  - Wikipedia: https://en.wikipedia.org/wiki/Gröbner_basis
- **C/Python Implementations**:
  - **Macaulay2** (C++ core, GPL): https://github.com/Macaulay2/M2
  - **SINGULAR** (C/C++): Computer algebra system
  - **SymPy** (Python): groebner() function
  - **SageMath** (Python): Gröbner basis computations

#### 2b. F4 Algorithm **[WIDELY USED]**
- **Type**: Symbolic/Exact (matrix-based)
- **Description**: Replaces many S-polynomial reductions with single large matrix reduction
- **Core Algorithm**:
  1. Select degree d, compute all S-polynomials of degree ≤ d
  2. Build Macaulay matrix M containing all polynomials and their multiples
  3. Perform Gaussian elimination on M (row echelon form)
  4. Extract new basis elements from non-zero rows
  5. Increment d and repeat until no new polynomials
  6. Key optimization: symbolic preprocessing to identify pivot columns
- **Complexity**: More efficient than Buchberger in practice
- **Applications**:
  - **msolve** (fastest open-source implementation, used in Oscar.jl, SageMath)
  - Cryptography (algebraic attacks on ciphers)
  - Computational biology (biochemical networks)
  - Engineering (polynomial system solving)
  - Parallel/GPU implementations available
- **Key References**:
  - Faugère, J.C. (1999): "A new efficient algorithm for computing Gröbner bases (F4)", Journal of Pure and Applied Algebra, 139(1-3), 61-88
  - Faugère, J.C. (2002): "A new efficient algorithm for computing Gröbner bases without reduction to zero (F5)", ACM SIGSAM International Symposium on Symbolic and Algebraic Computation
  - Berthomieu, J., Eder, C., Safey El Din, M. (2021): "msolve: A Library for Solving Polynomial Systems", ISSAC 2021
- **C/Python Implementations**:
  - **msolve** (C, GPL-2.0): https://github.com/algebraic-solving/msolve
  - **FGb** (C library): Commercial, interface available for SageMath
  - **Macaulay2** (C++): Includes F4 implementation
  - **Parallel F4** (C++): Multithreaded implementations on GitHub

#### 2c. F5 Algorithm
- **Type**: Symbolic/Exact
- **Description**: Improves F4 with signature-based criterion to avoid zero reductions
- **Complexity**: Nearly optimal for regular sequences
- **Key References**:
  - Faugère, J.C. (2002): "A new efficient algorithm for computing Gröbner bases without reduction to zero (F5)", ISSAC 2002
- **C/Python Implementations**:
  - **msolve** (C): Includes F5-like optimizations
  - **Macaulay2** (C++): F5 variants

#### 2d. FGLM Algorithm
- **Type**: Symbolic/Exact (change of ordering)
- **Description**: Converts Gröbner basis from one monomial ordering to another
- **Core Algorithm**:
  1. Input: Gröbner basis G for ordering O₁ (e.g., DRL - degree reverse lexicographic)
  2. Build multiplication matrices for variables using normal form computations
  3. Compute normal forms of monomials in target ordering O₂ (e.g., LEX - lexicographic)
  4. Use linear algebra to find linear dependencies
  5. Extract Gröbner basis for ordering O₂
  6. For 2 variables: LEX ordering gives triangular system {g₁(y), g₂(x,y)}
- **Complexity**: Polynomial in number of solutions for zero-dimensional systems
- **Key References**:
  - Faugère, J.C., Gianni, P., Lazard, D., Mora, T. (1993): "Efficient Computation of Zero-Dimensional Gröbner Bases by Change of Ordering", Journal of Symbolic Computation, 16(4), 329-344
- **C/Python Implementations**:
  - **msolve** (C): Includes FGLM implementation
  - **Macaulay2** (C++): FGLM algorithm
  - **SageMath** (Python): change_ring() with FGLM

### 3. Regular Chains / Triangular Decomposition
- **Type**: Symbolic/Exact
- **Description**: Decomposes polynomial system into triangular form
- **Complexity**: Varies by method
- **Key References**:
  - Aubry, P., Maza, M. Moreno (1999): "Triangular Sets for Solving Polynomial Systems: a Comparative Implementation of Four Methods", Journal of Symbolic Computation, 28(1-2), 125-154
  - Lazard, D. (1992): "Solving zero-dimensional algebraic systems", Journal of Symbolic Computation, 13(2), 117-131
- **C/Python Implementations**:
  - **RegularChains** (Maple library, C core)
  - **Epsilon** (C library)
  - **SageMath** (Python): triangular_decomposition()

### 4. Rational Univariate Representation (RUR)
- **Type**: Symbolic/Exact
- **Description**: Represents solutions using separating variable and rational functions
- **Core Algorithm**:
  1. Find linear combination u = a·x + b·y that separates all solutions
  2. Compute minimal polynomial f(u) of degree D (number of solutions)
  3. Express x and y as rational functions: x = g(u)/f'(u), y = h(u)/f'(u)
  4. Solutions obtained by: find roots uᵢ of f(u), then (xᵢ, yᵢ) = (g(uᵢ)/f'(uᵢ), h(uᵢ)/f'(uᵢ))
  5. Compact representation: store only 3 univariate polynomials instead of D points
- **Complexity**: Efficient for zero-dimensional systems
- **Key References**:
  - Rouillier, F. (1999): "Solving Zero-Dimensional Systems Through the Rational Univariate Representation", Applied Algebra in Engineering, Communication and Computing, 9(9), 433-461
  - Rouillier, F., Zimmermann, P. (2004): "Efficient isolation of polynomial's real roots", Journal of Computational and Applied Mathematics, 162(1), 33-50
- **C/Python Implementations**:
  - **msolve** (C): Includes RUR computation
  - **SageMath** (Python): rational_univariate_representation()
  - **Macaulay2** (C++): RUR methods

### 5. Wu's Method (Characteristic Set)
- **Type**: Symbolic/Exact
- **Description**: Chinese method for solving polynomial systems via characteristic sets
- **Complexity**: Varies
- **Key References**:
  - Wu, W.T. (1984): "Basic principles of mechanical theorem proving in elementary geometries", Journal of Automated Reasoning, 2, 221-252
  - Wikipedia: https://en.wikipedia.org/wiki/Wu%27s_method_of_characteristic_set
- **C/Python Implementations**:
  - **Epsilon** (C library)
  - **SageMath** (Python): Some Wu method implementations
  - **Maple**: CharSets package

## NUMERICAL METHODS

### 6. Homotopy Continuation **[WIDELY USED]**
- **Type**: Numerical/Approximate
- **Description**: Tracks solution paths from start system to target system
- **Core Algorithm**:
  1. Define homotopy H(x, t) = (1-t)·G(x) + t·F(x) where G is start system, F is target
  2. Start with known solutions of G(x) = 0 at t=0
  3. Track solution paths as t varies from 0 to 1 using predictor-corrector:
     - Predictor: Euler/Runge-Kutta step along path
     - Corrector: Newton iteration to stay on path
  4. Solutions at t=1 are solutions to F(x) = 0
  5. Handle path crossing, singularities, paths at infinity
- **Complexity**: Polynomial in number of paths (Bézout bound)
- **Applications**:
  - **Engineering**: Kinematics, mechanism design, robotics
  - **Chemistry**: Chemical equilibrium, reaction networks
  - **Physics**: Scattering amplitudes, particle physics
  - **Computer vision**: 3D reconstruction, camera calibration
  - **Numerical algebraic geometry**: Finding all isolated solutions
  - **Optimization**: Global optimization via polynomial systems
- **Advantages**:
  - Finds **all** isolated complex solutions (not just real)
  - Robust (doesn't require good initial guess)
  - Parallel-friendly (independent path tracking)
  - Handles positive-dimensional solution sets
- **Key References**:
  - Morgan, A. (1987): "Solving Polynomial Systems Using Continuation for Engineering and Scientific Problems", SIAM
  - Sommese, A.J., Wampler, C.W. (2005): "The Numerical Solution of Systems of Polynomials Arising in Engineering and Science", World Scientific
  - Verschelde, J. (1999): "Algorithm 795: PHCpack: A general-purpose solver for polynomial systems by homotopy continuation", ACM Transactions on Mathematical Software, 25(2), 251-276
  - Bates, D.J., Hauenstein, J.D., Sommese, A.J., Wampler, C.W. (2013): "Numerically Solving Polynomial Systems with Bertini", SIAM
  - Breiding, P., Timme, S. (2018): "HomotopyContinuation.jl: A Package for Homotopy Continuation in Julia", International Congress on Mathematical Software, 458-465
- **C/Python Implementations**:
  - **PHCpack** (C, Ada): https://github.com/janverschelde/PHCpack
  - **Bertini** (C): https://bertini.nd.edu/
  - **HomotopyContinuation.jl** (Julia, but has C backend): https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl
  - **phcpy** (Python interface to PHCpack): Available via pip

### 7. Newton-Raphson Method
- **Type**: Numerical/Approximate (local)
- **Description**: Iterative method for finding roots, requires good initial guess
- **Core Algorithm**:
  1. Start with initial guess x⁽⁰⁾
  2. Compute Jacobian matrix J(x) = [∂fᵢ/∂xⱼ]
  3. Iterate: x⁽ᵏ⁺¹⁾ = x⁽ᵏ⁾ - J(x⁽ᵏ⁾)⁻¹ · F(x⁽ᵏ⁾)
  4. Stop when ||F(x⁽ᵏ⁾)|| < tolerance or ||x⁽ᵏ⁺¹⁾ - x⁽ᵏ⁾|| < tolerance
  5. Only finds one root per initial guess
  6. Requires good initial guess for convergence
- **Complexity**: Quadratic convergence near solution
- **Key References**:
  - Standard numerical analysis textbooks
  - Press, W.H. et al. (2007): "Numerical Recipes: The Art of Scientific Computing", 3rd Edition
- **C/Python Implementations**:
  - **SciPy** (Python): scipy.optimize.fsolve, scipy.optimize.root
  - **NumPy** (Python/C): Basic implementations
  - **GSL** (C): GNU Scientific Library multiroot solvers
  - **NLopt** (C): Nonlinear optimization library

### 8. Interval Arithmetic / Subdivision Methods
- **Type**: Numerical with guarantees
- **Description**: Uses interval arithmetic to guarantee solution enclosures
- **Core Algorithm**:
  1. Start with initial box B = [x₁, x₂] × [y₁, y₂]
  2. Evaluate f(B) and g(B) using interval arithmetic
  3. If 0 ∉ f(B) or 0 ∉ g(B), discard box (no solution)
  4. If box is small enough and contains solution, store it
  5. Otherwise, subdivide box into smaller boxes
  6. Apply contractor operators (interval Newton, box consistency) to shrink boxes
  7. Recursively process sub-boxes
  8. Result: list of boxes guaranteed to contain all solutions

#### 8a. Bernstein Basis Subdivision Method
- **Type**: Subdivision with convex hull property
- **Description**: Exploits convex hull property of Bernstein basis for root isolation
- **Core Algorithm**:
  1. Convert polynomial p(x) from power basis to Bernstein basis on interval [a,b]
  2. Bernstein coefficients b₀, b₁, ..., bₙ form control points
  3. **Convex Hull Property**: Graph of p(x) lies within convex hull of control points
  4. **Root Exclusion Test**:
     - If all bᵢ > 0 or all bᵢ < 0 (no sign change), no roots in [a,b]
     - Uses variation-diminishing property (Descartes' rule for Bernstein basis)
  5. **Subdivision**: If sign changes exist, subdivide at midpoint using de Casteljau algorithm
  6. Recursively process sub-intervals until roots are isolated
  7. **Advantages**: Better numerical stability than power basis, natural geometric interpretation

- **Convex Hull Computation Algorithms** (for control points):
  - **Graham Scan**: O(n log n) - Sort points by angle, scan to build hull
  - **Jarvis March (Gift Wrapping)**: O(nh) where h = hull vertices - Good for small hulls
  - **Andrew's Monotone Chain**: O(n log n) - Sort by x-coordinate, build upper/lower chains
  - **QuickHull**: O(n log n) average, O(n²) worst - Divide and conquer approach
  - **Chan's Algorithm**: O(n log h) - Combines Graham scan and Jarvis march

  For Bernstein polynomial root finding, typically use:
  - **Simple bounding box** (axis-aligned) for quick exclusion tests
  - **Linear programming** to find intersection of convex hull with x-axis
  - **Sign variation counting** (most common) - just check sign changes in coefficients

#### 8b. ProjectedPolyhedral Method (This Project)
- **Type**: Subdivision with direction-by-direction projection
- **Description**: Computes root bounding box by projecting graph control points direction-by-direction
- **Core Algorithm** (for each direction i):
  1. For each equation, project all graph control points to 2D (coordinate i + function value)
  2. For each equation, compute convex hull of these 2D points
  3. For each equation, intersect convex hull with horizontal axis (function value = 0)
  4. For each equation, project intersection to 1D to get an interval
  5. Intersect all intervals from all equations to get the bound in direction i
- **Advantages**:
  - Simpler than full graph-space convex hull (GraphHull method)
  - More modular: each direction processed independently
  - Works for multivariate systems (1D and 2D implemented)
  - Numerically stable with Bernstein basis
- **Implementation**: Available in this polynomial-solver project (C++11)
- **References**: See `src/solver.cpp` and `tests/test_projected_polyhedral.cpp` in this repository

#### 8c. GraphHull Method (This Project)
- **Type**: Subdivision with exact convex hull in graph space
- **Description**: Computes exact convex hull of graph control points in R^{n+1}, intersects with hyperplane x_{n+1} = 0
- **Core Algorithm**:
  1. For each equation, compute convex hull of graph control points in R^{n+1}
  2. For each equation, intersect hull with hyperplane x_{n+1} = 0 (where function value = 0)
  3. For each equation, project to R^n parameter space
  4. Intersect all projected polyhedra from all equations
- **Advantages**:
  - Exact for linear functions (machine epsilon error: 2.22×10⁻¹⁶)
  - Tighter bounds than ProjectedPolyhedral in some cases
  - Rigorous geometric foundation
- **Limitations**: Only implemented for 1D and 2D systems (graphs in R² and R³)
- **Implementation**: Available in this polynomial-solver project (C++11)
- **References**: See `src/solver.cpp` and `tests/test_linear_graph_hull.cpp` in this repository

- **Complexity**: Varies
- **Key References**:
  - Moore, R.E. (1966): "Interval Analysis", Prentice-Hall
  - Neumaier, A. (1990): "Interval Methods for Systems of Equations", Cambridge University Press
  - **Bernstein Basis Specific**:
    - Spencer, M.R. (1994): "Polynomial Real Root Finding in Bernstein Form", PhD Thesis, Brigham Young University
    - Mourrain, B., Pavone, J.P. (2009): "Subdivision methods for solving polynomial equations", Journal of Symbolic Computation
    - Lane, J.M., Riesenfeld, R.F. (1980): "A Theoretical Development for the Computer Generation and Display of Piecewise Polynomial Surfaces"
    - de Casteljau, P. (1959): "Outillages méthodes calcul" (de Casteljau algorithm)
    - Farouki, R.T. (2012): "The Bernstein polynomial basis: A centennial retrospective", Computer Aided Geometric Design
  - **Convex Hull Algorithms**:
    - Graham, R.L. (1972): "An Efficient Algorithm for Determining the Convex Hull of a Finite Planar Set"
    - Jarvis, R.A. (1973): "On the Identification of the Convex Hull of a Finite Set of Points in the Plane"
    - Andrew, A.M. (1979): "Another Efficient Algorithm for Convex Hulls in Two Dimensions"
    - Barber, C.B., Dobkin, D.P., Huhdanpaa, H. (1996): "The Quickhull Algorithm for Convex Hulls"
    - Chan, T.M. (1996): "Optimal output-sensitive convex hull algorithms in two and three dimensions"
- **C/Python Implementations**:
  - **MPFI** (C): Multiple Precision Floating-point Interval library
  - **pyinterval** (Python): Interval arithmetic in Python
  - **mpmath** (Python): Arbitrary-precision interval arithmetic
  - **Arb** (C): Ball arithmetic library
  - **Bernstein Subdivision**:
    - **SciPy** (Python): scipy.spatial.ConvexHull for convex hull computation
    - **CGAL** (C++): Computational Geometry Algorithms Library with convex hull
    - **Qhull** (C): General dimension convex hull (used by SciPy)
    - Custom implementations in computer graphics libraries (Bézier curve subdivision)

### 9. Cylindrical Algebraic Decomposition (CAD)
- **Type**: Symbolic-Numeric hybrid
- **Description**: Decomposes space into cells where polynomials have constant sign
- **Core Algorithm** (Collins' algorithm for 2 variables):
  1. **Projection phase**: Compute projection polynomials from f, g
     - Resultants, discriminants, leading coefficients
     - Project to 1D to get critical x-values
  2. **Base phase**: Find all roots of projection polynomials (sample points on x-axis)
  3. **Lifting phase**: For each x-interval, solve for y
     - Evaluate polynomials at sample x-values
     - Find y-roots and construct cells
  4. Result: Decomposition of ℝ² into cells where f, g have constant sign
  5. Solutions are in cells where both f=0 and g=0
- **Complexity**: Doubly exponential
- **Key References**:
  - Collins, G.E. (1975): "Quantifier elimination for real closed fields by cylindrical algebraic decomposition", Automata Theory and Formal Languages, Lecture Notes in Computer Science, 33, 134-183
  - Basu, S., Pollack, R., Roy, M.F. (2006): "Algorithms in Real Algebraic Geometry", Springer
- **C/Python Implementations**:
  - **QEPCAD** (C): Quantifier Elimination by Partial CAD
  - **Mathematica**: Built-in CAD functions
  - **Maple**: CAD implementations
  - **SageMath** (Python): Interfaces to CAD tools

## HYBRID METHODS

### 10. Numerical Algebraic Geometry
- **Type**: Hybrid (symbolic preprocessing + numerical solving)
- **Description**: Combines symbolic and numerical techniques
- **Key References**:
  - Sommese, A.J., Wampler, C.W. (2005): "The Numerical Solution of Systems of Polynomials"
- **C/Python Implementations**:
  - **Bertini** (C): Numerical algebraic geometry software
  - **PHCpack** (C, Ada): Includes numerical algebraic geometry methods
  - **HomotopyContinuation.jl** (Julia): Modern NAG implementation
  - **Macaulay2** (C++): NAG package

