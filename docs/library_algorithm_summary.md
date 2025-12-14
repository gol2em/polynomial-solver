# Polynomial Solving Libraries and Algorithms Summary

This document provides a comprehensive mapping of polynomial solving algorithms to their implementations across various libraries and tools.

## Table of Contents
1. [Most Widely Used Algorithms](#most-widely-used-algorithms)
2. [Univariate Polynomial Solvers](#univariate-polynomial-solvers)
3. [Multivariate/Bivariate Polynomial System Solvers](#multivariatebivariate-polynomial-system-solvers)
4. [Specialized Libraries](#specialized-libraries)
5. [Computer Algebra Systems](#computer-algebra-systems)
6. [Degeneracy Handling](#degeneracy-handling)

---

## MOST WIDELY USED ALGORITHMS

This section highlights algorithms marked as **[WIDELY USED]** in the detailed documentation, indicating broad adoption across major libraries and practical applications.

### Univariate Polynomial Solving

1. **Quadratic Formula** (Degree 2)
   - **Applications**: Physics, engineering, graphics (ray-sphere intersection), finance
   - **Libraries**: All major systems (NumPy, MATLAB, Mathematica, etc.)
   - **Why widely used**: Simple, exact, O(1) complexity

2. **Companion Matrix Method**
   - **Applications**: NumPy/SciPy default, MATLAB roots(), control theory, signal processing
   - **Libraries**: NumPy, LAPACK, Eigen, MATLAB
   - **Why widely used**: Finds all roots (real + complex), leverages mature LAPACK eigenvalue solvers, industry standard

3. **Newton-Raphson Method**
   - **Applications**: Most common numerical root-finding, optimization, scientific computing, machine learning
   - **Libraries**: SciPy, GSL, MATLAB, NumPy
   - **Why widely used**: Quadratic convergence, simple implementation, generalizes to multivariate systems

4. **Aberth Method**
   - **Applications**: MPSolve (arbitrary precision), high-performance polynomial solving, parallel/GPU computing
   - **Libraries**: MPSolve, mpmath
   - **Why widely used**: Finds all roots simultaneously, cubic convergence, parallel-friendly

5. **Sturm's Theorem**
   - **Applications**: Root isolation, computer algebra systems, symbolic computation, theorem proving
   - **Libraries**: SymPy, SageMath, Macaulay2, Mathematica
   - **Why widely used**: Exact root counting in intervals, fundamental for symbolic methods

### Multivariate/Bivariate Polynomial System Solving

1. **Resultant Method (Sylvester/Bezout)**
   - **Applications**: Computer algebra systems, computational geometry, computer graphics (spoly), robotics, algebraic geometry
   - **Libraries**: SymPy, SageMath, Mathematica, Maple, Macaulay2, spoly
   - **Why widely used**: Direct elimination method, converts 2D to 1D problem, well-understood theory

2. **Gröbner Basis Methods (Buchberger, F4, F5)**
   - **Applications**: Computer algebra systems, algebraic geometry, robotics, cryptography, theorem proving
   - **Libraries**: msolve (F4), Macaulay2, SINGULAR, SymPy, SageMath, Mathematica, Maple
   - **Why widely used**: General-purpose method for polynomial ideals, handles arbitrary systems, foundation of computational algebra

3. **Homotopy Continuation**
   - **Applications**: Engineering (kinematics, robotics), chemistry, physics, computer vision, numerical algebraic geometry
   - **Libraries**: PHCpack, Bertini, HomotopyContinuation.jl
   - **Why widely used**: Finds **all** isolated complex solutions, robust (no initial guess needed), parallel-friendly, handles positive-dimensional sets

### Why These Algorithms Dominate

**Companion Matrix Method** (univariate):
- Leverages decades of LAPACK development
- Single method handles all polynomial degrees
- Numerically stable for moderate-degree polynomials
- Default in NumPy, MATLAB

**Gröbner Bases** (multivariate):
- Most general symbolic method
- Handles arbitrary polynomial systems
- Foundation of computational algebra
- Implemented in all major CAS

**Homotopy Continuation** (multivariate numerical):
- Only numerical method guaranteed to find all isolated solutions
- Doesn't require initial guess (unlike Newton-Raphson)
- Mature implementations (PHCpack since 1999, Bertini since 2006)
- Parallel-friendly architecture

---

## UNIVARIATE POLYNOMIAL SOLVERS

### Major Libraries

#### **msolve** (C)
- **URL**: https://msolve.lip6.fr / https://github.com/algebraic-solving/msolve
- **Language**: C
- **License**: GPL-2.0
- **Algorithms**:
  - Gröbner Basis (F4 algorithm variant)
  - Real root isolation
  - Multi-threading support
- **Features**:
  - Rational coefficients and prime fields
  - Dimension and degree computation
  - Multi-threaded linear algebra
  - Multi-modular computations
- **Integration**: Used in Oscar.jl, SageMath, AlgebraicSolving.jl
- **Status**: Active (2,434 commits, last updated 2025)

#### **MPSolve** (C)
- **URL**: https://numpi.dm.unipi.it/software/mpsolve
- **Language**: C
- **Algorithms**:
  - Aberth Method (simultaneous approximation)
  - Secular equation solver
- **Features**:
  - Arbitrary precision arithmetic
  - Parallel computation
  - All roots simultaneously
- **Status**: Mature, widely used

#### **NumPy/SciPy** (Python)
- **URL**: https://numpy.org, https://scipy.org
- **Language**: Python (C backend)
- **Algorithms**:
  - Companion Matrix Method (via numpy.linalg.eigvals)
  - Laguerre's Method (numpy.roots uses companion matrix)
- **Features**:
  - Fast for moderate-degree polynomials
  - Leverages LAPACK
- **Status**: Industry standard

#### **GSL (GNU Scientific Library)** (C)
- **URL**: https://www.gnu.org/software/gsl/
- **Language**: C
- **Algorithms**:
  - Companion Matrix Method
  - Balanced QR reduction
- **Features**:
  - Robust, well-tested
  - Part of GNU ecosystem
- **Status**: Mature, stable

#### **MPFI / Arb** (C)
- **URL**: http://perso.ens-lyon.fr/nathalie.revol/software.html (MPFI)
- **URL**: https://arblib.org/ (Arb)
- **Language**: C
- **Algorithms**:
  - Interval Newton Method
  - Interval arithmetic
- **Features**:
  - Guaranteed bounds
  - Arbitrary precision
- **Status**: Active

### Specialized GitHub Projects

#### **ZhepeiWang/Root-Finder** (C++)
- **URL**: https://github.com/ZhepeiWang/Root-Finder
- **Language**: C++
- **Year**: 2019
- **Algorithms**:
  - **Real Roots Isolation Method** (geometric bounds-based)
  - Cauchy's Bound
  - **Kojima's Bound** (tighter than Cauchy)
  - **Fujiwara's Bound** (theoretically best)
- **Performance**: 5% faster than Jenkins-Traub for 8th-order polynomials
- **Features**:
  - Can find roots in specific intervals
  - More stable than Companion Matrix
  - Recommended for up to 32-order polynomials
- **Status**: Research project

#### **aureooms-research/roots** (Python)
- **URL**: https://github.com/aureooms-research/roots
- **Language**: Python
- **Algorithms**:
  - **Thom Encoding** (symbolic root representation)
  - Sturm sequences
- **Features**:
  - Symbolic root comparison without numerical computation
  - Exact symbolic manipulation
- **Status**: Research project

#### **krlu/Abel** (Scala)
- **URL**: https://github.com/krlu/Abel
- **Language**: Scala
- **Algorithms**:
  - **Kronecker Factorization**
  - Integer factorization-based methods
- **Features**:
  - Exact factorization over integers
  - Computer algebra applications
- **Status**: Educational/research project

#### **ClecioJung/studying-c** (C)
- **URL**: https://github.com/ClecioJung/studying-c
- **Language**: C
- **Algorithms**:
  - Newton-Raphson
  - Bisection
  - Various numerical methods
- **Features**: Educational implementations
- **Status**: Learning resource

---

## MULTIVARIATE/BIVARIATE POLYNOMIAL SYSTEM SOLVERS

### Major Computer Algebra Systems

#### **msolve** (C) - See above
- **Multivariate Features**:
  - Gröbner bases for 0-dimensional ideals
  - Real solution isolation (rational coefficients)
  - Parametrization (prime fields)
  - One-block elimination orders
- **Degeneracy Handling**: Dimension detection, degree computation

#### **Macaulay2**
- **URL**: https://macaulay2.com/
- **Language**: C++
- **Algorithms**:
  - F4 Algorithm
  - Gröbner Basis (Buchberger, F4)
  - Regular Chains
- **Features**: Research-grade CAS
- **Status**: Active

#### **SINGULAR**
- **URL**: https://www.singular.uni-kl.de/
- **Language**: C/C++
- **Algorithms**:
  - Gröbner Basis (optimized implementations)
  - Standard bases
  - Syzygy computation
- **Features**: Specialized for commutative algebra
- **Status**: Active

#### **CoCoA**
- **URL**: http://cocoa.dima.unige.it/
- **Language**: C++
- **Algorithms**:
  - Gröbner Basis
  - Buchberger's Algorithm
- **Features**: Educational and research
- **Status**: Active

### Numerical Solvers for Multivariate Systems

#### **PHCpack**
- **URL**: http://homepages.math.uic.edu/~jan/PHCpack/phcpack.html
- **Language**: Ada/C
- **Algorithms**:
  - **Homotopy Continuation** (projected polyhedral, total degree)
  - Predictor-corrector methods
- **Features**:
  - Finds all isolated complex solutions
  - Parallel computation
  - Mixed volume computation
- **Status**: Mature, widely used

#### **Bertini**
- **URL**: https://bertini.nd.edu/
- **Language**: C/C++
- **Algorithms**:
  - **Homotopy Continuation**
  - Adaptive precision
- **Features**:
  - Numerical algebraic geometry
  - Witness sets
  - Positive-dimensional solution sets
- **Status**: Active

#### **HomotopyContinuation.jl** (Julia)
- **URL**: https://www.juliahomotopycontinuation.org/
- **Language**: Julia
- **Algorithms**:
  - Homotopy Continuation
  - Monodromy methods
- **Features**:
  - Modern, fast implementation
  - Parameter homotopies
- **Status**: Active, modern

### Specialized GitHub Projects for Bivariate Systems

#### **mollnn/spoly** (C++)
- **URL**: https://github.com/mollnn/spoly
- **Language**: C++
- **Year**: 2024 (SIGGRAPH 2024)
- **Paper**: "Specular Polynomials"
- **Algorithms**:
  - **Resultant Method** (bivariate to univariate reduction)
  - **Hidden Variable Resultant Method**
  - **Univariate Matrix Polynomial Determinant**
  - Laplacian expansion
  - Bisection solver
- **Application**: Specular path finding in rendering (glints, caustics)
- **Features**:
  - Converts bivariate polynomial systems to univariate matrix polynomials
  - Finds zeros of determinant
  - GPU-friendly (CUDA implementation available)
  - Completely deterministic (no Newton iteration)
- **Degeneracy Handling**: Handles specular constraints in rendering
- **Status**: Research project, 85 stars

#### **louislegrain/bivariate-polynomial-equation-solver** (Python)
- **URL**: https://github.com/louislegrain/bivariate-polynomial-equation-solver
- **Language**: Python
- **Year**: 2024
- **Algorithms**: Method for finding solutions to bivariate polynomial equations
- **Status**: Small research project (2 stars)

### Parallel/GPU Implementations

#### **Jbowman353/parallel-groebner** (Python)
- **URL**: https://github.com/Jbowman353/parallel-groebner
- **Language**: Python (GPU backend)
- **Algorithms**:
  - **Parallel Gröbner Basis on GPU**
  - Buchberger's Algorithm (parallelized)
- **Features**: GPU acceleration for Gröbner basis computation
- **Status**: Research project

#### **JohnS0819/Faugere-s-F4-algorithm-over-finite-fields-with-multithreaded-matrix-reduction** (C++)
- **URL**: https://github.com/JohnS0819/Faugere-s-F4-algorithm-over-finite-fields-with-multithreaded-matrix-reduction
- **Language**: C++
- **Algorithms**:
  - **F4 Algorithm with Multithreading**
  - Parallel matrix reduction
- **Features**: Optimized for finite fields
- **Status**: Research project

#### **wery0/GroebnerBasisLibrary** (C++)
- **URL**: https://github.com/wery0/GroebnerBasisLibrary
- **Language**: C++
- **Algorithms**: Gröbner Basis computation
- **Status**: Library project

#### **SunQpark/Multivariate_Polynomials** (C++)
- **URL**: https://github.com/SunQpark/Multivariate_Polynomials
- **Language**: C++
- **Algorithms**: Buchberger's Algorithm
- **Status**: Educational project

### JavaScript/Web Implementations

#### **antimatter15/groebner.js** (JavaScript)
- **URL**: https://github.com/antimatter15/groebner.js
- **Language**: JavaScript
- **Algorithms**: Buchberger's Algorithm
- **Features**: Browser-based Gröbner basis computation
- **Status**: Educational/demo project

---

## SPECIALIZED LIBRARIES

### Polynomial Pseudo-Division

#### **zhaojin1997/Polynomial_pseudo_division** (C++)
- **URL**: https://github.com/zhaojin1997/Polynomial_pseudo_division
- **Language**: C++
- **Algorithms**: Polynomial pseudo-division for large integer coefficients
- **Features**: Handles large coefficients efficiently
- **Status**: Research project

### Subdivision-Based Methods (Bernstein Basis)

#### **polynomial-solver** (This Project - C++)
- **URL**: https://github.com/gol2em/polynomial-solver
- **Language**: C++11
- **Algorithms**:
  - **ProjectedPolyhedral (PP)**: Direction-by-direction projection of graph control points
  - **GraphHull**: Exact convex hull in graph space R^{n+1}
  - Bernstein basis subdivision
  - De Casteljau algorithm for subdivision
- **Features**:
  - Multivariate polynomial support (1D and 2D systems)
  - High-precision result refinement (Newton's method)
  - Condition number estimation
  - Multiplicity detection
  - Degeneracy detection
  - Machine epsilon precision for linear systems (2.22×10⁻¹⁶)
  - Direct contraction to minimize error accumulation
- **Root Bounding Methods**:
  - **ProjectedPolyhedral**: Simpler, modular, direction-by-direction
  - **GraphHull**: Exact for linear functions, tighter bounds
  - **None**: Uniform subdivision (baseline)
- **Status**: Active development
- **License**: TBD

### Multivariate Orthogonal Polynomials

#### **JuliaApproximation/MultivariateOrthogonalPolynomials.jl** (Julia)
- **URL**: https://github.com/JuliaApproximation/MultivariateOrthogonalPolynomials.jl
- **Language**: Julia
- **Algorithms**: Multivariate orthogonal polynomial methods
- **Features**: Approximation theory applications
- **Status**: Active

### Recursive Polynomial Adaptive Sampling

#### **Hudson-257/Recursive-Polynomial-Adaptive-Sampling** (Python)
- **URL**: https://github.com/Hudson-257/Recursive-Polynomial-Adaptive-Sampling
- **Language**: Python
- **Algorithms**: Adaptive sampling for polynomial systems
- **Status**: Research project

---

## COMPUTER ALGEBRA SYSTEMS

### **SymPy** (Python)
- **URL**: https://www.sympy.org/
- **Language**: Python
- **Algorithms**:
  - Symbolic root finding (solve, solveset)
  - Gröbner Basis (Buchberger)
  - Resultants
  - Sturm sequences
- **Features**: Pure Python, symbolic computation
- **Status**: Active, widely used

### **SageMath** (Python)
- **URL**: https://www.sagemath.org/
- **Language**: Python (wraps many C/C++ libraries)
- **Algorithms**:
  - Uses msolve for real root isolation
  - Gröbner bases (via Singular)
  - Variety computation
- **Integration**: Wraps msolve, Singular, Macaulay2, etc.
- **Status**: Active, comprehensive

### **Oscar.jl** (Julia)
- **URL**: https://oscar-system.github.io/Oscar.jl
- **Language**: Julia
- **Algorithms**:
  - Uses msolve for solving polynomial systems
  - Gröbner bases
  - Real solutions via msolve
- **Integration**: Wraps msolve and other libraries
- **Status**: Active, modern

### **Maple**
- **URL**: https://www.maplesoft.com/
- **Language**: Proprietary
- **Algorithms**:
  - Comprehensive polynomial solving
  - Gröbner bases
  - Cylindrical Algebraic Decomposition
- **Status**: Commercial, mature

### **Mathematica**
- **URL**: https://www.wolfram.com/mathematica/
- **Language**: Proprietary
- **Algorithms**:
  - NSolve, Solve, Reduce
  - Gröbner bases
  - CAD
- **Status**: Commercial, mature

---

## DEGENERACY HANDLING

### Techniques Found in Libraries

1. **Dimension Detection** (msolve, Macaulay2, Singular)
   - Detects if system has dimension > 0
   - Computes Hilbert function/polynomial

2. **Multiplicity Handling**
   - GCD-based multiple root detection
   - Deflation techniques (Bertini)

3. **Singular Locus Computation**
   - Jacobian-based methods
   - Available in Macaulay2, Singular

4. **Regularization Methods**
   - Perturbation techniques
   - Used in homotopy continuation (PHCpack, Bertini)

5. **Zero-Dimensional Ideal Detection**
   - Leading monomial analysis
   - Gröbner basis properties

### Libraries with Explicit Degeneracy Support

- **msolve**: Dimension and degree computation, detects 0-dimensional ideals
- **Bertini**: Deflation for singular solutions
- **PHCpack**: Handles positive-dimensional solution sets
- **Macaulay2/Singular**: Comprehensive singular locus tools
- **mollnn/spoly**: Handles degenerate specular constraints in rendering

---

## SUMMARY TABLE: ALGORITHMS BY LIBRARY

| Algorithm | msolve | PHCpack | Bertini | Macaulay2 | Singular | SymPy | NumPy | MPSolve | Root-Finder |
|-----------|--------|---------|---------|-----------|----------|-------|-------|---------|-------------|
| Gröbner Basis (F4) | ✓ | | | ✓ | ✓ | | | | |
| Gröbner Basis (Buchberger) | | | | ✓ | ✓ | ✓ | | | |
| Homotopy Continuation | | ✓ | ✓ | | | | | | |
| Real Root Isolation | ✓ | | | | | | | | ✓ |
| Resultant Method | | | | ✓ | ✓ | ✓ | | | |
| Companion Matrix | | | | | | | ✓ | | |
| Aberth Method | | | | | | | | ✓ | |
| Kojima's Bound | | | | | | | | | ✓ |
| Thom Encoding | | | | | | | | | |
| CAD | | | | ✓ | | | | | |
| Multi-threading | ✓ | ✓ | ✓ | | | | | ✓ | |
| GPU Support | | | | | | | | | |
| Degeneracy Handling | ✓ | ✓ | ✓ | ✓ | ✓ | | | | |

---

## RECOMMENDATIONS BY USE CASE

### For Univariate Polynomials:
- **General purpose**: NumPy/SciPy (fast, easy)
- **High precision**: MPSolve, Arb
- **Interval arithmetic**: MPFI, Arb
- **Specific intervals**: ZhepeiWang/Root-Finder
- **Symbolic**: SymPy

### For Multivariate Systems (0-dimensional):
- **Rational coefficients**: msolve (fastest)
- **Numerical solutions**: PHCpack, Bertini, HomotopyContinuation.jl
- **Symbolic**: Macaulay2, Singular, SymPy
- **Integrated environment**: SageMath, Oscar.jl

### For Positive-Dimensional Systems:
- **Numerical**: Bertini, PHCpack
- **Symbolic**: Macaulay2, Singular

### For Bivariate Systems:
- **Resultant-based**: mollnn/spoly (GPU-friendly)
- **General**: msolve, Macaulay2

### For Parallel/GPU Computation:
- **GPU Gröbner**: Jbowman353/parallel-groebner
- **Multithreaded F4**: JohnS0819/Faugere-s-F4-algorithm
- **Multithreaded general**: msolve, PHCpack, Bertini

### For Degeneracy Handling:
- **Dimension detection**: msolve, Macaulay2, Singular
- **Deflation**: Bertini
- **Regularization**: PHCpack

---

## REFERENCES

1. **msolve**: Berthomieu, J., Eder, C., Safey El Din, M. (2021). "msolve: A Library for Solving Polynomial Systems." ISSAC 2021.
2. **PHCpack**: Verschelde, J. (1999). "Algorithm 795: PHCpack." ACM TOMS.
3. **Bertini**: Bates, D.J., et al. (2013). "Numerically Solving Polynomial Systems with Bertini."
4. **Specular Polynomials**: Fan, Z., et al. (2024). "Specular Polynomials." SIGGRAPH 2024.
5. **Real Roots Isolation**: Wang, Z. (2019). "Root-Finder" (GitHub project).

---

*Last updated: 2025-12-08*


