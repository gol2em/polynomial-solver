# Runtime Precision Control - Design Solutions

## The Problem

**User's concern**: "If I want to use flexible precisions, like I will implement a function to estimate how many bits are needed, how do I change it without refactor the code and minimize compile time for the user?"

### Current Private Template Design

```cpp
// include/detail/polynomial_impl.h (header-only template)
namespace detail {
    template<typename T>
    class PolynomialImpl {
        std::vector<T> bernstein_coeffs_;
        // ... algorithms
    };
}

// include/polynomial_mp.h (public API)
class PolynomialMP {
    detail::PolynomialImpl<mpreal> impl_;  // Fixed type at compile time!
};
```

**Problem**: The type `mpreal` is fixed at compile time. If you want to change precision at runtime (e.g., from 100 bits to 500 bits), you need to:
1. Recompile with different precision
2. Or create multiple instantiations (slow compilation)

## Solution 1: Runtime Precision with MPFR Backend (Recommended)

### Key Insight

With MPFR backend, `mpreal` **already supports runtime precision changes**!

```cpp
#include <boost/multiprecision/mpfr.hpp>
using mpreal = boost::multiprecision::mpfr_float;

// Create with different precisions at runtime
mpreal x(100);  // 100 bits precision
mpreal y(500);  // 500 bits precision

// Change precision dynamically
x.precision(500);  // Now 500 bits

// Set default precision for all new mpreal objects
mpreal::default_precision(256);  // All new mpreal will have 256 bits
```

### Implementation

**No code changes needed!** The template design already supports this:

```cpp
// User code - no recompilation needed
#include "polynomial_mp.h"

// Estimate precision needed
int bits_needed = estimatePrecision(problem_size);

// Set precision before creating objects
mpreal::default_precision(bits_needed);

// Create polynomial with estimated precision
PolynomialMP poly(degrees, coeffs);  // Uses bits_needed precision

// Later, increase precision if needed
mpreal::default_precision(bits_needed * 2);
PolynomialMP poly_refined(degrees, coeffs_refined);
```

**Advantages**:
- ✅ No recompilation
- ✅ No template instantiation overhead
- ✅ Change precision at any time
- ✅ Works with existing design
- ✅ Zero code refactoring

### Example: Adaptive Precision Solver

```cpp
class AdaptivePrecisionSolver {
public:
    std::vector<mpreal> solve(const PolynomialSystem& system) {
        // Start with low precision
        int precision = 100;
        mpreal::default_precision(precision);
        
        while (true) {
            // Solve with current precision
            auto result = solveLowLevel(system);
            
            // Check if precision is sufficient
            if (isPrecisionSufficient(result)) {
                return result;
            }
            
            // Increase precision and retry
            precision *= 2;
            mpreal::default_precision(precision);
            std::cout << "Increasing precision to " << precision << " bits\n";
        }
    }
    
private:
    int estimatePrecision(double condition_number) {
        // Your estimation algorithm
        return static_cast<int>(std::log2(condition_number) * 2 + 100);
    }
};
```

## Solution 2: Precision Policy Pattern (Advanced)

For more control, use a precision policy:

```cpp
// include/precision_policy.h
class PrecisionPolicy {
public:
    virtual int getPrecision(const std::string& operation) const = 0;
    virtual ~PrecisionPolicy() = default;
};

class AdaptivePrecisionPolicy : public PrecisionPolicy {
    int base_precision_;
    std::map<std::string, int> operation_precision_;
    
public:
    AdaptivePrecisionPolicy(int base) : base_precision_(base) {}
    
    int getPrecision(const std::string& operation) const override {
        auto it = operation_precision_.find(operation);
        return (it != operation_precision_.end()) ? it->second : base_precision_;
    }
    
    void setOperationPrecision(const std::string& op, int prec) {
        operation_precision_[op] = prec;
    }
};

// Usage
AdaptivePrecisionPolicy policy(100);  // Base 100 bits
policy.setOperationPrecision("newton_refinement", 500);  // 500 bits for refinement
policy.setOperationPrecision("subdivision", 200);  // 200 bits for subdivision

// In your code
mpreal::default_precision(policy.getPrecision("newton_refinement"));
auto refined = refineRoot(root);
```

## Solution 3: Precision Context Manager (Python-like)

```cpp
// include/precision_context.h
class PrecisionContext {
    int old_precision_;
    
public:
    PrecisionContext(int new_precision) 
        : old_precision_(mpreal::default_precision()) {
        mpreal::default_precision(new_precision);
    }
    
    ~PrecisionContext() {
        mpreal::default_precision(old_precision_);
    }
    
    // Prevent copying
    PrecisionContext(const PrecisionContext&) = delete;
    PrecisionContext& operator=(const PrecisionContext&) = delete;
};

// Usage (RAII pattern)
{
    PrecisionContext ctx(500);  // Temporarily use 500 bits
    auto result = computeSomething();
}  // Automatically restores old precision
```

## Solution 4: Per-Object Precision (Most Flexible)

Add precision parameter to constructors:

```cpp
// include/polynomial_mp.h
class PolynomialMP {
    int precision_;
    
public:
    // Constructor with explicit precision
    PolynomialMP(const std::vector<unsigned int>& degrees,
                 const std::vector<mpreal>& coeffs,
                 int precision = 0)  // 0 = use default
        : precision_(precision > 0 ? precision : mpreal::default_precision()) {
        
        // Ensure all coefficients have correct precision
        for (const auto& c : coeffs) {
            if (c.precision() != precision_) {
                // Convert or warn
            }
        }
    }
    
    // Evaluate with guaranteed precision
    mpreal evaluate(const std::vector<mpreal>& params) const {
        PrecisionContext ctx(precision_);  // Use object's precision
        return evaluateImpl(params);
    }
};
```

## Comparison of Solutions

| Solution | Runtime Flexibility | Compile Time | Code Changes | Complexity |
|----------|-------------------|--------------|--------------|------------|
| **Solution 1: MPFR Runtime** | ✅ Excellent | ✅ Fast | ✅ None | ✅ Simple |
| Solution 2: Policy Pattern | ✅ Excellent | ✅ Fast | ⚠️ Moderate | ⚠️ Moderate |
| Solution 3: Context Manager | ✅ Good | ✅ Fast | ⚠️ Small | ✅ Simple |
| Solution 4: Per-Object | ✅ Excellent | ✅ Fast | ⚠️ Moderate | ⚠️ Moderate |

## Recommendation

**Use Solution 1 (MPFR Runtime Precision) as the primary approach.**

It requires:
- ✅ **Zero code changes** to the template design
- ✅ **Zero recompilation** when changing precision
- ✅ **Zero template instantiation overhead**
- ✅ Works perfectly with the private template design

**Add Solution 3 (Context Manager)** as a convenience utility for users.

**Add Solution 4 (Per-Object Precision)** if users need fine-grained control.

## Why This Works with Private Templates

The private template design is **type-based**, not **precision-based**:

```cpp
// Template is instantiated once for type 'mpreal'
template class PolynomialImpl<mpreal>;

// But mpreal itself supports runtime precision!
mpreal x(100);   // 100 bits
mpreal y(500);   // 500 bits
// Both are the same type 'mpreal', so same template instantiation
```

**Key point**: The template is instantiated for the **type** (`mpreal`), not for a specific precision. Since `mpreal` supports runtime precision, you get flexibility without recompilation.

## Implementation Plan

### Phase 2 Addition: Add Precision Utilities

When implementing Phase 2 (Type Definition Headers), add:

```cpp
// include/high_precision_types.h
#ifdef USE_MPFR_BACKEND
    using mpreal = boost::multiprecision::mpfr_float;
    
    // Convenience functions
    inline void setPrecision(int bits) {
        mpreal::default_precision(bits);
    }
    
    inline int getPrecision() {
        return mpreal::default_precision();
    }
#endif
```

### Phase 5 Addition: Add Context Manager

```cpp
// include/precision_context.h
#ifdef ENABLE_HIGH_PRECISION
class PrecisionContext {
    // ... (as shown in Solution 3)
};
#endif
```

## Example: Complete Adaptive Workflow

```cpp
#include "polynomial_solver.h"
#include "precision_context.h"

int main() {
    // Estimate precision needed
    int estimated_bits = estimatePrecisionNeeded(problem);
    
    // Set global precision
    setPrecision(estimated_bits);
    
    // Create and solve with estimated precision
    PolynomialMP poly(degrees, coeffs);
    auto roots = solver.solve(poly);
    
    // Refine specific roots with higher precision
    for (auto& root : roots) {
        PrecisionContext ctx(estimated_bits * 2);  // Double precision
        root = refineRoot(root, poly);
    }
    
    return 0;
}
```

**No recompilation needed when changing `estimated_bits`!**

## Important: Backend Differences

### MPFR Backend (Recommended for Runtime Precision)

```cpp
using mpreal = boost::multiprecision::mpfr_float;

// ✅ Runtime precision control
mpreal::default_precision(100);  // 100 bits
mpreal::default_precision(500);  // 500 bits - no recompilation!

// ✅ Per-object precision
mpreal x(0, 100);  // 100 bits
mpreal y(0, 500);  // 500 bits
```

**Advantages**:
- ✅ Arbitrary precision at runtime
- ✅ No recompilation needed
- ✅ Perfect for adaptive algorithms
- ✅ Fastest performance

### cpp_dec_float Backend (Fixed Precision)

```cpp
using mpreal = boost::multiprecision::cpp_dec_float_50;  // Fixed 50 digits

// ❌ Cannot change precision at runtime
// mpreal::default_precision(100);  // Does not exist!

// To change precision, must use different type:
using mpreal100 = boost::multiprecision::cpp_dec_float_100;  // Fixed 100 digits
```

**Limitations**:
- ❌ Precision fixed at compile time
- ❌ Need recompilation to change precision
- ❌ Need multiple template instantiations for multiple precisions

**Workaround for cpp_dec_float**:

If you must use cpp_dec_float and need multiple precisions, you can:

```cpp
// Define multiple precision levels at compile time
using mpreal50 = boost::multiprecision::cpp_dec_float_50;
using mpreal100 = boost::multiprecision::cpp_dec_float_100;

// Instantiate templates for each precision
template class PolynomialImpl<mpreal50>;
template class PolynomialImpl<mpreal100>;

// User selects at runtime
if (estimated_bits <= 166) {  // 50 decimal digits ≈ 166 bits
    PolynomialImpl<mpreal50> poly(degrees, coeffs50);
} else {
    PolynomialImpl<mpreal100> poly(degrees, coeffs100);
}
```

**Cost**: Longer compile time (multiple instantiations), but still no recompilation when switching between predefined precisions.

### Quadmath Backend (Fixed Precision)

```cpp
using quad = __float128;

// ❌ Fixed 128-bit precision
// Cannot change at runtime
```

**Limitation**: Always 128 bits (33-36 decimal digits).

### Recommendation

**For adaptive precision algorithms**: Use MPFR backend (`-DUSE_MPFR=ON -DUSE_GMP=ON`)

**For fixed precision**: Any backend works (cpp_dec_float is header-only, easier to deploy)

## Summary

The private template design **does not increase compile time** when changing precision at runtime, **if you use MPFR backend**.

**Key insight**:
- Template instantiation happens once per **type**
- MPFR's `mpreal` type supports **runtime precision**
- Therefore: One instantiation, infinite precisions! ✅

**Compile time comparison**:

| Scenario | MPFR Backend | cpp_dec_float Backend |
|----------|--------------|----------------------|
| Single precision | 1 instantiation | 1 instantiation |
| Runtime precision change | 0 recompilation | N/A (not supported) |
| Multiple precisions | 1 instantiation | N instantiations |
| Adaptive algorithm | ✅ Perfect | ⚠️ Workaround needed |

**Conclusion**: The private template design + MPFR backend = **zero compile-time overhead for runtime precision changes**.

