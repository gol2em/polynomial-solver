# Direct Contraction: Visual Comparison

## Current Approach (Incremental Restriction)

```
Initial State:
┌─────────────────────────────────────────┐
│ Node at depth 0                         │
│ Global box: [0, 1]                      │
│ Polynomials: p₀ (original)              │
│ Error: ε₀ = 0                           │
└─────────────────────────────────────────┘

Iteration 1: Contract to [0.1, 0.9]
┌─────────────────────────────────────────┐
│ Bounding box gives: [0.1, 0.9] (local)  │
│ Update global box: [0.1, 0.9]           │
│ Restrict: p₁ = p₀.restrict(0.1, 0.9)    │  ← 2 subdivisions, error ε₁
│ Error: ε₁ ≈ n·ε·||p₀||                   │
└─────────────────────────────────────────┘

Iteration 2: Contract to [0.2, 0.8] (in local [0,1])
┌─────────────────────────────────────────┐
│ Bounding box gives: [0.2, 0.8] (local)  │
│ Maps to global: [0.18, 0.82]            │
│ Restrict: p₂ = p₁.restrict(0.2, 0.8)    │  ← 2 subdivisions, error ε₂
│ Error: ε₂ ≈ ε₁ + n·ε·||p₁||              │
└─────────────────────────────────────────┘

Iteration 3: Contract to [0.3, 0.7] (in local [0,1])
┌─────────────────────────────────────────┐
│ Bounding box gives: [0.3, 0.7] (local)  │
│ Maps to global: [0.372, 0.628]          │
│ Restrict: p₃ = p₂.restrict(0.3, 0.7)    │  ← 2 subdivisions, error ε₃
│ Error: ε₃ ≈ ε₂ + n·ε·||p₂||              │
│       ≈ 3·n·ε·||p₀|| (accumulates!)      │
└─────────────────────────────────────────┘

After k iterations:
┌─────────────────────────────────────────┐
│ Global box: [A, B]                      │
│ Polynomials: pₖ                         │
│ Error: εₖ ≈ k·n·ε·||p₀||                │  ← LINEAR GROWTH!
└─────────────────────────────────────────┘
```

**Problem**: Error grows linearly with number of contractions!

---

## Proposed Approach (Direct Restriction)

```
Initial State:
┌─────────────────────────────────────────┐
│ Node at depth 0                         │
│ Global box: [0, 1]                      │
│ Polynomials: p (current)                │
│ Original: p₀ (stored, never modified)   │
│ Error: ε = 0                            │
└─────────────────────────────────────────┘

Iteration 1: Contract to [0.1, 0.9]
┌─────────────────────────────────────────┐
│ Bounding box gives: [0.1, 0.9] (local)  │
│ Update global box: [0.1, 0.9]           │
│ Restrict: p = p₀.restrict(0.1, 0.9)     │  ← 2 subdivisions from ORIGINAL
│ Error: ε ≈ n·ε·||p₀||                    │
└─────────────────────────────────────────┘

Iteration 2: Contract to [0.2, 0.8] (in local [0,1])
┌─────────────────────────────────────────┐
│ Bounding box gives: [0.2, 0.8] (local)  │
│ Maps to global: [0.18, 0.82]            │
│ Restrict: p = p₀.restrict(0.18, 0.82)   │  ← 2 subdivisions from ORIGINAL
│ Error: ε ≈ n·ε·||p₀|| (same!)            │
└─────────────────────────────────────────┘

Iteration 3: Contract to [0.3, 0.7] (in local [0,1])
┌─────────────────────────────────────────┐
│ Bounding box gives: [0.3, 0.7] (local)  │
│ Maps to global: [0.372, 0.628]          │
│ Restrict: p = p₀.restrict(0.372, 0.628) │  ← 2 subdivisions from ORIGINAL
│ Error: ε ≈ n·ε·||p₀|| (same!)            │
└─────────────────────────────────────────┘

After k iterations:
┌─────────────────────────────────────────┐
│ Global box: [A, B]                      │
│ Polynomials: p = p₀.restrict(A, B)      │
│ Original: p₀ (unchanged)                │
│ Error: ε ≈ n·ε·||p₀||                    │  ← CONSTANT!
└─────────────────────────────────────────┘
```

**Benefit**: Error is independent of number of contractions!

---

## Key Difference Illustrated

### Current (Incremental):
```
p₀ ──restrict──> p₁ ──restrict──> p₂ ──restrict──> p₃ ──restrict──> ... ──> pₖ
   ε₁            ε₂            ε₃            ε₄                        εₖ

Total error: ε₁ + ε₂ + ε₃ + ... + εₖ ≈ k·n·ε·||p₀||
```

### Proposed (Direct):
```
p₀ ──restrict──> p (iteration 1)
   ε₁

p₀ ──restrict──> p (iteration 2)
   ε₂

p₀ ──restrict──> p (iteration 3)
   ε₃

...

p₀ ──restrict──> p (iteration k)
   εₖ

Total error: εₖ ≈ n·ε·||p₀|| (independent of k!)
```

---

## Computational Cost Comparison

### Current Approach:
```
Iteration 1: Compute bounding box + restrict (2 subdivisions)
Iteration 2: Compute bounding box + restrict (2 subdivisions)
Iteration 3: Compute bounding box + restrict (2 subdivisions)
...
Iteration k: Compute bounding box + restrict (2 subdivisions)

Total: k × (bounding box + 2 subdivisions)
```

### Proposed Approach:
```
Iteration 1: Compute bounding box + restrict from original (2 subdivisions)
Iteration 2: Compute bounding box + restrict from original (2 subdivisions)
Iteration 3: Compute bounding box + restrict from original (2 subdivisions)
...
Iteration k: Compute bounding box + restrict from original (2 subdivisions)

Total: k × (bounding box + 2 subdivisions)
```

**Cost: EXACTLY THE SAME!**

The key insight: Both approaches do 2 subdivisions per contraction. The difference is:
- Current: Restrict from current [0,1] to [a,b]
- Proposed: Restrict from original [0,1] to global [A,B]

Both are 2 subdivisions, but restricting from original eliminates error accumulation!

---

## Memory Usage Comparison

### Current Approach:
```
SubdivisionNode:
  - box_lower, box_upper: 2×dim doubles
  - depth: 1 unsigned int
  - polys: num_equations × polynomial

Memory per node: ~(2×dim + num_equations×poly_size)
```

### Proposed Approach:
```
SubdivisionNode:
  - box_lower, box_upper: 2×dim doubles
  - depth: 1 unsigned int
  - polys: num_equations × polynomial
  - original_polys: num_equations × polynomial  ← NEW

Memory per node: ~(2×dim + 2×num_equations×poly_size)
```

**Memory increase**: 2× polynomial storage

**For degree 10, 2D, 2 equations**:
- Polynomial size: ~11×11 = 121 doubles
- Current: 2×2 + 2×121 = 246 doubles ≈ 2 KB
- Proposed: 2×2 + 4×121 = 488 doubles ≈ 4 KB
- Increase: 2 KB per node

**For 1000 nodes**: 2 MB total (negligible!)

---

## Summary

| Aspect | Current | Proposed | Change |
|--------|---------|----------|--------|
| **Error** | k·n·ε·||p₀|| | n·ε·||p₀|| | **k× better** |
| **Cost** | 2 subdivisions/contraction | 2 subdivisions/contraction | **Same** |
| **Memory** | 1× polynomials | 2× polynomials | **2× more** |
| **Complexity** | Simple | Simple | **Same** |

**Verdict**: Pure accuracy improvement with negligible memory cost!

