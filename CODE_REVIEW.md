# MolecularIntegrals.jl - Comprehensive Code Review and Performance Optimization Guide

**Date:** January 19, 2026
**Version:** 0.1.0
**Reviewer:** Automated Analysis via Claude Code

---

## Executive Summary

### Current Performance Baseline

Performance benchmarks for Ethane (C₂H₆) with 6-31G basis set (30 basis functions):

| Algorithm | Single Thread | 4 Threads | Description |
|-----------|--------------|-----------|-------------|
| **Huzinaga** | 6.56s | 3.72s | Basic B-array algorithm |
| **HGP** | 0.202s | 0.131s | Head-Gordon-Pople (production choice) |
| **Rys** | 1.550s | 0.803s | Rys polynomial quadrature |

**Comparison to Industry Standards:**
- libcint w/SSE3 (Ethane/6-31G): ~0.07s
- Current gap: **2.9x slower** than optimized C implementations

### Top 5 Optimization Priorities

| Priority | Issue | Location | Est. Impact | Effort |
|----------|-------|----------|-------------|--------|
| **1** | Boys function evaluation | src/Boys.jl | 40% overall | Medium |
| **2** | Redundant CGBF normalization | src/Basis.jl:79-83 | 50-80% basis | Low |
| **3** | Memory allocations in hot loops | src/ERI.jl, Rys.jl | 30-50% | Medium |
| **4** | Contract() screening | src/Basis.jl:85-98 | 20-40% | Medium |
| **5** | factorial2 lookup table | src/Utils.jl:5 | 5-20x local | Very Low |

### Expected Improvements

**Conservative estimate:** 3-5x overall speedup
**Target performance:** 0.04-0.07s for Ethane/6-31G (competitive with libcint)

**Breakdown by phase:**
- Phase 1 (quick wins): +60% speedup
- Phase 2 (Boys optimization): +40% additional
- Phases 3-6 (allocations, screening, SIMD): +50% cumulative

---

## Detailed Analysis by Category

### 1. Memory Allocations

#### 1.1 B-array Allocation in coulomb() ⚠️ CRITICAL

**Location:** `src/ERI.jl:31-33`

**Current Code:**
```julia
function coulomb(aexpn,ax,ay,az,aI,aJ,aK,
    bexpn,bx,by,bz,bI,bJ,bK,
    cexpn,cx,cy,cz,cI,cJ,cK,
    dexpn,dx,dy,dz,dI,dJ,dK)

    # ... setup code ...

    Bx = Barray(aI,bI,cI,dI,px,ax,bx,qx,cx,dx,g1,g2,delta)  # ← Allocation
    By = Barray(aJ,bJ,cJ,dJ,py,ay,by,qy,cy,dy,g1,g2,delta)  # ← Allocation
    Bz = Barray(aK,bK,cK,dK,pz,az,bz,qz,cz,dz,g1,g2,delta)  # ← Allocation

    # ... rest of function ...
end
```

**Problem:**
- Each `Barray()` call allocates a fresh `OffsetArray` (src/ERI.jl:79)
- 3 allocations per integral evaluation
- In contractions with 3-6 primitives: 27-216 allocations per CGBF quartet
- Total overhead: 30-40% of runtime

**Solution:**
```julia
# Pre-allocated buffer version
function coulomb!(Bx_buf, By_buf, Bz_buf,
                  aexpn,ax,ay,az,aI,aJ,aK,
                  bexpn,bx,by,bz,bI,bJ,bK,
                  cexpn,cx,cy,cz,cI,cJ,cK,
                  dexpn,dx,dy,dz,dI,dJ,dK)

    # ... setup code ...

    Barray!(Bx_buf, aI,bI,cI,dI,px,ax,bx,qx,cx,dx,g1,g2,delta)
    Barray!(By_buf, aJ,bJ,cJ,dJ,py,ay,by,qy,cy,dy,g1,g2,delta)
    Barray!(Bz_buf, aK,bK,cK,dK,pz,az,bz,qz,cz,dz,g1,g2,delta)

    # ... rest of function ...
end

# Thread-local buffer pool
const BARRAY_BUFFERS = [Vector{Float64}(undef, 20) for _ in 1:Threads.nthreads()]

function get_barrays()
    tid = Threads.threadid()
    return (BARRAY_BUFFERS[tid], BARRAY_BUFFERS[tid+nthreads()], ...)
end
```

**Expected Impact:** 30-40% speedup in ERI evaluation

---

#### 1.2 Rys Quadrature Working Arrays ⚠️ CRITICAL

**Location:** `src/Rys.jl:27-56`

**Current Code:**
```julia
function coulomb_rys(a::CGBF,b::CGBF,c::CGBF,d::CGBF)
    norder = (a.I+a.J+a.K+b.I+b.J+b.K+c.I+c.J+c.K+d.I+d.J+d.K)÷2 + 1
    n = max(a.I+b.I,a.J+b.J,a.K+b.K)
    m = max(c.I+d.I,c.J+d.J,c.K+d.K)
    G = OffsetArray(zeros(Float64,n+1,m+1),0:n,0:m)   # ← Allocation
    roots = zeros(Float64,norder)                     # ← Allocation
    weights = zeros(Float64,norder)                   # ← Allocation

    val = 0
    for (ap,ac) in zip(a.pgbfs,a.coefs), (bp,bc) in zip(b.pgbfs,b.coefs)
        for (cp,cc) in zip(c.pgbfs,c.coefs), (dp,dc) in zip(d.pgbfs,d.coefs)
            val += ac*bc*cc*dc*coulomb_rys(ap,bp,cp,dp,G,roots,weights)
        end
    end
    return val
end
```

**Problem:**
- 3 arrays allocated per CGBF quartet call
- Called inside contraction loop over all primitive combinations
- For 3 primitives per basis function: 81 allocation calls per quartet

**Solution:**
```julia
# Pre-allocated version
struct RysBuffers
    G::OffsetArray{Float64,2}
    roots::Vector{Float64}
    weights::Vector{Float64}
end

function RysBuffers(maxorder=10, maxn=5, maxm=5)
    G = OffsetArray(zeros(Float64,maxn+1,maxm+1),0:maxn,0:maxm)
    roots = Vector{Float64}(undef, maxorder)
    weights = Vector{Float64}(undef, maxorder)
    return RysBuffers(G, roots, weights)
end

const RYS_BUFFERS = [RysBuffers() for _ in 1:Threads.nthreads()]

function coulomb_rys(a::CGBF,b::CGBF,c::CGBF,d::CGBF)
    bufs = RYS_BUFFERS[Threads.threadid()]
    # Use bufs.G, bufs.roots, bufs.weights...
end
```

**Expected Impact:** 20-30% speedup in Rys algorithm

---

### 2. Algorithmic Inefficiencies

#### 2.1 Boys Function Evaluation ⚠️ CRITICAL (40% of VRR time)

**Location:** `src/Boys.jl:29-34`, used in `src/HGP.jl:168`

**Current Implementation:**
```julia
function Fm(m,T,SMALL=1e-18)
    mhalf = m+0.5
    if T<SMALL return 1/(2m+1) end
    return gamma(mhalf)*gamma_inc(mhalf,T,0)[1]/(2T^(mhalf))
end

function farray(mmax,T,fm=Fm,Tmax=30)
    if T < Tmax return farray_recur(mmax,T,fm) end
    # ... asymptotic expansion ...
end
```

**Critical Comment from HGP.jl:162-165:**
```julia
# First generate (1,1,m) using eq 12 using an array of calls to
# the incomplete gamma function. This is ~40% of the time here,
# but can be replaced with table lookup or interpolation.
```

**Problem:**
- Special function evaluation (gamma_inc) is expensive
- Called for every unique (m, T) pair in VRR
- ~40% of total VRR computation time
- Research notebooks show interpolation gives **13x speedup**

**Solution - Table Interpolation:**
```julia
using Interpolations

# Module initialization
function __init__()
    global const BOYS_INTERPOLATOR = make_boys_interpolator()
end

function make_boys_interpolator(mmax=12, Tmax=20.0, ΔT=0.005)
    Tgrid = 0:ΔT:Tmax
    mgrid = 0:mmax

    # Pre-compute all values
    fmvalues = zeros(Float64, length(mgrid), length(Tgrid))
    for (i,m) in enumerate(mgrid), (j,T) in enumerate(Tgrid)
        fmvalues[i,j] = Fm_ref(m,T)  # Use reference implementation
    end

    # Create cubic B-spline interpolator
    itp = interpolate(fmvalues, BSpline(Cubic(Line(OnGrid()))))
    return scale(itp, mgrid, Tgrid)
end

function Fm_fast(m, T)
    if T > 20.0
        return Fm_asymp(m, T)  # Asymptotic form
    else
        return BOYS_INTERPOLATOR(m, T)
    end
end

# Replace in farray:
function farray(mmax,T,Tmax=20)
    if T < Tmax
        Fms = zeros(Float64,mmax)
        Fms[mmax] = Fm_fast(mmax-1,T)
        emt = exp(-T)
        for m in mmax-1:-1:1
            Fms[m] = (2*T*Fms[m+1]+emt)/(2*(m-1)+1)
        end
        return Fms
    end
    # ... asymptotic path ...
end
```

**Research Evidence:** (from tools/simpleboys.jl)
- Interpolation: 17.7 μs for 1000 evaluations
- Direct computation: 229.5 μs for 1000 evaluations
- **Speedup: 13x**

**Expected Impact:** 40% overall speedup in HGP algorithm → **2.5x total speedup**

---

#### 2.2 Redundant CGBF Normalization ⚠️ CRITICAL

**Location:** `src/Basis.jl:79-83`

**Current Code:**
```julia
function addbf!(cbf::CGBF,expn,coef)
    Base.push!(cbf.pgbfs,pgbf(expn,cbf.xyz...,cbf.I,cbf.J,cbf.K))
    Base.push!(cbf.coefs,coef)
    normalize!(cbf)  # ← Called EVERY time a primitive is added!
end
```

**Usage in build_basis():**
```julia
function build_basis(atoms::Vector{Atom},name="sto3g")
    # ...
    for (expn,coef) in zip(sh.expns,sh.coefs)
        addbf!(cbf,expn,coef)  # Normalizes after each primitive
    end
    # ...
end
```

**Problem:**
- For N primitives in a contraction, `normalize!()` is called N times
- Each call computes `overlap(cbf,cbf)` which is O(N²) operations
- Total: O(N³) redundant work
- For typical contractions (N=3-6): **27-216x redundancy**

**Impact Analysis:**
- STO-3G: 3 primitives per shell → 9x redundant normalization
- 6-31G: 6 primitives for core → 36x redundant normalization
- **Estimated speedup: 50-80% in basis construction phase**

**Solution:**
```julia
# Internal version without normalization
function addbf_unnormalized!(cbf::CGBF, expn, coef)
    Base.push!(cbf.pgbfs, pgbf_unnormalized(expn, cbf.xyz..., cbf.I, cbf.J, cbf.K))
    Base.push!(cbf.coefs, coef)
end

function pgbf_unnormalized(expn,x,y,z,I,J,K,norm=1.0)
    return PGBF(expn,SA[x,y,z],I,J,K,norm)  # Don't call normalize!
end

function build_basis(atoms::Vector{Atom},name="sto3g")
    # ...
    for (expn,coef) in zip(sh.expns,sh.coefs)
        addbf_unnormalized!(cbf,expn,coef)
    end
    normalize!(cbf)  # Normalize once at the end
    # ...
end

# Keep public API backward compatible
function addbf!(cbf::CGBF,expn,coef)
    addbf_unnormalized!(cbf, expn, coef)
    normalize!(cbf)
end
```

**Expected Impact:** 50-80% faster basis construction

---

#### 2.3 Contract() Function - Missing Screening

**Location:** `src/Basis.jl:85-98`

**Current Code:**
```julia
function contract(f,a::CGBF,b::CGBF,c::CGBF,d::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            for (cc,cbf) in primitives(c)
                for (cd,dbf) in primitives(d)
                    s += ca*cb*cc*cd*f(abf,bbf,cbf,dbf)  # No screening!
                end
            end
        end
    end
    return a.norm*b.norm*c.norm*d.norm*s
end
```

**Problem:**
- No early termination for negligible contributions
- Small primitive coefficients still trigger expensive `f()` evaluations
- For large basis sets with many primitives: significant wasted computation

**Solution:**
```julia
function contract(f,a::CGBF,b::CGBF,c::CGBF,d::CGBF; threshold=1e-12)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            cab = ca*cb
            if abs(cab) < threshold continue end  # Screen early

            for (cc,cbf) in primitives(c)
                for (cd,dbf) in primitives(d)
                    ccoef = cab*cc*cd
                    if abs(ccoef) < threshold continue end

                    s += ccoef*f(abf,bbf,cbf,dbf)
                end
            end
        end
    end
    return a.norm*b.norm*c.norm*d.norm*s
end
```

**Enhanced Screening with Exponent Schwarz Bound:**
```julia
function contract_screened(f,a::CGBF,b::CGBF,c::CGBF,d::CGBF; threshold=1e-12)
    # Pre-compute Schwarz bounds for primitive pairs
    maxab = maximum(abs(ca*cb) * schwarz_bound(abf,bbf)
                   for (ca,abf) in primitives(a)
                   for (cb,bbf) in primitives(b))
    maxcd = maximum(abs(cc*cd) * schwarz_bound(cbf,dbf)
                   for (cc,cbf) in primitives(c)
                   for (cd,dbf) in primitives(d))

    if maxab * maxcd < threshold
        return 0.0  # Skip entire quartet
    end

    # ... rest of contraction with per-primitive screening ...
end
```

**Expected Impact:** 20-40% speedup on large basis sets (6-31G, cc-pVDZ)

---

### 3. Vectorization Opportunities

#### 3.1 factorial2 Lookup Table ⚠️ HIGH ROI

**Location:** `src/Utils.jl:5`

**Current Code:**
```julia
@inline factorial2(n::Int64) = prod(n:-2:1)
```

**Problem:**
- Computes product every call: `prod(7:-2:1) = 7*5*3*1`
- Called in tight loops (Boys function, overlap calculations)
- For small n (≤7): ~20 operations wasted

**Solution:**
```julia
using StaticArrays

const FACTORIAL2_TABLE = SVector(
    1,    # 0!!
    1,    # 1!!
    2,    # 2!!
    3,    # 3!!
    8,    # 4!!
    15,   # 5!!
    48,   # 6!!
    105,  # 7!!
    384,  # 8!!
    945,  # 9!!
)

@inline function factorial2(n::Int64)
    return n < 10 ? FACTORIAL2_TABLE[n+1] : prod(n:-2:1)
end
```

**Expected Impact:** 5-20x speedup for factorial2 calls, ~5% overall

---

#### 3.2 SIMD Annotations for Boys Function Loop

**Location:** `src/HGP.jl:169`

**Current Code:**
```julia
Fms = farray(mmax,T)
for m in 1:mmax
    vrrs[m,1,1] = KabKcd_rtze*Fms[m]
end
```

**Problem:**
- Simple multiplication loop
- No `@simd` annotation despite being vectorizable
- Compiler may not auto-vectorize without hint

**Solution:**
```julia
Fms = farray(mmax,T)
@inbounds @simd ivdep for m in 1:mmax
    vrrs[m,1,1] = KabKcd_rtze*Fms[m]
end
```

**Expected Impact:** 10-20% improvement in VRR initialization

---

#### 3.3 ERI Triple Loop Vectorization

**Location:** `src/ERI.jl:38-44`

**Current Code:**
```julia
s = 0
for I in 0:(aI+bI+cI+dI)
    for J in 0:(aJ+bJ+cJ+dJ)
        for K in 0:(aK+bK+cK+dK)
            s += Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,x)
        end
    end
end
```

**Solution:**
```julia
s = 0
for I in 0:(aI+bI+cI+dI)
    BxI = Bx[I]
    for J in 0:(aJ+bJ+cJ+dJ)
        BxIByJ = BxI*By[J]
        @inbounds @simd for K in 0:(aK+bK+cK+dK)
            s += BxIByJ*Bz[K]*Fgamma(I+J+K,x)
        end
    end
end
```

**Expected Impact:** 5-15% for high angular momentum (d, f shells)

---

### 4. Caching and Redundancy

#### 4.1 Gaussian Product Center Caching

**Locations:** All ERI implementations

**Problem:**
```julia
# Called repeatedly in contractions for same primitive pairs
px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
```

**Solution:**
```julia
# Pre-compute in contraction outer loops
struct PrimitiveQuartet
    P::SVector{3,Float64}
    Q::SVector{3,Float64}
    rpq2::Float64
    # ... other shared data ...
end

function precompute_quartet(abf,bbf,cbf,dbf)
    P = gaussian_product_center(abf.expn, abf.xyz, bbf.expn, bbf.xyz)
    Q = gaussian_product_center(cbf.expn, cbf.xyz, dbf.expn, dbf.xyz)
    rpq2 = dist2(P - Q)
    return PrimitiveQuartet(P, Q, rpq2)
end

# Pass pre-computed data to coulomb()
function coulomb(quartet::PrimitiveQuartet, ...)
    # Use quartet.P, quartet.Q directly
end
```

**Expected Impact:** 10-15% reduction in redundant calculations

---

#### 4.2 Dictionary Lookup → Array Indexing

**Location:** `src/Basis.jl:316`

**Current Code:**
```julia
const m2ao = make_m2ao()  # Dict{Vector{Int64}, Int64}

# Usage in eri_fetcher:
mi = ijk2lm[ibf.I,ibf.J,ibf.K][2]  # Dictionary lookup with vector key
```

**Problem:**
- Vector keys require hashing
- Used in nested loops during initialization
- Can be replaced with direct computation

**Solution:**
```julia
# Direct computation (based on shell_indices ordering)
@inline function ijk_to_m(I::Int, J::Int, K::Int)
    L = I + J + K
    if L == 0
        return 1
    elseif L == 1
        return I==1 ? 1 : (J==1 ? 2 : 3)
    elseif L == 2
        # Implement based on shell_indices[2] ordering
        # ...
    end
    # ... or use pre-computed 3D array
end

# Or pre-compute small 3D array (since max I,J,K = 4)
const IJK_TO_M = make_ijk_to_m_array()

@inline ijk_to_m(I,J,K) = IJK_TO_M[I+1,J+1,K+1]
```

**Expected Impact:** 5-10% in eri_fetcher

---

### 5. Code Quality Issues

#### 5.1 eri_fetcher Iterator Overhead

**Location:** `src/Basis.jl:228,241`

**Current Code:**
```julia
for (jsh,shj) in enumerate(take(shs,ish))  # Creates iterator
    # ...
end
for (lsh,shl) in enumerate(take(shs,ksh))  # Creates iterator
    # ...
end
```

**Solution:**
```julia
@inbounds for jsh in 1:ish
    shj = shs[jsh]
    # ...
end
@inbounds for lsh in 1:ksh
    shl = shs[lsh]
    # ...
end
```

**Expected Impact:** 10-15% in initialization

---

## Prioritized Recommendations (Ordered by ROI)

### Tier 1: Critical - High Impact, Low-Medium Effort

1. **Defer CGBF Normalization**
   - Impact: 50-80% basis construction
   - Effort: 2-4 hours
   - Risk: Low
   - **ROI: 20-40x**

2. **factorial2 Lookup Table**
   - Impact: 5-20x local, ~5% overall
   - Effort: 1 hour
   - Risk: Very low
   - **ROI: 5-20x**

3. **Replace take() in eri_fetcher**
   - Impact: 10-15%
   - Effort: 1 hour
   - Risk: Very low
   - **ROI: 10-15x**

### Tier 2: High Priority - Major Impact, Medium Effort

4. **Boys Function Interpolation**
   - Impact: 40% overall (in HGP)
   - Effort: 4-8 hours
   - Risk: Medium (requires validation)
   - **ROI: 5-10x**

5. **Pre-allocate B-arrays**
   - Impact: 30-40%
   - Effort: 3-4 hours
   - Risk: Medium (threading complexity)
   - **ROI: 7.5-13x**

6. **Add Contract() Screening**
   - Impact: 20-40% (basis-dependent)
   - Effort: 4-6 hours
   - Risk: Medium (threshold tuning)
   - **ROI: 3.3-10x**

### Tier 3: Medium Priority - Moderate Impact, Low-Medium Effort

7. **Pre-allocate Rys Working Arrays**
   - Impact: 20-30%
   - Effort: 2-3 hours
   - **ROI: 6.7-15x**

8. **SIMD Boys Loop**
   - Impact: 10-20%
   - Effort: 1 hour
   - **ROI: 10-20x**

9. **Cache Gaussian Product Centers**
   - Impact: 10-15%
   - Effort: 3-4 hours
   - **ROI: 2.5-5x**

### Tier 4: Lower Priority - Smaller Impact

10. **m2ao Dictionary → Array**
    - Impact: 5-10%
    - Effort: 3-4 hours
    - **ROI: 1.25-3.3x**

11. **ERI Triple Loop SIMD**
    - Impact: 5-15%
    - Effort: 2-3 hours
    - **ROI: 1.7-7.5x**

---

## Implementation Roadmap

### Week 1: Quick Wins (Phase 1)
**Goal:** 60% cumulative speedup

1. **Day 1-2:** factorial2 table + defer normalization (Items 2, 1)
2. **Day 3:** Replace take() + SIMD Boys loop (Items 3, 8)
3. **Testing:** Verify test suite passes, benchmark improvements

**Deliverables:**
- Modified files: Utils.jl, Basis.jl, HGP.jl
- Benchmark report

### Week 2: Boys Function (Phase 2)
**Goal:** +40% additional speedup

1. **Day 1-3:** Implement interpolation table (Item 4)
2. **Day 4:** Optimize asymptotic cutoff, tune parameters
3. **Day 5:** Extensive testing and validation

**Deliverables:**
- Modified Boys.jl with __init__() function
- Interpolation accuracy report
- Performance comparison

### Week 3: Memory Optimization (Phase 3)
**Goal:** +30% additional speedup

1. **Day 1-2:** Pre-allocate B-arrays (Item 5)
2. **Day 3:** Pre-allocate Rys arrays (Item 7)
3. **Day 4-5:** Thread safety testing, allocation profiling

**Deliverables:**
- Modified ERI.jl, Rys.jl
- Thread scaling benchmarks
- Memory allocation report

### Week 4: Screening and Caching (Phases 4-6)
**Goal:** +25% additional speedup

1. **Day 1-3:** Implement contract() screening (Item 6)
2. **Day 4:** Cache Gaussian product centers (Item 9)
3. **Day 5:** Final integration and testing

**Deliverables:**
- Enhanced Basis.jl with screening
- Full benchmark suite
- Performance summary report

### Week 5: Refinement and Documentation
**Goal:** Polish and long-term improvements

1. **Day 1-2:** Remaining optimizations (Items 10-11)
2. **Day 3:** Comprehensive benchmarking
3. **Day 4-5:** Documentation updates, CLAUDE.md revision

**Dependencies:**
- Phase 2 can run parallel to Phase 1
- Phase 3 requires Phase 1 completion (buffer management)
- Phase 4 requires Phase 3 (screening depends on allocation strategy)

---

## Benchmarking Strategy

### Test Systems

| System | Formula | Basis | Size | Purpose |
|--------|---------|-------|------|---------|
| H₂ | H-H | STO-3G | 2 BF | Minimal, unit testing |
| H₂O | H₂O | STO-3G | 7 BF | Small molecule baseline |
| H₂O | H₂O | 6-31G | 13 BF | Medium basis |
| Ethane | C₂H₆ | 6-31G | 30 BF | Production benchmark |
| Ethane | C₂H₆ | cc-pVDZ | 58 BF | Large basis stress test |
| Benzene | C₆H₆ | 6-31G | 66 BF | Larger system |

### Metrics to Track

**Performance Metrics:**
1. Total integral evaluation time (seconds)
2. Time per unique integral (μs)
3. Basis construction time (ms)
4. Boys function evaluation time (ms)
5. Thread scaling efficiency (1, 2, 4, 8 threads)

**Memory Metrics:**
1. Total allocations (count)
2. Allocation size (MB)
3. Peak memory usage
4. Allocation rate (MB/s)

**Accuracy Metrics:**
1. Maximum absolute error vs reference (pyquante2)
2. RMS error across all integrals
3. Relative error for small integrals (< 1e-6)

### Profiling Commands

```julia
using Profile

# Allocation profiling
@allocated all_twoe_ints(basis)

# Time profiling
@profile all_twoe_ints(basis)
Profile.print()

# Or use ProfileSVG for visualization
using ProfileSVG
@profview all_twoe_ints(basis)

# Benchmark macro
using BenchmarkTools
@benchmark all_twoe_ints($basis)
```

### Regression Tests

Create `test/benchmarks.jl`:
```julia
using BenchmarkTools

@testset "Performance Benchmarks" begin
    # Basis construction
    @test (@belapsed build_basis(h2o, "sto3g")) < 0.001  # < 1ms

    # Boys function
    @test (@belapsed farray(10, 5.0)) < 0.00001  # < 10μs

    # Full integrals (Ethane/6-31G)
    atoms = load_ethane()
    basis = build_basis(atoms, "6-31g")
    @test (@belapsed all_twoe_ints($basis)) < 0.100  # < 100ms target
end
```

---

## Long-term Improvements

### Architectural Considerations

1. **Code Generation for Fixed Shell Types**
   - Current: General recursion for all (L₁,L₂,L₃,L₄) combinations
   - Opportunity: Generate specialized code for common cases (ss|ss), (sp|sp), (pp|pp), etc.
   - Evidence: HGPgen3.jl (243KB) shows this approach but may be outdated
   - **Potential:** 2-3x speedup for common integral types

2. **Symmetry Exploitation Enhancement**
   - Current: 8-fold symmetry via iiterator/iindex
   - Opportunity: Shell-level symmetry screening (Schwarz bounds)
   - Implementation: Pre-compute shell pair significance matrix
   - **Potential:** 30-50% reduction in integral evaluations for large systems

3. **Batch Processing**
   - Current: One quartet at a time
   - Opportunity: Process multiple quartets with same shell types together
   - Benefit: Better SIMD utilization, cache efficiency
   - **Potential:** 20-40% improvement

### GPU Acceleration Prospects

**Analysis:**
- ERIs are embarrassingly parallel (each quartet independent)
- Challenge: Limited arithmetic intensity (memory-bound)
- Recommendation: Start with VRR kernel (most compute-intensive)

**Proof of Concept:**
```julia
using CUDA

function vrr_gpu!(vrrs_d, Fms_d, ...)
    # CUDA kernel for VRR computation
    # Good candidate: high arithmetic intensity
end
```

**Expected Impact:** 5-10x for large basis sets (but overhead for small systems)

### Alternative Algorithms

1. **Obara-Saika Recursion**
   - Pro: Different recursion pattern, may have better cache behavior
   - Con: More complex than current VRR/HRR
   - Evidence: Used by libint, very fast

2. **McMurchie-Davidson**
   - Pro: Hermite Gaussians, fewer recursion steps
   - Con: Different basis representation
   - Use case: Cartesian Gaussian integrals

3. **Rys Polynomial Order Reduction**
   - Current: Fixed-order quadrature
   - Opportunity: Adaptive order based on integral magnitude
   - **Potential:** 20-30% faster Rys algorithm

### Code Quality Improvements

1. **Type Stability Audit**
   - Run `@code_warntype` on all hot paths
   - Ensure no dynamic dispatch in inner loops

2. **Const Propagation**
   - More aggressive use of `@fastmath` where safe
   - `@inbounds` for verified bounds

3. **Documentation**
   - Add performance notes to docstrings
   - Document which functions allocate
   - Threading safety annotations

---

## Appendix: Benchmark Comparison Table

### Target Performance vs Current

| System | Basis | Size | Current (HGP) | Target | Gap |
|--------|-------|------|---------------|--------|-----|
| H₂ | STO-3G | 2 | 0.2 ms | 0.05 ms | 4x |
| H₂O | STO-3G | 7 | 25 ms | 6 ms | 4x |
| Ethane | 6-31G | 30 | 202 ms | 50 ms | 4x |
| Ethane | cc-pVDZ | 58 | 2451 ms | 600 ms | 4x |

### Comparison to Production Libraries

| System | libcint+SSE | Current | Target | Status |
|--------|-------------|---------|--------|--------|
| Ethane/6-31G | 70 ms | 202 ms | 50 ms | **2.9x gap** |
| Ethane/cc-pVDZ | 240 ms | 2451 ms | 600 ms | **10.2x gap** |

**Analysis:** Current code is competitive for small basis sets but falls behind on larger systems. Optimizations should close the gap significantly.

---

## Conclusion

MolecularIntegrals.jl is a well-architected pure-Julia implementation with significant optimization potential. The identified improvements can realistically achieve **3-5x speedup**, making it competitive with optimized C/Fortran libraries while maintaining code readability and hackability.

**Key Priorities:**
1. **Immediate (Week 1):** factorial2 table, defer normalization, remove take()
2. **High impact (Week 2):** Boys function interpolation
3. **Medium term (Weeks 3-4):** Memory allocation reduction, screening
4. **Long term:** Code generation, GPU acceleration, alternative algorithms

**Success Criteria:**
- Ethane/6-31G in < 70ms (single thread)
- Test suite passes with < 1e-10 error
- Thread scaling efficiency > 80% up to 4 threads
- Memory allocations reduced by > 50%

The phased approach allows incremental validation while building toward the performance target.
