# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MolecularIntegrals.jl is a Julia package for computing fast, hackable one- and two-electron integrals for computational chemistry. The goal is to provide **readable, pure-Julia implementations** that are competitive with C/Fortran libraries like libint and libcint.

The package implements multiple algorithms for electron repulsion integrals (ERIs):
- **Huzinaga (ERI.jl)**: Basic algorithm using B-arrays and Boys functions
- **HGP (HGP.jl)**: Head-Gordon-Pople recurrence relations - typically fastest for small-medium systems
- **Rys (Rys.jl)**: Rys polynomial quadrature - competitive for larger basis sets

## Running Tests

```bash
# Run all tests
julia --project=. test/runtests.jl

# Or using Julia's package manager
julia --project=.
julia> using Pkg
julia> Pkg.test()
```

Tests are defined in `test/runtests.jl` and use Julia's `@testset` framework. The test file defines standard test molecules (h2, h2o) at the top that are used throughout.

## Benchmarking

Performance benchmarks are in `test/timing.jl`:

```bash
# Run timing comparisons
julia --project=. test/timing.jl
```

This compares the three ERI methods (Huzinaga, HGP, Rys) across different basis sets and thread counts.

## Python Wrapper

The `python/` directory contains a Python wrapper via juliacall that provides a pyquante2-compatible interface.

### Python Environment Setup

The Python wrapper uses `uv` for environment management:

```bash
cd python
uv sync --extra dev    # Install all dependencies including pyquante2
uv run python examples/compare_with_pyquante2.py
uv run jupyter notebook examples/test_ints.ipynb
```

**Important:** The wrapper requires Julia dependencies to be specified in `python/juliapkg.json` because juliacall manages its own Julia environment separate from any local Julia installations.

Current required Julia dependencies in juliapkg.json:
- StaticArrays
- OffsetArrays
- DataStructures
- SpecialFunctions
- Interpolations

### Python API

The wrapper in `python/molecular_integrals.py` exports:
- `coulomb_repulsion()` - pyquante2-compatible primitive ERI function
- `Atom`, `Basis`, `PGBF`, `CGBF` - basis set classes
- `MolecularIntegrals` - high-level interface with `coulomb()`, `all_twoe_ints()`, `overlap()`, `kinetic()`

## Architecture

### Core Data Structures

**Atoms (Atoms.jl)**
- `Atom(atno, xyz::MVector{3,Float64})` - atomic number and coordinates (in Bohr)
- Coordinates use StaticArrays.MVector for efficient small vector operations

**Basis Functions (Basis.jl)**
- `PGBF` - Primitive Gaussian: exp(-α*r²) * x^I * y^J * z^K, with normalization
- `CGBF` - Contracted Gaussian: linear combination of PGBFs with coefficients
- `Shell` - Group of CGBFs sharing same center, exponents, and L quantum number
- `Basis` - Complete basis set with shell bookkeeping (ishell/mshell mappings)

**Shell Indexing**
- `shell_indices[L]` maps L quantum number to Cartesian (I,J,K) powers
- Example: L=1 (p-orbitals) → [(1,0,0), (0,1,0), (0,0,1)] for px, py, pz
- `m2ao` and `ao2m` provide bidirectional mapping between Cartesian indices and sequential AO numbering

### Integral Symmetry and Indexing

ERIs have 8-fold permutational symmetry: (ij|kl) = (ji|kl) = (ij|lk) = ...

**Compact storage:**
- `iiterator(n)` generates unique (i,j,k,l) indices with i≤j, k≤l, ij≤kl
- `iindex(i,j,k,l)` maps 4D indices to 1D packed array using triangle indexing
- Result: n×n×n×n tensor stored in ~n⁴/8 elements

**Thread parallelization:**
```julia
Threads.@threads for (i,j,k,l) in collect(iiterator(n))
    ints2e[iindex(i,j,k,l)] = ERI(bfs[i],bfs[j],bfs[k],bfs[l])
end
```

### Algorithm Implementations

**Basic Algorithm (ERI.jl):**
- Computes B-arrays via Bterm recurrence over angular momentum
- Uses Boys function Fγ(m,T) = ∫₀¹ t^(2m) exp(-Tt²) dt via incomplete gamma
- Triple nested loop over B-array indices with Boys function evaluation

**HGP (HGP.jl):**
- Separates into vertical (VRR) and horizontal (HRR) recurrence relations
- VRR: builds (a0|c0) integrals from (s-type|s-type)
- HRR: builds (ab|cd) from (a0|cd) using transfer relations
- `cvrr()` contracts primitive integrals at VRR stage (optimal efficiency)
- `chrr()` applies HRR to contracted integrals for final (ab|cd)

**Rys Quadrature (Rys.jl):**
- Uses Rys polynomial roots/weights for numerical integration
- Pre-allocates G arrays and root/weight buffers to avoid allocations
- Good for high angular momentum where recurrence relations become expensive

### Boys Function Evaluation (Boys.jl)

Multiple strategies for Fₘ(T):
- `Fm()`: Uses SpecialFunctions.gamma_inc (incomplete gamma function)
- `farray_recur()`: Downward recursion Fₘ = (2TFₘ₊₁ + e^(-T))/(2m+1)
- `farray()`: Switches to asymptotic expansion for T>30
- `Fm_asymp()`: √(π/2) × (2m-1)!! / (2T)^(m+0.5)

### Performance Considerations

**StaticArrays usage:**
- Atom coordinates and intermediate vectors use SVector/MVector (size 3)
- Avoids heap allocations for small fixed-size arrays
- Note: ERI.jl line 5 documents that StaticArrays didn't help for core ERI loops

**OffsetArrays:**
- B-arrays and other recurrence buffers use 0-based indexing to match mathematical notation
- Example: `B = OffsetArray(zeros(Float64, Imax), 0:(Imax-1))`

**Threading:**
- `all_twoe_ints()` uses `Threads.@threads` over unique integral indices
- Launch Julia with `julia -t 4` to use 4 threads

**Contraction strategies:**
- HGP contracts at VRR stage before HRR for optimal performance
- Avoids computing unnecessary primitive integrals

## Key Implementation Details

**Normalization:**
- PGBFs auto-normalize on creation via `normalize!()`
- CGBFs renormalize after each primitive is added with `addbf!()`
- Python wrapper must handle norms separately (pyquante2 compatibility)

**Gaussian Product Centers:**
- P = (αA + βB)/(α+β) for product of Gaussians centered at A and B
- Used throughout ERI calculations for reduced centers

**Basis Set Data (Data.jl):**
- Contains basis set definitions (STO-3G, 6-31G, cc-pVDZ, etc.)
- Accessed via `basis_data[name][atno]` → list of (symbol, [(expn, coef), ...])
- `build_basis(atoms, "sto3g")` constructs full basis from atom list

## Common Patterns

**Adding new basis sets:** Extend `basis_data` dictionary in Data.jl

**Testing new ERI methods:**
1. Implement function with signature `f(::CGBF, ::CGBF, ::CGBF, ::CGBF) -> Float64`
2. Add `all_twoe_ints(bfs, ERI=your_method)` call
3. Compare against reference: `@test isapprox(all_twoe_ints(bfs, ERI=new_method), all_twoe_ints(bfs), rtol=1e-7)`

**Performance profiling:**
```julia
using Profile
@profile all_twoe_ints(bfs)
Profile.print()
# Or use test/profile.jl for profiling specific methods
```

## Notes on Generated Code

- HGPgen.jl, HGPgen2.jl, HGPgen3.jl contain generated code for specific angular momentum cases
- HGPgen3.jl is 243KB - likely unrolled loops for common (ss|ss), (sp|sp), etc. combinations
- These provide optimized paths avoiding general recursion for frequently encountered integrals
