# MolecularIntegrals.jl

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.com/rpmuller/MolecularIntegrals.jl.svg?branch=master)](https://travis-ci.com/rpmuller/MolecularIntegrals.jl) | [![codecov.io](http://codecov.io/github/rpmuller/MolecularIntegrals.jl/coverage.svg?branch=master)](http://codecov.io/github/rpmuller/MolecularIntegrals.jl?branch=master) |


The goal of MolecularIntegrals.jl is to supply **fast** and
**hackable** one- and two-electron integrals for computational chemistry calculations. 
There are many excellent molecular integral packages available to Julia programmers, but few are written *in* Julia. This project will explore how fast we can make these integrals while maintaining a readable and hackable code base.

The code is released under a [MIT License](LICENSE.md).

## Starting points
MolecularIntegrals.jl will do its best to leverage the wonderful work done by many similar projects:

- [libints](https://github.com/evaleev/libint) and its Julia bindings [Lints.jl](https://github.com/FermiQC/Lints.jl)
- [Pyscf](https://github.com/pyscf/pyscf), the [libcint](https://github.com/sunqm/libcint) package.
- [JuliaChem.jl](https://github.com/davpoolechem/JuliaChem.jl)'s [JERI bindings](https://github.com/davpoolechem/JuliaChem.jl/tree/development/deps/src)
- [PyQuante](https://github.com/rpmuller/pyquante2)'s [python](https://github.com/rpmuller/pyquante2/tree/master/pyquante2/ints) and [c/cython](https://github.com/rpmuller/pyquante2/tree/master/cython) integrals, and the experimental [Julia version](https://github.com/rpmuller/pyquante2/tree/master/julia)


## Timing comparison
As a starting point to motivate the development, we will consider Table 5 from the paper [Libcint: an efficient general integral library for Gaussian basis functions](http://arxiv.org/abs/1412.0649) written by Qiming Sun, the author of Libcint and Pyscf, reporting the time required to compute electron repulsion integrals for ethane using different programs and different basis sets.

| Basis       | size | Psi4   | Molpro | Libcint w/o | w/SSE3 |
| ----------- | ---- | ------ | ------ | ----------- | ------ |
| 6-31G       | 30   | 0.10   | 0.09   | 0.09    | 0.07   |
| 6-311G**    | 72   | 0.64   | 0.49   | 0.49    | 0.41   |
| ANO         | 238  | 2527.6 | 51.13  | 53.59   | 37.78  |
| cc-pVDZ     | 58   | 0.45   | 0.34   | 0.24    | 0.21   |
| aug-cc-pVDZ | 100  | 1.87   | 1.18   | 1.02    | 0.85   |
| cc-pVTZ     | 144  | 4.98   | 4.82   | 2.65    | 2.05   |
| aug-cc-pVTZ | 230  | 26.03  | 23.12  | 12.40   | 9.27   |
| cc-pVQZ     | 290  | 81.24  | 65.12  | 31.51   | 22.60  |
| aug-cc-pVQZ | 436  | 444.23 | 324.29 | 151.04  | 107.24 |

### Preliminary Julia Timings for Ethane/6-31G
Compare to roughly 0.1 sec albeit on completely different hardware. This is mostly just to set a crude benchmark to compare coding improvements to.

| Basis   | size   | Huz    | HGP    |
| ------- | ------ | ------ | ------ |
| sto-3G  | 16     | 1.72   | 0.53   |
| 6-31G   | 30     | 6.62   | 2.12   |
| cc-pVDZ | 58     |        |        |

The HGP results hopefully still have lots of room for speedups.