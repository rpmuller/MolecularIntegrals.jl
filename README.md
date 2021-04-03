# MolecularIntegrals.jl

The goal of MolecularIntegrals.jl is to supply **fast** and
**hackable** one- and two-electron integrals. There are a lot of
excellent molecular integral packages in Julia, but few are written
*in* Julia. This project will explore how fast we can make these integrals while maintaining a readable and hackable code base.

## Starting points
MolecularIntegrals.jl will do its best to leverage the wonderful work done by many similar projects:

- [libints](https://github.com/evaleev/libint) and its Julia bindings [Lints.jl](https://github.com/FermiQC/Lints.jl)
- [Pyscf](https://github.com/pyscf/pyscf), the [libcint](https://github.com/sunqm/libcint) package, and the Julia bindings ???
- JuliaChem.jl's JERI bindings
- PyQuante's python and c integrals, and the experimental Julia version

As a starting point to motivate the development, we will consider Table 5 from the paper [Libcint: an efficient general integral library for Gaussian basis functions](http://arxiv.org/abs/1412.0649) written by Qiming Sun, the author of Libcint and Pyscf, reporting the time required to compute electron repulsion integrals for ethane using different programs and different basis sets.

| Basis       | size | Psi4   | Molpro | Libcint w/o | w/SSE3 |
|-------------|------|--------|--------|-------------|-- -----|
| 6-31G       | 30   | 0.10   | 0.09   | 0.09    | 0.07   |
| 6-311G**    | 72   | 0.64   | 0.49   | 0.49    | 0.41   |
| ANO         | 238  | 2527.6 | 51.13  | 53.59   | 37.78  |
| cc-pVDZ     | 58   | 0.45   | 0.34   | 0.24    | 0.21   |
| aug-cc-pVDZ | 100  | 1.87   | 1.18   | 1.02    | 0.85   |
| cc-pVTZ     | 144  | 4.98   | 4.82   | 2.65    | 2.05   |
| aug-cc-pVTZ | 230  | 26.03  | 23.12  | 12.40   | 9.27   |
| cc-pVQZ     | 290  | 81.24  | 65.12  | 31.51   | 22.60  |
| aug-cc-pVQZ | 436  | 444.23 | 324.29 | 151.04  | 107.24 |
