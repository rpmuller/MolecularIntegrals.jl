# MolecularIntegrals.jl

The goal of MolecularIntegrals.jl is to supply **fast** and
**hackable** one- and two-electron integrals for computational
chemistry calculations.  There are many excellent molecular integral
packages available to Julia programmers, but few are written *in*
Julia. This project will explore how fast we can make these integrals
while maintaining a readable and hackable code base.

# One Electron Integrals

# Two Electron Integrals
MolecularIntegrals.jl supports a slower method based on Huzinaga's work,
and a faster set of integrals based on Head-Gordon and Pople's work [^HGP].

## Huzinaga Method

## Head-Gordan/Pople Recurrance Relations
```@docs
vrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
```

```@docs
hrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
```

## Other Methods
- Rys quadrature

# References
[^HGP]: A method for two-electron Gaussian integral and integral derivative
      evaluation using recurrence relations. Martin Head-Gordon and John
      A. Pople. JCP, 89 (9), 5777, 1988.
[^G]: Molecular Integrals Over Gaussian Basis Functions. Peter M. W. Gill. Adv.
      Q. Chem., 25, 141 (1994).
[^GP]: The Prism Algorithm for Two-Electron Integrals. Peter M. W. Gill and John
      A. Pople. IJQC, 40, 753 (1991).


