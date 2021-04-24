# MolecularIntegrals.jl

The goal of MolecularIntegrals.jl is to supply **fast** and
**hackable** one- and two-electron integrals for computational
chemistry calculations.  There are many excellent molecular integral
packages available to Julia programmers, but few are written *in*
Julia. This project will explore how fast we can make these integrals
while maintaining a readable and hackable code base.

MolecularIntegrals.jl strives to leverage the excellent work done by similar projects:

- [libints](https://github.com/evaleev/libint) and its Julia bindings [Lints.jl](https://github.com/FermiQC/Lints.jl)
- [Pyscf](https://github.com/pyscf/pyscf), the [libcint](https://github.com/sunqm/libcint) package.
- [JuliaChem.jl](https://github.com/davpoolechem/JuliaChem.jl)'s [JERI bindings](https://github.com/davpoolechem/JuliaChem.jl/tree/development/deps/src)
- [PyQuante](https://github.com/rpmuller/pyquante2)'s [python](https://github.com/rpmuller/pyquante2/tree/master/pyquante2/ints) and [c/cython](https://github.com/rpmuller/pyquante2/tree/master/cython) integrals, and the experimental [Julia version](https://github.com/rpmuller/pyquante2/tree/master/julia)



# One Electron Integrals
```@docs
overlap
kinetic
nuclear_attraction
```

# Two Electron Integrals
MolecularIntegrals.jl supports a slower method based on Huzinaga's work,
and a faster set of integrals based on Head-Gordon and Pople's work [^HGP].

The interface to the Huzinaga integrals is
```@docs
coulomb
```
These are useful for checking the values of other integrals, and have a 
simple, mature interface.

## Head-Gordon/Pople recurrence Relations
Integrals are computed using Head-Gordon and Pople's[^HGP] recurrence relations using
vertical (VRR) and horizontal (HRR) recurrence relations. In the notation of [^HGP], vertical recurrence relations 
construct integrals of the form [a0|c0] from kernels [00|00]^m via eq 6. Horizontal recurrence relations 
construct integrals of the form [ab|cd] from kernels [a0|c0].

### Vertical recurrence relations
MolecularIntegrals.jl returns VRRs either as a 2d array using `vrr_array`, or as 
a 6d array using `vrr_widearray`:
```@docs
vrr_array
vrr_widearray
```

### Horizontal recurrence relations
```@docs
hrr_array
hrr_dict
```


# References
[^HGP]: [A method for two-electron Gaussian integral and integral derivative
      evaluation using recurrence relations](https://doi.org/10.1063/1.455553). 
      Martin Head-Gordon and John A. Pople. JCP, 89 (9), 5777, 1988.
