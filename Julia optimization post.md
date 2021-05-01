[tl;dr] I'd like to write a fast code to evaluate electron repulsion integrals in Julia, and would welcome input on how to optimize [the code I already have](https://github.com/rpmuller/MolecularIntegrals.jl/blob/master/src/HGP.jl) and how to use metaprogramming in a clever way to generate even faster code.

The fast codes currently in use generate code that is then compiled, using either C++ template metaprogramming [1] or using a common lisp script to write C code [2]. Codes written in modern languages like Python [3] are many times slower, and now typically just wrap one of the fast packages [1,2].

I think there's an opportunity to write a fast package in Julia that would also be readable and hackable. I've put together working code in the package MolecularIntegrals.jl [4], which is much faster than the pure python code. But it's still slower than Libint or Libcint.

The issue is the metaprogramming. These programs make use of recurrence relations [5] to generate integrals for higher angular momentum basis functions in terms of the lower angular momentum integrals, which are ultimate worked on in terms of incomplete error functions. But there are enough conditionals in the code that it's much faster to use metaprogramming techniques to write out a source file that can then be compiled into fast code.

Like I said, I have fairly decent Julia code in [5] (see the [HGP.jl](https://github.com/rpmuller/MolecularIntegrals.jl/blob/master/src/HGP.jl) file). I'd love to have advice on how to speed the code. And I'd also like some advice on how to use Julia's metaprogramming functionality to generate code in a readable way.

I've written some hand-optimized Julia code in the [HGPgen.jl](https://github.com/rpmuller/MolecularIntegrals.jl/blob/master/src/HGPgen.jl) file that doesn't have any loops or conditionals in it. This code is 10-20 times faster than the normal Julia code. I'm pretty confident I could write a Julia function that would generate a Julia file that I could then compile. But that seems like a waste of Julia's metaprogramming capabilities. Shouldn't there be a way to auto-generate these functions directly, rather than writing and compiling a file?

I'm just going to dive in and start writing routines to generate other routines, and hopefully something clever will suggest itself along the way. But I'd love to hear anyone's ideas about other ways to do this well, particularly if there are examples of other packages that have done so.

References
[1] [Libint: high-performance library for computing Gaussian integrals in quantum mechanics](https://github.com/evaleev/libint)
[2] [General GTO integrals for quantum chemistry](https://github.com/sunqm/libcint)
[3] [PyQuante](https://github.com/rpmuller/pyquante2)
[4] [MolecularIntegrals.jl](https://github.com/rpmuller/MolecularIntegrals.jl)
[5] [A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations](https://doi.org/10.1063/1.455553). Martin Head-Gordon and John A. Pople. JCP, 89 (9), 5777, 1988.