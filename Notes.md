# Notes on changes to MolecularIntegrals.jl
## 2021-05-02
Attempted to pass in global arrays explicitly to local scope, i.e.:
```
function vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D
    # Pass in global arrays explicitly to local scope
    ,shift_index=shift_index,shift_direction=shift_direction
    )
```
Didn't have any impact on timings:
global ethane/6-31G/cc-pvdz=1.11/16.1, local 1.12,16.6

Can we build vrr around a new data structure? 
vrrs is a dict whose i,j element is the matrix of [iao,jao]. 
Maybe I need m's in there.

Some progress at implementing this, (HGP.jl/vrr_shell_dict) but it's 
not trivial. Also playing/thinking a bit with resurrecting the 
recursive version of the code, but to work with shell matrices
instead. The two endeavors are somewhat related, in that you're 
going to need to work with multiple shells at once, which requires
having access to them in some useful form.

One thing that could be much more efficient is the return of the
vrr and hrr routines. vrr returns a nao(ash+bsh) * nao(csh+dsh) matrix,
when it really just needs the nao(ash):nao(ash+bsh) * nao(csh):nao(csh+dsh)
quarter.

Similarly, hrr copies in all of these terms from vrr to hrr[:,1,:,1], but
it really only needs that last quarter of them, and then it just needs 
to return the nao(ash) * nao(bsh) * nao(csh) * nao(dsh) wedge.

Not only that, but when it's computing the b shell or the d shell, it computes
all of the corresponding a and c values.

I guess the right way to do this is to:
1. Return only the parts from hrr that we really need
2. Compute only these values in hrr
3. Return only the necessary pieces of vrr
4. Compute only the necessary pieces in vrr

One nice thing about the dictionary of matrices format is that it
makes it easier to return these pieces.

It's worth noting that instead of the dictionary of matrices, 
the full matrix has all the pieces one wants, and matrix views
make it easy to pass the memory behind a slice of the matrix
to the recursive routines, for example.

---------1---------2---------3---------4---------5---------6---------7---------8---------9
