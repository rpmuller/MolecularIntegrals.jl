# Notes on changes to MolecularIntegrals.jl
## 2021-05-22
Passed in pre-allocated space to vrr for some big savings when calling multiple times.

Interpolation of the Boys function looked like it would be a big win, but it was 50% slower 
in reality.

Maybe I spoke too soon. Getting much faster results now by defining the interpolation object
to be `const`. But I'm losing some accuracy in the resulting integrals, and I don't know why.
I thought that an earlier version didn't loose accuracy, but I don't know what changed.

## 2021-05-21
I'm going to push further Fgamma improvements until after the 0.1.0 release. Here's what I'd like to accomplish before then:
- Timing comparisons with Libint and Libcint
- Intelligent way of handling contraction coefficients

After the 0.1.0 release, I'm planning to post a message for help to the Julia discourse channel to speed up VRR further.
- Thus far, use of things like @inbounds and 

## 2021-05-08
Played around a little bit with the generated code dispatch table, and
even [asked a question on the Julia list](https://discourse.julialang.org/t/how-do-i-make-a-dispatch-table-using-multiple-dispatch-instead-of-dict/60784), 
but didn't succeed in making the code fast.

Working on a Pluto notebook for calculating `Fgamma` faster. 

Experimented with LoopVectorization.jl/@avx to speed up the inner loops in vrr.
It uniformly slowed everything down, which I think means that these just aren't
big enough. Removed.

Also tried to see whether @inbounds would speed things up. It didn't. Removed.

Guess I really need to make Fgamma faster. Moved the stock Fgamma code to using 
the asymptotic version when x is large. Sped things up maybe 10%.

## 2021-05-07
Streamlined the codegen a bit to reduce flops, but not a huge 
difference in the performance.

There is something wrong with the f-shell generated code. Seems to
infinite loop.

Actually, there isn't a problem with the f-shell, it just takes forever to compile, since those routines are >1kloc.

For now, removing the generated code, and trying to speed 
the rest.

Actually got the regular non-generated code to be much faster 
than the generated code. Much of this came from specifying
a type for ao2m[].

## 2021-05-06
Got the code generation working for vrr. Doesn't flag any bugs
that the vrr code didn't already flag, but it's certainly not
tested sufficiently. Code is in `tools/generate_vrr.jl` and
`HGP/vrr_autogen`.

## 2021-05-05
Finished the vrr_dd code, but it's untested and I assume has lots
of mistakes the code. It might be easier to write the autogen code
than to debug this code.

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
