# TODO and notes

## separate VRR.jl and HRR.jl code

## Delete HGPother.jl

## Rename HGP2.jl to HGPold.jl

## Write vrr_array and vrr_wide_array versions
- These are vrr5 and vrr1, respectively
- Comment the code to indicate that we're supporting multiple
    interfaces for convenience, and that that the speed is 
    roughly equivalent

## Write hrr_dict and hrr_array versions
- These are the old version of hrr1 and hrr5, respectively
- Should be able to call *either* vrr_array or vrr_wide
- Comment code to indicate that we're supporting multiple
    interfaces for convenience, and that the hrr_array is
    faster, but the hrr_dict can only return the requisite
    terms.
- Prune out unnecessary results from hrr_dict (when != ashell, != bshell, etc.)



## Write test functions in runtests.jl
- vrr_test(ashell,cshell) code:
    - Test all elements of vrr_array and vrr_wide_array against ERI.jl code
    - Random A,C, aex, cex?
- hrr_test(ashell,bshell,cshell,dshell)
    - Test all elements of hrr_array and hrr_dict against ERI.jl code
    - Random A,B,C,D???





## Comment ERI.jl code
- Add comments to ERI code to the effect that it is slow code that is very
    likely to be correct that is kept around for reference and testing.

## Move contracted routines to CHGP.jl
- Test chrr
- Write out a standard integral record using this code

# Someday changes
- Make bf I,J,K and amplitude() calls arrays [I,J,K]
- Decide on a better m2ao: are the arguments arrays or tuples. Are they both?

## Write/test special purpose code 
- Expand the m-loops in sp routines
- Time special purpose sp routines 
    Special purpose code is 20x faster:
    vrr5(1,1)   33.262 μs (384 allocations: 15.94 KiB)
    vrr_pp      1.586 μs (16 allocations: 2.66 KiB)

- Test special purpose sp routines
- Write/time/test ds,sd,dp,pd,dd
- Figure out metaprogramming to generate code for the special hand-codedcases.

## Create more structured testing of integral code using coulomb(), which should work now.
- Consider something more like fuzz testing.


## Implement other molecule methods:
- nocc, nclosed, nopen, nup, ndown, stoich, center!

