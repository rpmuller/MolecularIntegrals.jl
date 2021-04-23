# TODO and notes

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

