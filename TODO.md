# TODO and notes

## New strategy for vrr/hrr:
1. Optimize loop structure for vrr/hrr
2. Time special purpose routines
3. Test special purpose routines

7. Decide on a better m2ao: are the arguments arrays or tuples. Are they both?

## chrr
- Test chrr

## Write out a standard integral record using this code


## Only computing, returning symmetric ERI pairs (i>j,k>l) ij>kl
Think about ordering an integral call (ij,kl) such that i>j, k>l, and ij>kl.
Also, put in warning when bsh>ash or dsh>csh in hrr2?

# Someday changes
## Make bf IJK and amplitude() calls to arrays
## Write/test special purpose code for VRRs: ds,sd,dp,pd,dd
## Create more structured testing of integral code using coulomb(), which should work now.
- Consider something more like fuzz testing.
## Figure out metaprogramming to generate code for the special hand-codedcases.

## Implement other molecule methods:
- nocc, nclosed, nopen, nup, ndown, stoich, center!

