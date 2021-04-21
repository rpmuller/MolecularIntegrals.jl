# TODO and notes

## New strategy for vrr5/hrr5:
[X] Write hrr5
[X] Benchmark vrr5/hrr5 against vrr/hrr
1. Optimize loop structure for vrr/hrr
2. Write special purpose routines for ss,sp,ps,pp,ds,sd,dp,pd,dd; only call general 
code when you exceed one of these cases, but make sure the general case still works
for all of these.
3. Create more structured testing of integral code using coulomb(), which should work now.
Consider something more like fuzz testing.
4. Delete old vrr/hrr implementations, since the vrr5/hrr5 scheme is the most efficient way to store integrals.
5. Write out a standard integral record using this code
6. Figure out metaprogramming to generate code for the special cases.
7. Decide on a better m2ao: are the arguments arrays or tuples. Are they both?

## chrr
- Finish chrr
- Figure out a way to share code between chrr and hrr

## Only computing, returning symmetric ERI pairs (i>j,k>l) ij>kl
Think about ordering an integral call (ij,kl) such that i>j, k>l, and ij>kl.
Also, put in warning when bsh>ash or dsh>csh in hrr2?

# Someday changes
## Make bf IJK and amplitude() calls to arrays

## Implement other molecule methods:
- nocc, nclosed, nopen, nup, ndown, stoich, center!

