# TODO and notes

## New strategy for vrr5/hrr5:
X. Keep existing calling order, but just change (ax,ay,az) indices to integers via m2ao.
Get the code working and debugged for all existing tests. 
2. Write special purpose routines for ss,sp,ps,pp,ds,sd,dp,pd,dd; only call general 
code when you exceed one of these cases, but make sure the general case still works
for all of these.
3. Create more structured testing of integral code using coulomb(), which should work now.
Consider something more like fuzz testing.
4. Delete old vrr/hrr implementations, since the vrr5/hrr5 scheme is the most efficient way to store integrals.
5. Write out a standard integral record using this code
6. Figure out metaprogramming to generate code for the special cases.
7. At some point (maybe earlier than step 7), move vrr to static arrays.
8. Decide on a better m2ao: are the arguments arrays or tuples. Are they both?
X. Fix problems with mmax: do I need to increment this by one?

## chrr
- Finish chrr
- Figure out a way to share code between chrr and hrr

## Figure out method to generate integral records
May have to generate an indexing array that points to the
different shells, primitive, and contracted functions.

## Only computing, returning symmetric ERI pairs (i>j,k>l) ij>kl
Think about ordering an integral call (ij,kl) such that i>j, k>l, and ij>kl.
Also, put in warning when bsh>ash or dsh>csh in hrr2?


# Someday changes
## Move basis functions from bf.x,bf.y,bf.z to bf.xyz
- [X] Do CGBFs 
- Then amplitude(CGBF)
- [X] then PGBFs
- Then amplitude(CGBF)
- [X] Scan for .y, .z
- Also consider a similar move for bf.I, bf.J, bf.K.

## Implement other molecule methods:
- nocc, nclosed, nopen, nup, ndown, stoich, center!

