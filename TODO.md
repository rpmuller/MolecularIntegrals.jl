# TODO and notes

## chrr
- Finish chrr
- Figure out a way to share code between chrr and hrr

## Figure out method to generate integral records
May have to generate an indexing array that points to the
different shells, primitive, and contracted functions.

## Reconcile coulomb(px,s,s,s) with hrr calls: 
While working on the new vrr code psss(), I found a discrepancy comparing to coulomb 
that I originally assumed was a mistake in psss(), but which I later found matched
vrr for this code. Which means that it is likely that the following test fails:
@test MolecularIntegrals.vrr(1.0,0,0,0,1,0,0,1.0,0,0,0,1.0,0,0,0,0,0,0,1.0,0,0,0,0) ≈ coulomb(px,s,s,s)
I'm going to move forward with coding the vrr routines, but I'm flagging this as
something to investigate and fix later.

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

