# TODO and notes
- [ ] Reconcile coulomb(px,s,s,s) with psss(): see note below
- [ ] Note on the interface: currently the new vrr routines are returning multidimensional 
    arrays over angular momentum states and m. If mmax=0 is submitted, this will reduce the return value by 1,
    allowing, for example, a 3x1 maxtrix to be accessed as a 3 vector. Could cause problems with higher dimension
    returns (3x3x1, e.g.)
- [ ] The ultimate new vrr code will return a [ix,iy,iz,jx,jy,jz,m] matrix. Each term of this will be a primitive 
    integral. We will need to find a way to do the contraction over multiple primitives, potentially introducing another index that we'll contract over: [k,ix,iy,iz,jx,jy,jz,m]. Although an 8-dimensional array is large, none of the dimensions will have very many terms. Still 4^8 is a lot of terms. I wonder whether this will require a more intelligent data structure.
    - [ ] Could also use a SparseArray
    - [ ] Such a data structure could be using the ao l,m indices to map the I,J,K. There's preliminary code in Basis.jl 
    for these transformations.
- [ ] Think about ordering an integral call (ij,kl) such that i>j, k>l, and ij>kl.
    - [ ] Put in warning when bsh>ash or dsh>csh in hrr2?


 # TODO: reconcile coulomb(px,s,s,s) with psss():
 #   While working on the new vrr code psss(), I found a discrepancy comparing to coulomb 
 #   that I originally assumed was a mistake in psss(), but which I later found matched
 #   vrr for this code. Which means that it is likely that the following test fails:
 #   @test MolecularIntegrals.vrr(1.0,0,0,0,1,0,0,1.0,0,0,0,1.0,0,0,0,0,0,0,1.0,0,0,0,0) ≈ coulomb(px,s,s,s)
 #   I'm going to move forward with coding the vrr routines, but I'm flagging this as
 #   something to investigate and fix later.
