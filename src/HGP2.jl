# HGP2 - A (hopefully) fast, simple implementation of Head-Gordon, Pople's [^1]
# ERI recurrance relations.
#
# We're going to break this into several steps:
#
# 1. Primitive shell generation [ab,cd]
# 
# Given a shell-level description of a basis set ([ss,ss], [sp,ss], [dd,df], 
# including l [ll,ss]) generate all ERIs for the primitive basis functions in that shell.
# 
# 2. Contraction to (ab,cd)
# 
# 3. Integral Array Generation
# 
# Populate an integral record with the relevant terms for use in 
# an electronic structure theory code.
#
#
# 4. Future optimizations
# Gill's work on PRISM suggests [refs] being more flexible about when the basis function
# contraction is performed can reduce the operations count, but this will be simple and 
# likely fast.
# 
# There are many other recurrence relations to consider (MD [ref], LRL [ref], 
# OS [ref], Rys [refs]) but this should be a template for those others. In particular,
# several other schemes for the VRR have been proposed. For now, we'll just stick with the
# HGP version, since that is a well-written and reasonably self-contained paper.
# 
# Notation:
# Since I can't use [] or () in function names, I'm going to use a,b,c for primitive
# basis shells, and A,B,C for contracted shells. Therefore, we'll define routines
# like `ssss`, `psps`, `SSSS`, etc.
#
# 5. References
# [^1]: A method for two-electron Gaussian integral and integral derivative
#       evaluation using recurrence relations. Martin Head-Gordon and John
#       A. Pople. JCP, 89 (9), 5777, 1988.


# 0. Warm up
#    To see how well this works in practice, it might be useful to generate a few
#    simple cases and run them through steps 1-3 to see what's wrong with the plan.
#    Therefore, here are a few simple warm up exercises:
#
#    A. ssss and SSSS generation

"ssss - Calculate the [0]m terms for VRR equations for a range of m values"
function ssss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn,dxyz,ms)
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    qxyz = gaussian_product_center(cexpn,cxyz,dexpn,dxyz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    wxyz = gaussian_product_center(zeta,pxyz,eta,qxyz)
    rab2 = dist2(axyz-bxyz)
    rcd2 = dist2(cxyz-dxyz)
    rpq2 = dist2(pxyz-qxyz)
    T = zeta*eta/(zeta+eta)*rpq2
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    return Kab*Kcd/sqrt(zeta+eta)*[Fgamma(m,T) for m in ms]   # HGP eq 12
end
ssss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn,dxyz) = ssss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn,dxyz,[0])[1]

function psss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn,dxyz)
    sarray = ssss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn, dxyz,[0,1])
    @show sarray
    # Recalculate a number of terms from ssss. When the code works, be more
    # judicious in what we recalculate.
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    qxyz = gaussian_product_center(cexpn,cxyz,dexpn,dxyz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    wxyz = gaussian_product_center(zeta,pxyz,eta,qxyz)
    xsss = (pxyz[1]-axyz[1])*sarray[1] + (wxyz[1]-pxyz[1])*sarray[2]
    ysss = (pxyz[2]-axyz[2])*sarray[1] + (wxyz[2]-pxyz[2])*sarray[2]
    zsss = (pxyz[3]-axyz[3])*sarray[1] + (wxyz[3]-pxyz[3])*sarray[2]
    return [xsss,ysss,zsss]
end

#    B. pppp and PPPP generation
#    B'. Consider whether llll or LLLL is appropriate for sto-3g
#    C. Integral array generation for h2o/sto-3g, which should be do-able with the above.

# 1. Primitive shell generation [ab,cd]
#
#    A. [0]m generation
#    B. [a,c] generation (VRR)
#       Rewriting eq 6 from HGP with b=d=0 gives:
#           [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
#               + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)
#               + ci/2(zeta+eta)[a,c-1]m+1
#    C. [ab,cd] generation (HRR)
#

# 2. Contraction to (ab,cd)
# 
# 3. Integral Array Generation