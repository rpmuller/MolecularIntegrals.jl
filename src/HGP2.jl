"""
vrr_ss(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A and C
shells both have s-type angular momentum. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_ss(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    vrrs = zeros(Float64,1,1)

    P = gaussian_product_center(aexpn,A,bexpn,B)
    Q = gaussian_product_center(cexpn,C,dexpn,D)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    W = gaussian_product_center(zeta,P,eta,Q)
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2/ze
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    
    vrrs[1,1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    return vrrs
end

"""
vrr_ps(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A shell has
p-type angular momentum and the C shell has s-type a.m. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_ps(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 2
    vrrs = zeros(Float64,4,1,mmax)

    P = gaussian_product_center(aexpn,A,bexpn,B)
    Q = gaussian_product_center(cexpn,C,dexpn,D)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    W = gaussian_product_center(zeta,P,eta,Q)
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2/ze
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)

    vrrs[2,1,1] = (P[1]-A[1])*vrrs[1,1,1] + (W[1]-A[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (P[2]-A[2])*vrrs[1,1,1] + (W[2]-A[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (P[3]-A[3])*vrrs[1,1,1] + (W[3]-A[3])*vrrs[1,1,2]

    return vrrs[:,:,1]
end

"""
vrr_sp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A shell has
s-type angular momentum and the C shell has p-type a.m. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_sp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 2
    vrrs = zeros(Float64,1,4,mmax)

    P = gaussian_product_center(aexpn,A,bexpn,B)
    Q = gaussian_product_center(cexpn,C,dexpn,D)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    W = gaussian_product_center(zeta,P,eta,Q)
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2/ze
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)

    vrrs[1,2,1] = (Q[1]-C[1])*vrrs[1,1,1] + (W[1]-Q[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (Q[2]-C[2])*vrrs[1,1,1] + (W[2]-Q[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (Q[3]-C[3])*vrrs[1,1,1] + (W[3]-Q[3])*vrrs[1,1,2]

    return vrrs[:,:,1]
end

"""
vrr_pp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A and C
shells both have p-type angular momentum. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_pp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 3
    vrrs = zeros(Float64,4,4,mmax)

    P = gaussian_product_center(aexpn,A,bexpn,B)
    Q = gaussian_product_center(cexpn,C,dexpn,D)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    W = gaussian_product_center(zeta,P,eta,Q)
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2/ze
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)
    vrrs[1,1, 3] = Kab*Kcd*Fgamma(2,T)/sqrt(ze)

    vrrs[2,1,1] = (P[1]-A[1])*vrrs[1,1,1] + (W[1]-A[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (P[2]-A[2])*vrrs[1,1,1] + (W[2]-A[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (P[3]-A[3])*vrrs[1,1,1] + (W[3]-A[3])*vrrs[1,1,2]
    vrrs[1,2,1] = (Q[1]-C[1])*vrrs[1,1,1] + (W[1]-Q[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (Q[2]-C[2])*vrrs[1,1,1] + (W[2]-Q[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (Q[3]-C[3])*vrrs[1,1,1] + (W[3]-Q[3])*vrrs[1,1,2]
    vrrs[2,1,2] = (P[1]-A[1])*vrrs[1,1,2] + (W[1]-A[1])*vrrs[1,1,3]
    vrrs[3,1,2] = (P[2]-A[2])*vrrs[1,1,2] + (W[2]-A[2])*vrrs[1,1,3]
    vrrs[4,1,2] = (P[3]-A[3])*vrrs[1,1,2] + (W[3]-A[3])*vrrs[1,1,3]
    vrrs[1,2,2] = (Q[1]-C[1])*vrrs[1,1,2] + (W[1]-Q[1])*vrrs[1,1,3]
    vrrs[1,3,2] = (Q[2]-C[2])*vrrs[1,1,2] + (W[2]-Q[2])*vrrs[1,1,3]
    vrrs[1,4,2] = (Q[3]-C[3])*vrrs[1,1,2] + (W[3]-Q[3])*vrrs[1,1,3]

    vrrs[2,2,1] = (Q[1]-C[1])*vrrs[2,1,1] + (W[1]-Q[1])*vrrs[2,1,2] +
        0.5/ze*vrrs[1,1,2]
    vrrs[2,3,1] = (Q[2]-C[2])*vrrs[2,1,1] + (W[2]-Q[2])*vrrs[2,1,2] +
        0.5/ze*vrrs[1,1,2]
    vrrs[2,4,1] = (Q[3]-C[3])*vrrs[2,1,1] + (W[3]-Q[3])*vrrs[2,1,2] +
        0.5/ze*vrrs[1,1,2]

    vrrs[3,2,1] = (Q[1]-C[1])*vrrs[3,1,1] + (W[1]-Q[1])*vrrs[3,1,2] +
        0.5/ze*vrrs[1,1,2]
    vrrs[3,3,1] = (Q[2]-C[2])*vrrs[3,1,1] + (W[2]-Q[2])*vrrs[3,1,2] +
        0.5/ze*vrrs[1,1,2]
    vrrs[3,4,1] = (Q[3]-C[3])*vrrs[3,1,1] + (W[3]-Q[3])*vrrs[3,1,2] +
        0.5/ze*vrrs[1,1,2]

    vrrs[4,2,1] = (Q[1]-C[1])*vrrs[4,1,1] + (W[1]-Q[1])*vrrs[4,1,2] +
        0.5/ze*vrrs[1,1,2]
    vrrs[4,3,1] = (Q[2]-C[2])*vrrs[4,1,1] + (W[2]-Q[2])*vrrs[4,1,2] +
        0.5/ze*vrrs[1,1,2]
    vrrs[4,4,1] = (Q[3]-C[3])*vrrs[4,1,1] + (W[3]-Q[3])*vrrs[4,1,2] +
        0.5/ze*vrrs[1,1,2]
    return vrrs[:,:,1]
end