#HGPgen - code for hand-written and auto-generated routines

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

    PA = P-A
    WP = W-P
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)

    vrrs[2,1,1] = (PA[1])*vrrs[1,1,1] + (WP[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (PA[2])*vrrs[1,1,1] + (WP[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (PA[3])*vrrs[1,1,1] + (WP[3])*vrrs[1,1,2]

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

    QC = Q-C
    WQ = W-Q
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)

    vrrs[1,2,1] = (QC[1])*vrrs[1,1,1] + (WQ[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (QC[2])*vrrs[1,1,1] + (WQ[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (QC[3])*vrrs[1,1,1] + (WQ[3])*vrrs[1,1,2]

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
    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)
    vrrs[1,1, 3] = Kab*Kcd*Fgamma(2,T)/sqrt(ze)

    vrrs[2,1,1] = (PA[1])*vrrs[1,1,1] + (WP[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (PA[2])*vrrs[1,1,1] + (WP[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (PA[3])*vrrs[1,1,1] + (WP[3])*vrrs[1,1,2]
    vrrs[1,2,1] = (QC[1])*vrrs[1,1,1] + (WQ[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (QC[2])*vrrs[1,1,1] + (WQ[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (QC[3])*vrrs[1,1,1] + (WQ[3])*vrrs[1,1,2]
    vrrs[2,1,2] = (PA[1])*vrrs[1,1,2] + (WP[1])*vrrs[1,1,3]
    vrrs[3,1,2] = (PA[2])*vrrs[1,1,2] + (WP[2])*vrrs[1,1,3]
    vrrs[4,1,2] = (PA[3])*vrrs[1,1,2] + (WP[3])*vrrs[1,1,3]
    vrrs[1,2,2] = (QC[1])*vrrs[1,1,2] + (WQ[1])*vrrs[1,1,3]
    vrrs[1,3,2] = (QC[2])*vrrs[1,1,2] + (WQ[2])*vrrs[1,1,3]
    vrrs[1,4,2] = (QC[3])*vrrs[1,1,2] + (WQ[3])*vrrs[1,1,3]

    vrrs[2,2,1] = (QC[1])*vrrs[2,1,1] + (WQ[1])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,3,1] = (QC[2])*vrrs[2,1,1] + (WQ[2])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,4,1] = (QC[3])*vrrs[2,1,1] + (WQ[3])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[3,2,1] = (QC[1])*vrrs[3,1,1] + (WQ[1])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,3,1] = (QC[2])*vrrs[3,1,1] + (WQ[2])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,4,1] = (QC[3])*vrrs[3,1,1] + (WQ[3])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[4,2,1] = (QC[1])*vrrs[4,1,1] + (WQ[1])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,3,1] = (QC[2])*vrrs[4,1,1] + (WQ[2])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,4,1] = (QC[3])*vrrs[4,1,1] + (WQ[3])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    return vrrs[:,:,1]
end

"""
vrr_sd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A shell has
s-type angular momentum and the C shell has d-type a.m. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_sd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 3
    vrrs = zeros(Float64,1,10,mmax)

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

    QC = Q-C
    WQ = W-Q
    
    vrrs[1,2,1] = (QC[1])*vrrs[1,1,1] + (WQ[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (QC[2])*vrrs[1,1,1] + (WQ[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (QC[3])*vrrs[1,1,1] + (WQ[3])*vrrs[1,1,2]
    
    vrrs[1,2,2] = (QC[1])*vrrs[1,1,2] + (WQ[1])*vrrs[1,1,3]
    vrrs[1,3,2] = (QC[2])*vrrs[1,1,2] + (WQ[2])*vrrs[1,1,3]
    vrrs[1,4,2] = (QC[3])*vrrs[1,1,2] + (WQ[3])*vrrs[1,1,3]

    # [MVector(2,0,0),MVector(1,1,0),MVector(1,0,1),MVector(0,2,0),MVector(0,1,1),MVector(0,0,2)],

    vrrs[1,5,1] = (QC[1])*vrrs[1,2,1] + (WQ[1])*vrrs[1,2,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,6,1] = (QC[2])*vrrs[1,2,1] + (WQ[2])*vrrs[1,2,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,7,1] = (QC[3])*vrrs[1,2,1] + (WQ[3])*vrrs[1,2,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,8,1] = (QC[2])*vrrs[1,3,1] + (WQ[2])*vrrs[1,3,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,9,1] = (QC[3])*vrrs[1,3,1] + (WQ[3])*vrrs[1,3,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,10,1] =(QC[3])*vrrs[1,4,1] + (WQ[3])*vrrs[1,4,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])

    return vrrs[:,:,1]
end

"""
vrr_ds(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A shell has
c-type angular momentum and the C shell has s-type a.m. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_ds(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 3
    vrrs = zeros(Float64,10,1,mmax)

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
    PA = P-A
    WP = W-P
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)
    vrrs[1,1, 3] = Kab*Kcd*Fgamma(2,T)/sqrt(ze)

    vrrs[2,1,1] = (PA[1])*vrrs[1,1,1] + (WP[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (PA[2])*vrrs[1,1,1] + (WP[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (PA[3])*vrrs[1,1,1] + (WP[3])*vrrs[1,1,2]

    vrrs[2,1,2] = (PA[1])*vrrs[1,1,2] + (WP[1])*vrrs[1,1,3]
    vrrs[3,1,2] = (PA[2])*vrrs[1,1,2] + (WP[2])*vrrs[1,1,3]
    vrrs[4,1,2] = (PA[3])*vrrs[1,1,2] + (WP[3])*vrrs[1,1,3]

    # [MVector(2,0,0),MVector(1,1,0),MVector(1,0,1),MVector(0,2,0),MVector(0,1,1),MVector(0,0,2)],

    vrrs[5,1,1] = (PA[1])*vrrs[2,1,1] + (WP[1])*vrrs[2,1,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[6,1,1] = (PA[2])*vrrs[2,1,1] + (WP[2])*vrrs[2,1,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[7,1,1] = (PA[3])*vrrs[2,1,1] + (WP[3])*vrrs[2,1,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[8,1,1] = (PA[2])*vrrs[3,1,1] + (WP[2])*vrrs[3,1,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[9,1,1] = (PA[3])*vrrs[3,1,1] + (WP[3])*vrrs[3,1,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[10,1,1] =(PA[3])*vrrs[4,1,1] + (WP[3])*vrrs[4,1,2] +  1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])

    return vrrs[:,:,1]
end

"""
vrr_dp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A shell has
c-type angular momentum and the C shell has s-type a.m. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_dp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 4
    vrrs = zeros(Float64,10,4,mmax)

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

    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)
    vrrs[1,1, 3] = Kab*Kcd*Fgamma(2,T)/sqrt(ze)
    vrrs[1,1, 4] = Kab*Kcd*Fgamma(3,T)/sqrt(ze)

    vrrs[2,1,1] = (PA[1])*vrrs[1,1,1] + (WP[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (PA[2])*vrrs[1,1,1] + (WP[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (PA[3])*vrrs[1,1,1] + (WP[3])*vrrs[1,1,2]

    vrrs[2,1,2] = (PA[1])*vrrs[1,1,2] + (WP[1])*vrrs[1,1,3]
    vrrs[3,1,2] = (PA[2])*vrrs[1,1,2] + (WP[2])*vrrs[1,1,3]
    vrrs[4,1,2] = (PA[3])*vrrs[1,1,2] + (WP[3])*vrrs[1,1,3]

    vrrs[1,2,1] = (QC[1])*vrrs[1,1,1] + (WQ[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (QC[2])*vrrs[1,1,1] + (WQ[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (QC[3])*vrrs[1,1,1] + (WQ[3])*vrrs[1,1,2]
    
    vrrs[1,2,2] = (QC[1])*vrrs[1,1,2] + (WQ[1])*vrrs[1,1,3]
    vrrs[1,3,2] = (QC[2])*vrrs[1,1,2] + (WQ[2])*vrrs[1,1,3]
    vrrs[1,4,2] = (QC[3])*vrrs[1,1,2] + (WQ[3])*vrrs[1,1,3]

    vrrs[2,2,1] = (QC[1])*vrrs[2,1,1] + (WQ[1])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,3,1] = (QC[2])*vrrs[2,1,1] + (WQ[2])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,4,1] = (QC[3])*vrrs[2,1,1] + (WQ[3])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[3,2,1] = (QC[1])*vrrs[3,1,1] + (WQ[1])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,3,1] = (QC[2])*vrrs[3,1,1] + (WQ[2])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,4,1] = (QC[3])*vrrs[3,1,1] + (WQ[3])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[4,2,1] = (QC[1])*vrrs[4,1,1] + (WQ[1])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,3,1] = (QC[2])*vrrs[4,1,1] + (WQ[2])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,4,1] = (QC[3])*vrrs[4,1,1] + (WQ[3])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
 
    # [MVector(2,0,0),MVector(1,1,0),MVector(1,0,1),MVector(0,2,0),MVector(0,1,1),MVector(0,0,2)],

    vrrs[5,1,1] = (PA[1])*vrrs[2,1,1] + (WP[1])*vrrs[2,1,2] +  
            1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[6,1,1] = (PA[2])*vrrs[2,1,1] + (WP[2])*vrrs[2,1,2] +  
            1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[7,1,1] = (PA[3])*vrrs[2,1,1] + (WP[3])*vrrs[2,1,2] +  
            1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[8,1,1] = (PA[2])*vrrs[3,1,1] + (WP[2])*vrrs[3,1,2] +  
            1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[9,1,1] = (PA[3])*vrrs[3,1,1] + (WP[3])*vrrs[3,1,2] +  
            1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[10,1,1] =(PA[3])*vrrs[4,1,1] + (WP[3])*vrrs[4,1,2] +  
            1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])

    vrrs[5,2,1] = (PA[1])*vrrs[2,2,1] + (WP[1])*vrrs[2,2,2] +  
            1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[6,2,1] = (PA[2])*vrrs[2,2,1] + (WP[2])*vrrs[2,2,2] +  
            1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[7,2,1] = (PA[3])*vrrs[2,2,1] + (WP[3])*vrrs[2,2,2] +  
            1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[8,2,1] = (PA[2])*vrrs[3,2,1] + (WP[2])*vrrs[3,2,2] +  
            1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[9,2,1] = (PA[3])*vrrs[3,2,1] + (WP[3])*vrrs[3,2,2] +  
            1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[10,2,1] =(PA[3])*vrrs[4,2,1] + (WP[3])*vrrs[4,2,2] +  
            1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]

    vrrs[5,3,1] = (PA[1])*vrrs[2,3,1] + (WP[1])*vrrs[2,3,2] +  
            1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[6,3,1] = (PA[2])*vrrs[2,3,1] + (WP[2])*vrrs[2,3,2] +  
            1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[7,3,1] = (PA[3])*vrrs[2,3,1] + (WP[3])*vrrs[2,3,2] +  
            1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[8,3,1] = (PA[2])*vrrs[3,3,1] + (WP[2])*vrrs[3,3,2] +  
            1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[9,3,1] = (PA[3])*vrrs[3,3,1] + (WP[3])*vrrs[3,3,2] +  
            1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[10,3,1] =(PA[3])*vrrs[4,3,1] + (WP[3])*vrrs[4,3,2] +  
            1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]

    vrrs[5,4,1] = (PA[1])*vrrs[2,4,1] + (WP[1])*vrrs[2,4,2] +  
            1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[6,4,1] = (PA[2])*vrrs[2,4,1] + (WP[2])*vrrs[2,4,2] +  
            1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[7,4,1] = (PA[3])*vrrs[2,4,1] + (WP[3])*vrrs[2,4,2] +  
            1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[8,4,1] = (PA[2])*vrrs[3,4,1] + (WP[2])*vrrs[3,4,2] +  
            1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[9,4,1] = (PA[3])*vrrs[3,4,1] + (WP[3])*vrrs[3,4,2] +  
            1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[10,4,1] =(PA[3])*vrrs[4,4,1] + (WP[3])*vrrs[4,4,2] +  
            1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]

    return vrrs[:,:,1]
end

"""
vrr_pd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A shell has
c-type angular momentum and the C shell has s-type a.m. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_pd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 4
    vrrs = zeros(Float64,4,10,mmax)

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
    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)
    vrrs[1,1, 3] = Kab*Kcd*Fgamma(2,T)/sqrt(ze)
    vrrs[1,1, 4] = Kab*Kcd*Fgamma(3,T)/sqrt(ze)

    vrrs[2,1,1] = (PA[1])*vrrs[1,1,1] + (WP[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (PA[2])*vrrs[1,1,1] + (WP[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (PA[3])*vrrs[1,1,1] + (WP[3])*vrrs[1,1,2]

    vrrs[2,1,2] = (PA[1])*vrrs[1,1,2] + (WP[1])*vrrs[1,1,3]
    vrrs[3,1,2] = (PA[2])*vrrs[1,1,2] + (WP[2])*vrrs[1,1,3]
    vrrs[4,1,2] = (PA[3])*vrrs[1,1,2] + (WP[3])*vrrs[1,1,3]

    vrrs[1,2,1] = (QC[1])*vrrs[1,1,1] + (WQ[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (QC[2])*vrrs[1,1,1] + (WQ[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (QC[3])*vrrs[1,1,1] + (WQ[3])*vrrs[1,1,2]
    
    vrrs[1,2,2] = (QC[1])*vrrs[1,1,2] + (WQ[1])*vrrs[1,1,3]
    vrrs[1,3,2] = (QC[2])*vrrs[1,1,2] + (WQ[2])*vrrs[1,1,3]
    vrrs[1,4,2] = (QC[3])*vrrs[1,1,2] + (WQ[3])*vrrs[1,1,3]

    vrrs[2,2,1] = (QC[1])*vrrs[2,1,1] + (WQ[1])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,3,1] = (QC[2])*vrrs[2,1,1] + (WQ[2])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,4,1] = (QC[3])*vrrs[2,1,1] + (WQ[3])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[3,2,1] = (QC[1])*vrrs[3,1,1] + (WQ[1])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,3,1] = (QC[2])*vrrs[3,1,1] + (WQ[2])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,4,1] = (QC[3])*vrrs[3,1,1] + (WQ[3])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[4,2,1] = (QC[1])*vrrs[4,1,1] + (WQ[1])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,3,1] = (QC[2])*vrrs[4,1,1] + (WQ[2])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,4,1] = (QC[3])*vrrs[4,1,1] + (WQ[3])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
 
    # [MVector(2,0,0),MVector(1,1,0),MVector(1,0,1),MVector(0,2,0),MVector(0,1,1),MVector(0,0,2)],
    vrrs[1,5,1] = (QC[1])*vrrs[1,2,1] + (WQ[1])*vrrs[1,2,2] +
          1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,6,1] = (QC[2])*vrrs[1,2,1] + (WQ[2])*vrrs[1,2,2] +
          1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,7,1] = (QC[3])*vrrs[1,2,1] + (WQ[3])*vrrs[1,2,2] +
          1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,8,1] = (QC[2])*vrrs[1,3,1] + (WQ[2])*vrrs[1,3,2] +
          1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,9,1] = (QC[3])*vrrs[1,3,1] + (WQ[3])*vrrs[1,3,2] +
          1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,10,1] =(QC[3])*vrrs[1,4,1] + (WQ[3])*vrrs[1,4,2] +
         1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])


    vrrs[2,5,1] = (QC[1])*vrrs[2,2,1] + (WQ[1])*vrrs[2,2,2] +
          1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,6,1] = (QC[2])*vrrs[2,2,1] + (WQ[2])*vrrs[2,2,2] +
          1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,7,1] = (QC[3])*vrrs[2,2,1] + (WQ[3])*vrrs[2,2,2] +
          1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,8,1] = (QC[2])*vrrs[2,3,1] + (WQ[2])*vrrs[2,3,2] +
          1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,9,1] = (QC[3])*vrrs[2,3,1] + (WQ[3])*vrrs[2,3,2] +
          1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,10,1] =(QC[3])*vrrs[2,4,1] + (WQ[3])*vrrs[2,4,2] +
         1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]


    vrrs[3,5,1] = (QC[1])*vrrs[3,2,1] + (WQ[1])*vrrs[3,2,2] +
          1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,6,1] = (QC[2])*vrrs[3,2,1] + (WQ[2])*vrrs[3,2,2] +
          1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,7,1] = (QC[3])*vrrs[3,2,1] + (WQ[3])*vrrs[3,2,2] +
          1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,8,1] = (QC[2])*vrrs[3,3,1] + (WQ[2])*vrrs[3,3,2] +
          1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,9,1] = (QC[3])*vrrs[3,3,1] + (WQ[3])*vrrs[3,3,2] +
          1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,10,1] =(QC[3])*vrrs[3,4,1] + (WQ[3])*vrrs[3,4,2] +
         1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]


    vrrs[4,5,1] = (QC[1])*vrrs[4,2,1] + (WQ[1])*vrrs[4,2,2] +
          1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,6,1] = (QC[2])*vrrs[4,2,1] + (WQ[2])*vrrs[4,2,2] +
          1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,7,1] = (QC[3])*vrrs[4,2,1] + (WQ[3])*vrrs[4,2,2] +
          1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,8,1] = (QC[2])*vrrs[4,3,1] + (WQ[2])*vrrs[4,3,2] +
          1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,9,1] = (QC[3])*vrrs[4,3,1] + (WQ[3])*vrrs[4,3,2] +
          1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,10,1] =(QC[3])*vrrs[4,4,1] + (WQ[3])*vrrs[4,4,2] +
         1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]

    return vrrs[:,:,1]
end

"""
vrr_dd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.

This is a hand-written routine specific to when the A shell has
d-type angular momentum and the C shell has d-type a.m. It is meant to be
a model for machine generated angular momentum specific code.
"""
function vrr_dd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = 5
    vrrs = zeros(Float64,10,10,mmax)

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

    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q
    
    vrrs[1,1, 1] = Kab*Kcd*Fgamma(0,T)/sqrt(ze)
    vrrs[1,1, 2] = Kab*Kcd*Fgamma(1,T)/sqrt(ze)
    vrrs[1,1, 3] = Kab*Kcd*Fgamma(2,T)/sqrt(ze)
    vrrs[1,1, 4] = Kab*Kcd*Fgamma(3,T)/sqrt(ze)
    vrrs[1,1, 5] = Kab*Kcd*Fgamma(4,T)/sqrt(ze)

    vrrs[2,1,1] = (PA[1])*vrrs[1,1,1] + (WP[1])*vrrs[1,1,2]
    vrrs[3,1,1] = (PA[2])*vrrs[1,1,1] + (WP[2])*vrrs[1,1,2]
    vrrs[4,1,1] = (PA[3])*vrrs[1,1,1] + (WP[3])*vrrs[1,1,2]

    vrrs[2,1,2] = (PA[1])*vrrs[1,1,2] + (WP[1])*vrrs[1,1,3]
    vrrs[3,1,2] = (PA[2])*vrrs[1,1,2] + (WP[2])*vrrs[1,1,3]
    vrrs[4,1,2] = (PA[3])*vrrs[1,1,2] + (WP[3])*vrrs[1,1,3]

    vrrs[1,2,1] = (QC[1])*vrrs[1,1,1] + (WQ[1])*vrrs[1,1,2]
    vrrs[1,3,1] = (QC[2])*vrrs[1,1,1] + (WQ[2])*vrrs[1,1,2]
    vrrs[1,4,1] = (QC[3])*vrrs[1,1,1] + (WQ[3])*vrrs[1,1,2]
    
    vrrs[1,2,2] = (QC[1])*vrrs[1,1,2] + (WQ[1])*vrrs[1,1,3]
    vrrs[1,3,2] = (QC[2])*vrrs[1,1,2] + (WQ[2])*vrrs[1,1,3]
    vrrs[1,4,2] = (QC[3])*vrrs[1,1,2] + (WQ[3])*vrrs[1,1,3]

    vrrs[2,2,1] = (QC[1])*vrrs[2,1,1] + (WQ[1])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,3,1] = (QC[2])*vrrs[2,1,1] + (WQ[2])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[2,4,1] = (QC[3])*vrrs[2,1,1] + (WQ[3])*vrrs[2,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[3,2,1] = (QC[1])*vrrs[3,1,1] + (WQ[1])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,3,1] = (QC[2])*vrrs[3,1,1] + (WQ[2])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[3,4,1] = (QC[3])*vrrs[3,1,1] + (WQ[3])*vrrs[3,1,2] + 0.5/ze*vrrs[1,1,2]

    vrrs[4,2,1] = (QC[1])*vrrs[4,1,1] + (WQ[1])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,3,1] = (QC[2])*vrrs[4,1,1] + (WQ[2])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
    vrrs[4,4,1] = (QC[3])*vrrs[4,1,1] + (WQ[3])*vrrs[4,1,2] + 0.5/ze*vrrs[1,1,2]
 
    # [MVector(2,0,0),MVector(1,1,0),MVector(1,0,1),MVector(0,2,0),MVector(0,1,1),MVector(0,0,2)],
    vrrs[1,5,1] = (QC[1])*vrrs[1,2,1] + (WQ[1])*vrrs[1,2,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,6,1] = (QC[2])*vrrs[1,2,1] + (WQ[2])*vrrs[1,2,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,7,1] = (QC[3])*vrrs[1,2,1] + (WQ[3])*vrrs[1,2,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,8,1] = (QC[2])*vrrs[1,3,1] + (WQ[2])*vrrs[1,3,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,9,1] = (QC[3])*vrrs[1,3,1] + (WQ[3])*vrrs[1,3,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[1,10,1] =(QC[3])*vrrs[1,4,1] + (WQ[3])*vrrs[1,4,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])

    vrrs[2,5,1] = (QC[1])*vrrs[2,2,1] + (WQ[1])*vrrs[2,2,2] +
              1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,6,1] = (QC[2])*vrrs[2,2,1] + (WQ[2])*vrrs[2,2,2] +
              1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,7,1] = (QC[3])*vrrs[2,2,1] + (WQ[3])*vrrs[2,2,2] +
              1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,8,1] = (QC[2])*vrrs[2,3,1] + (WQ[2])*vrrs[2,3,2] +
              1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,9,1] = (QC[3])*vrrs[2,3,1] + (WQ[3])*vrrs[2,3,2] +
              1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[2,10,1] =(QC[3])*vrrs[2,4,1] + (WQ[3])*vrrs[2,4,2] +
              1/(2*zeta)*(vrrs[2,1,1]-eta/ze*vrrs[2,1,2]) + 0.5/ze*vrrs[1,1,2]

    vrrs[3,5,1] = (QC[1])*vrrs[3,2,1] + (WQ[1])*vrrs[3,2,2] +
              1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,6,1] = (QC[2])*vrrs[3,2,1] + (WQ[2])*vrrs[3,2,2] +
              1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,7,1] = (QC[3])*vrrs[3,2,1] + (WQ[3])*vrrs[3,2,2] +
              1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,8,1] = (QC[2])*vrrs[3,3,1] + (WQ[2])*vrrs[3,3,2] +
              1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,9,1] = (QC[3])*vrrs[3,3,1] + (WQ[3])*vrrs[3,3,2] +
              1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[3,10,1] =(QC[3])*vrrs[3,4,1] + (WQ[3])*vrrs[3,4,2] +
              1/(2*zeta)*(vrrs[3,1,1]-eta/ze*vrrs[3,1,2]) + 0.5/ze*vrrs[1,1,2]

    vrrs[4,5,1] = (QC[1])*vrrs[4,2,1] + (WQ[1])*vrrs[4,2,2] +
              1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,6,1] = (QC[2])*vrrs[4,2,1] + (WQ[2])*vrrs[4,2,2] +
              1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,7,1] = (QC[3])*vrrs[4,2,1] + (WQ[3])*vrrs[4,2,2] +
              1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,8,1] = (QC[2])*vrrs[4,3,1] + (WQ[2])*vrrs[4,3,2] +
              1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,9,1] = (QC[3])*vrrs[4,3,1] + (WQ[3])*vrrs[4,3,2] +
              1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[4,10,1] =(QC[3])*vrrs[4,4,1] + (WQ[3])*vrrs[4,4,2] +
              1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]
 
    vrrs[5,1,1] = (PA[1])*vrrs[2,1,1] + (WP[1])*vrrs[2,1,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[6,1,1] = (PA[2])*vrrs[2,1,1] + (WP[2])*vrrs[2,1,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[7,1,1] = (PA[3])*vrrs[2,1,1] + (WP[3])*vrrs[2,1,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[8,1,1] = (PA[2])*vrrs[3,1,1] + (WP[2])*vrrs[3,1,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[9,1,1] = (PA[3])*vrrs[3,1,1] + (WP[3])*vrrs[3,1,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])
    vrrs[10,1,1] =(PA[3])*vrrs[4,1,1] + (WP[3])*vrrs[4,1,2] +
              1/(2*zeta)*(vrrs[1,1,1]-eta/ze*vrrs[1,1,2])

    vrrs[5,2,1] = (PA[1])*vrrs[2,2,1] + (WP[1])*vrrs[2,2,2] +
              1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[6,2,1] = (PA[2])*vrrs[2,2,1] + (WP[2])*vrrs[2,2,2] +
              1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[7,2,1] = (PA[3])*vrrs[2,2,1] + (WP[3])*vrrs[2,2,2] +
              1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[8,2,1] = (PA[2])*vrrs[3,2,1] + (WP[2])*vrrs[3,2,2] +
              1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[9,2,1] = (PA[3])*vrrs[3,2,1] + (WP[3])*vrrs[3,2,2] +
              1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[10,2,1] =(PA[3])*vrrs[4,2,1] + (WP[3])*vrrs[4,2,2] +
              1/(2*zeta)*(vrrs[1,2,1]-eta/ze*vrrs[1,2,2]) + 0.5/ze*vrrs[1,1,2]

    vrrs[5,3,1] = (PA[1])*vrrs[2,3,1] + (WP[1])*vrrs[2,3,2] +
              1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[6,3,1] = (PA[2])*vrrs[2,3,1] + (WP[2])*vrrs[2,3,2] +
              1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[7,3,1] = (PA[3])*vrrs[2,3,1] + (WP[3])*vrrs[2,3,2] +
              1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[8,3,1] = (PA[2])*vrrs[3,3,1] + (WP[2])*vrrs[3,3,2] +
              1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[9,3,1] = (PA[3])*vrrs[3,3,1] + (WP[3])*vrrs[3,3,2] +
              1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[10,3,1] =(PA[3])*vrrs[4,3,1] + (WP[3])*vrrs[4,3,2] +
              1/(2*zeta)*(vrrs[1,3,1]-eta/ze*vrrs[1,3,2]) + 0.5/ze*vrrs[1,1,2]

    vrrs[5,4,1] = (PA[1])*vrrs[2,4,1] + (WP[1])*vrrs[2,4,2] +
              1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[6,4,1] = (PA[2])*vrrs[2,4,1] + (WP[2])*vrrs[2,4,2] +
              1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[7,4,1] = (PA[3])*vrrs[2,4,1] + (WP[3])*vrrs[2,4,2] +
              1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[8,4,1] = (PA[2])*vrrs[3,4,1] + (WP[2])*vrrs[3,4,2] +
              1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[9,4,1] = (PA[3])*vrrs[3,4,1] + (WP[3])*vrrs[3,4,2] +
              1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]
    vrrs[10,4,1] =(PA[3])*vrrs[4,4,1] + (WP[3])*vrrs[4,4,2] +
              1/(2*zeta)*(vrrs[1,4,1]-eta/ze*vrrs[1,4,2]) + 0.5/ze*vrrs[1,1,2]

    vrrs[4,5,1] = (QC[1])*vrrs[4,2,1] + (WQ[1])*vrrs[4,2,2] +
              1/(2*zeta)*(vrrs[4,1,1]-eta/ze*vrrs[4,1,2]) + 0.5/ze*vrrs[1,1,2]

    return vrrs[:,:,1]
end
