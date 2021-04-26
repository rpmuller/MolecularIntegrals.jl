using SpecialFunctions
using OffsetArrays

export overlap, kinetic, nuclear_attraction

"""
overlap(a::PGBF,b::PGBF)

Compute the overlap between primitive Gaussian basis functions `a` and `b`.
"""
function overlap(a::PGBF,b::PGBF)
    return a.norm*b.norm*overlap(a.expn,a.xyz,a.I,a.J,a.K,
    b.expn,b.xyz,b.I,b.J,b.K)
end

"""
overlap(a::CGBF,b::CGBF)

Compute the overlap between contracted Gaussian basis functions `a` and `b`.
"""
overlap(a::CGBF,b::CGBF) = contract(overlap,a,b)

"""
overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)

Compute the overlap between primitive Gaussian basis functions 
defined by centers `ax,ay,az`, `bx,by,bz`, powers `aI,aJ,aK`
`bI,bJ,bK`, and exponents `aexpn,bexpn`.
"""
function overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    return overlap(aexpn,[ax,ay,az],aI,aJ,aK,bexpn,[bx,by,bz],bI,bJ,bK)
end

"""
overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK)

Compute the overlap between primitive Gaussian basis functions 
defined by centers `axyz`, `bxyz`, powers `aI,aJ,aK`
`bI,bJ,bK`, and exponents `aexpn,bexpn`.
"""
function overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK)
    gamma = aexpn+bexpn
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    pa = pxyz-axyz
    pb = pxyz-bxyz
    rab2 = dist2(axyz-bxyz) 
    pre = (pi/gamma)^1.5*exp(-aexpn*bexpn*rab2/gamma)
    wx = overlap1d(aI,bI,pa[1],pb[1],gamma)
    wy = overlap1d(aJ,bJ,pa[2],pb[2],gamma)
    wz = overlap1d(aK,bK,pa[3],pb[3],gamma)
    return pre*wx*wy*wz
end

"""
gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)

Compute the Gaussian function center defined by the 
product of Gaussian functions at  `ax,ay,az`, `bx,by,bz`,
and exponents `aexpn,bexpn`.
"""
gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz) = gaussian_product_center(aexpn,[ax,ay,az],bexpn,[bx,by,bz])

"""
gaussian_product_center(a::PGBF,b::PGBF)

Compute the Gaussian function center defined by the 
product of Gaussian functions `a` and `b`.
"""
gaussian_product_center(a::PGBF,b::PGBF) = gaussian_product_center(a.expn,a.xyz,b.expn,b.xyz)

"""
gaussian_product_center(aexpn,axyz,bexpn,bxyz)

Compute the Gaussian function center defined by the 
product of Gaussian functions at  `axyz`, `bxyz`,
and exponents `aexpn,bexpn`.
"""
function gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    ab = aexpn+bexpn
    a = aexpn/ab
    b = bexpn/ab
    return a*axyz+b*bxyz
end

function overlap1d(la,lb,ax,bx,gamma)
    total = 0
    for i in 0:div(la+lb,2)
        total += binomial_prefactor(2i,la,lb,ax,bx)*factorial2(2i-1)/(2gamma)^i
    end
    return total
end

binomial_prefactor(s,ia,ib,xpa,xpb) = sum(binomial_kernel(ia,s-t,xpa)*binomial_kernel(ib,t,xpb) for t in 0:s if (s-ia) <= t <= ib)
binomial_kernel(i,t,x) = binomial(i,t)x^(i-t)

"""
kinetic(a::PGBF,b::PGBF)

Compute the kinetic energy between primitive Gaussian basis functions `a` and `b`.
"""
function kinetic(a::PGBF,b::PGBF)
    return a.norm*b.norm*kinetic(a.expn,a.xyz,a.I,a.J,a.K,
                                b.expn,b.xyz,b.I,b.J,b.K)
end

"""
kinetic(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)

Compute the kinetic energy between primitive Gaussian basis functions 
defined by centers `ax,ay,az`, `bx,by,bz`, powers `aI,aJ,aK`
`bI,bJ,bK`, and exponents `aexpn,bexpn`.
"""
function kinetic(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    return kinetic(aexpn,[ax,ay,az],aI,aJ,aK,bexpn,[bx,by,bz],bI,bJ,bK)
end

"""
kinetic(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK)

Compute the kinetic energy between primitive Gaussian basis functions 
defined by centers `ax,yz`, `bxyz`, powers `aI,aJ,aK`
`bI,bJ,bK`, and exponents `aexpn,bexpn`.
"""
function kinetic(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK)
    overlap0 = overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK)
    overlapx1 = overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI+2,bJ,bK)
    overlapy1 = overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ+2,bK)
    overlapz1 = overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK+2)
    overlapx2 = overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI-2,bJ,bK)
    overlapy2 = overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ-2,bK)
    overlapz2 = overlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK-2)
    term0 = bexpn*(2*(bI+bJ+bK)+3)*overlap0
    term1 = -2*(bexpn^2)*(overlapx1+overlapy1+overlapz1)
    term2 = -0.5*(bI*(bI-1)*overlapx2+bJ*(bJ-1)*overlapy2+bK*(bK-1)*overlapz2)
    return term0+term1+term2
end

"""
kinetic(a::CGBF,b::CGBF)

Compute the kinetic energy between contracted Gaussian basis functions `a` and `b`.
"""
kinetic(a::CGBF,b::CGBF) = contract(kinetic,a,b)

function Aterm(i,r,u,l1,l2,ax,bx,cx,gamma)
    term1 = (-1)^i*binomial_prefactor(i,l1,l2,ax,bx)
    term2 = (-1)^u*factorial(i)*cx^(i-2r-2u)
    term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
    return term1*term2*term3
end

function Aarray(l1,l2,a,b,c,g)
    Imax = l1+l2+1
    A = OffsetArray(zeros(Float64,Imax),0:(Imax-1))
    for i in 0:(Imax-1)
        for r in 0:div(i,2)
            for u in 0:div(i-2r,2)
                I = i-2r-u
                A[I] += Aterm(i,r,u,l1,l2,a,b,c,g)
            end
        end
    end
    return A
end

"""
nuclear_attraction(a::PGBF,b::PGBF,cxyz)

Compute the nuclear attraction energy between primitive Gaussian basis functions `a` and `b`
and center `cxyz`.
"""
function nuclear_attraction(a::PGBF,b::PGBF,cxyz)
    return a.norm*b.norm*nuclear_attraction(a.expn,a.xyz,a.I,a.J,a.K,
                                            b.expn,b.xyz,b.I,b.J,b.K,cxyz)
end

"""
nuclear_attraction(a::PGBF,b::PGBF,c::Atom)

Compute the nuclear attraction energy between primitive Gaussian basis functions `a` and `b`
and atom `c`.
"""
nuclear_attraction(a::PGBF,b::PGBF,c::Atom) = c.atno*nuclear_attraction(a,b,c.xyz)

"""
nuclear_attraction(a::PGBF,b::PGBF,m::Vector{Atom})

Compute the sum of nuclear attraction energy between primitive Gaussian basis 
functions `a` and `b` and the vector of atoms in `m`.
"""
nuclear_attraction(a::PGBF,b::PGBF,m::Vector{Atom}) = sum([nuclear_attraction(a,b,c) for c in m])

"""
nuclear_attraction(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,cx,cy,cz)

Compute the nuclear attraction energy between primitive Gaussian basis functions 
defined by `ax,ay,az`, `bx,by,bz`, powers `aI,aJ,aK`, `bI,bJ,bK`, and exponents `aexpn,bexpn`
and center `cxyz`.
"""
function nuclear_attraction(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,cx,cy,cz)
    return nuclear_attraction(aexpn,[ax,ay,az],aI,aJ,aK,bexpn,[bx,by,bz],bI,bJ,bK,[cx,cy,cz])
end

"""
nuclear_attraction(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,cx,cy,cz)

Compute the nuclear attraction energy between primitive Gaussian basis functions 
defined by `ax,ay,az`, `bx,by,bz`, powers `aI,aJ,aK`, `bI,bJ,bK`, and exponents `aexpn,bexpn`
and center `cxyz`.
"""
function nuclear_attraction(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK,cxyz)
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    gamma = aexpn+bexpn
    pa = pxyz-axyz
    pb = pxyz-bxyz
    pc = pxyz-cxyz
    rab2 = dist2(axyz-bxyz)
    rcp2 = dist2(pc)
    Ax = Aarray(aI,bI,pa[1],pb[1],pc[1],gamma)
    Ay = Aarray(aJ,bJ,pa[2],pb[2],pc[2],gamma)
    Az = Aarray(aK,bK,pa[3],pb[3],pc[3],gamma)
    total = 0
    for I in 0:(aI+bI)
        for J in 0:(aJ+bJ)
            for K in 0:(aK+bK)
                total += Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma)
            end
        end
    end
    return -2pi*exp(-aexpn*bexpn*rab2/gamma)*total/gamma
end

"Boys Fgamma function, using the lower incomplete gamma function."
function Fgamma(m,x,SMALL=1e-18)
    mhalf = m+0.5
    x = max(x,SMALL) # Evidently needs underflow protection
    return 0.5*x^-mhalf*gammainc(mhalf,x)
end

"gammainc returns the lower incomplete gamma function"
gammainc(a,x) = gamma(a)*gamma_inc(a,x)[1]

"""
nuclear_attraction(a::CGBF,b::CGBF,cxyz)

Compute the nuclear attraction energy between contracted Gaussian 
basis functions `a` and `b` and center `cxyz`.
"""
function nuclear_attraction(a::CGBF,b::CGBF,cxyz)
    na(a,b) = nuclear_attraction(a,b,cxyz)
    contract(na,a,b)
end

"""
nuclear_attraction(a::CGBF,b::CGBF,c::Atom)

Compute the nuclear attraction energy between contracted Gaussian 
basis functions `a` and `b` and atom `c`.
"""
function nuclear_attraction(a::CGBF,b::CGBF,c::Atom)
    na(a,b) = nuclear_attraction(a,b,c)
    contract(na,a,b)
end

"""
nuclear_attraction(a::CGBF,b::CGBF,m::Vector{Atom})

Compute the sum of nuclear attraction energy between contracted Gaussian 
basis functions `a` and `b` and vector of atoms in `m`.
"""
function nuclear_attraction(a::CGBF,b::CGBF,m::Vector{Atom})
    na(a,b) = nuclear_attraction(a,b,m)
    contract(na,a,b)
end

"""
all_1e_ints(bfs::Vector{CGBF},mol::Vector{Atom})

Return the overlap, the kinetic energy, and the nuclear attraction
integrals between vector of contracted functions `bfs` and
atoms `mol`.
"""
function all_1e_ints(bfs::Vector{CGBF},mol::Vector{Atom})
    n = length(bfs)
    S = Array{Float64}(n,n)
    T = Array{Float64}(n,n)
    V = Array{Float64}(n,n)
    for (i,j) in pairs(n)
        a,b = bfs[i],bfs[j]
        S[i,j] = S[j,i] = overlap(a,b)
        T[i,j] = T[j,i] = kinetic(a,b)
        V[i,j] = V[j,i] = nuclear_attraction(a,b,mol)
    end
    return S,T,V
end
