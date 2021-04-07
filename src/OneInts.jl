using SpecialFunctions

export overlap, kinetic, nuclear_attraction

function overlap(a::PGBF,b::PGBF)
    return a.norm*b.norm*overlap(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

overlap(a::CGBF,b::CGBF) = contract(overlap,a,b)

function overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    gamma = aexpn+bexpn
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    rab2 = dist2(ax-bx,ay-by,az-bz) 
    pre = (pi/gamma)^1.5*exp(-aexpn*bexpn*rab2/gamma)
    wx = overlap1d(aI,bI,px-ax,px-bx,gamma)
    wy = overlap1d(aJ,bJ,py-ay,py-by,gamma)
    wz = overlap1d(aK,bK,pz-az,pz-bz,gamma)
    return pre*wx*wy*wz
end

function gaussian_product_center(a::PGBF,b::PGBF)
    return (a.expn*[a.x,a.y,a.z]+b.expn*[b.x,b.y,b.z])/(a.expn+b.expn)
end

function gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
         return (aexpn*[ax,ay,az]+bexpn*[bx,by,bz])/(aexpn+bexpn)    
end

function overlap1d(la,lb,ax,bx,gamma)
    total = 0
    for i in 0:div(la+lb,2)
        total += binomial_prefactor(2i,la,lb,ax,bx)*factorial2(2i-1)/(2gamma)^i
    end
    return total
end

#=function binomial_prefactor(s,ia,ib,xpa,xpb)
    total = 0
    for t in 0:s 
        if (s-ia) <= t <= ib
            total += binomial_kernel(ia,s-t,xpa)*binomial_kernel(ib,t,xpb) 
        end
    end
    return total
end=#
binomial_prefactor(s,ia,ib,xpa,xpb) = sum(binomial_kernel(ia,s-t,xpa)*binomial_kernel(ib,t,xpb) for t in 0:s if (s-ia) <= t <= ib)
binomial_kernel(i,t,x) = binomial(i,t)x^(i-t)

function kinetic(a::PGBF,b::PGBF)
    return a.norm*b.norm*kinetic(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
                                b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

function kinetic(aexpn,ax,ay,az,
                 aI,aJ,aK,bexpn,bx,
                 by,bz,bI,bJ,bK)
    overlap0 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    overlapx1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI+2,bJ,bK)
    overlapy1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ+2,bK)
    overlapz1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK+2)
    overlapx2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI-2,bJ,bK)
    overlapy2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ-2,bK)
    overlapz2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK-2)
    term0 = bexpn*(2*(bI+bJ+bK)+3)*overlap0
    term1 = -2*(bexpn^2)*(overlapx1+overlapy1+overlapz1)
    term2 = -0.5*(bI*(bI-1)*overlapx2+bJ*(bJ-1)*overlapy2+bK*(bK-1)*overlapz2)
    return term0+term1+term2
end

kinetic(a::CGBF,b::CGBF) = contract(kinetic,a,b)

function Aterm(i,r,u,l1,l2,ax,bx,cx,gamma)
    term1 = (-1)^i*binomial_prefactor(i,l1,l2,ax,bx)
    term2 = (-1)^u*factorial(i)*cx^(i-2r-2u)
    term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
    return term1*term2*term3
end

function Aarray(l1,l2,a,b,c,g)
    Imax = l1+l2+1
    A = zeros(Float64,Imax)
    for i in 0:(Imax-1)
        for r in 0:div(i,2)
            for u in 0:div(i-2r,2)
                I = i-2r-u+1
                A[I] += Aterm(i,r,u,l1,l2,a,b,c,g)
            end
        end
    end
    return A
end

function nuclear_attraction(aexpn,ax,ay,az,
                            aI,aJ,aK,
                            bexpn,bx,by,bz,
                            bI,bJ,bK,
                            cx,cy,cz)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    gamma = aexpn+bexpn
    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcp2 = dist2(cx-px,cy-py,cz-pz)
    Ax = Aarray(aI,bI,px-ax,px-bx,px-cx,gamma)
    Ay = Aarray(aJ,bJ,py-ay,py-by,py-cy,gamma)
    Az = Aarray(aK,bK,pz-az,pz-bz,pz-cz,gamma)
    total = 0
    for I in 0:(aI+bI)
        for J in 0:(aJ+bJ)
            for K in 0:(aK+bK)
                total += Ax[I+1]*Ay[J+1]*Az[K+1]*Fgamma(I+J+K,rcp2*gamma)
            end
        end
    end
    val=-2pi*exp(-aexpn*bexpn*rab2/gamma)*total/gamma
    return val
end

function nuclear_attraction(a::PGBF,b::PGBF,cx,cy,cz)
    return a.norm*b.norm*nuclear_attraction(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
                                            b.expn,b.x,b.y,b.z,b.I,b.J,b.K,cx,cy,cz)
end
nuclear_attraction(a::PGBF,b::PGBF,c::Atom) = c.atno*nuclear_attraction(a,b,c.x,c.y,c.z)
nuclear_attraction(a::PGBF,b::PGBF,m::Molecule) = sum([nuclear_attraction(a,b,c) for c in m.atomlist])

"Boys Fgamma function, using the lower incomplete gamma function."
function Fgamma(m,x,SMALL=1e-18)
    x = max(x,SMALL) # Evidently needs underflow protection
    return 0.5*x^(-m-0.5)*gammainc(m+0.5,x)
end

"gammainc returns the lower incomplete gamma function"
gammainc(a,x) = gamma(a)*gamma_inc(a,x)[1]

#= Commenting out old gammainc code:

function gammainc(a,x)
    # This is the series version of gamma from pyquante. For reasons I
    # don't get, it doesn't work around a=1. This works alright, but
    # is only a stopgap solution until Julia gets an incomplete gamma
    # function programmed
    if abs(a-1) < 1e-3
        println("Warning: gammainc_series is known to have problems for a ~ 1")
    end
    if x < (a+1.0)
        #Use the series representation
        gam,gln = gser(a,x)
    else 
        #Use continued fractions
        gamc,gln = gcf(a,x)
        gam = 1-gamc
    end
    return exp(gln)*gam
end

function gser(a,x,ITMAX=100,EPS=3e-9)
    # Series representation of Gamma. NumberRec sect 6.1.
    gln=loggamma(a)
    if x == 0
        return 0,gln
    end
    ap = a
    delt = s = 1/a
    for i in 1:ITMAX
        ap += 1
        delt *= (x/ap)
        s += delt
        if abs(delt) < abs(s)*EPS
            break
        end
    end
    return s*exp(-x+a*log(x)-gln),gln
end

function gcf(a,x,ITMAX=200,EPS=3e-9,FPMIN=1e-30)
    #Continued fraction representation of Gamma. NumRec sect 6.1"
    gln=loggamma(a)
    b=x+1.0-a
    c=1.0/FPMIN
    d=1.0/b
    h=d
    for i in 1:ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if abs(d) < FPMIN
            d=FPMIN
        end
        c=b+an/c
        if abs(c) < FPMIN
            c=FPMIN
        end
        d=1.0/d
        delt=d*c
        h=h*delt
        if abs(delt-1.) < EPS
            break
        end
    end
    gammcf = exp(-x+a*log(x)-gln)*h
    return gammcf,gln
end
=#

# Need a nested scope to squeeze this into the contract function
function nuclear_attraction(a::CGBF,b::CGBF,cx,cy,cz)
    na(a,b) = nuclear_attraction(a,b,cx,cy,cz)
    contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,c::Atom)
    na(a,b) = nuclear_attraction(a,b,c)
    contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,m::Molecule)
    na(a,b) = nuclear_attraction(a,b,m)
    contract(na,a,b)
end

function all_1e_ints(bfs::BasisSet,mol::Molecule)
    n = length(bfs.bfs)
    S = Array{Float64}(n,n)
    T = Array{Float64}(n,n)
    V = Array{Float64}(n,n)
    for (i,j) in pairs(n)
        a,b = bfs.bfs[i],bfs.bfs[j]
        S[i,j] = S[j,i] = overlap(a,b)
        T[i,j] = T[j,i] = kinetic(a,b)
        V[i,j] = V[j,i] = nuclear_attraction(a,b,mol)
    end
    return S,T,V
end
