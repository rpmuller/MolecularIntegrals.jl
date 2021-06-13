export coulomb

using OffsetArrays

#using StaticArrays # Surprisingly, didn't speed up at all
"""
    coulomb(aexpn,ax,ay,az,aI,aJ,aK,
        bexpn,bx,by,bz,bI,bJ,bK,
        cexpn,cx,cy,cz,cI,cJ,cK,
        dexpn,dx,dy,dz,dI,dJ,dK)

Compute the coulomb repulsion between four primitive Gaussian basis 
functions defined by their exponents, xyz coordinates, and IJK powers.        
"""
function coulomb(aexpn,ax,ay,az,aI,aJ,aK,
    bexpn,bx,by,bz,bI,bJ,bK,
    cexpn,cx,cy,cz,cI,cJ,cK,
    dexpn,dx,dy,dz,dI,dJ,dK)

    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)

    g1 = aexpn+bexpn
    g2 = cexpn+dexpn
    delta = 0.25*(1/g1+1/g2)
    
    Bx = Barray(aI,bI,cI,dI,px,ax,bx,qx,cx,dx,g1,g2,delta)
    By = Barray(aJ,bJ,cJ,dJ,py,ay,by,qy,cy,dy,g1,g2,delta)
    Bz = Barray(aK,bK,cK,dK,pz,az,bz,qz,cz,dz,g1,g2,delta)

    x = 0.25*rpq2/delta

    s = 0
    for I in 0:(aI+bI+cI+dI)
        for J in 0:(aJ+bJ+cJ+dJ)
            for K in 0:(aK+bK+cK+dK)
                s += Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,x)
            end
        end
    end
    return 2*pi^(2.5)/(g1*g2*sqrt(g1+g2))*exp(-aexpn*bexpn*rab2/g1)*exp(-cexpn*dexpn*rcd2/g2)*s
end


"""
    coulomb(a::PGBF,b::PGBF,c::PGBF,d::PGBF)

Compute the coulomb repulsion between four primitive Gaussian basis 
functions a,b,c,d.        
"""
function coulomb(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
    return a.norm*b.norm*c.norm*d.norm*coulomb(a.expn,a.xyz...,a.I,a.J,a.K,
        b.expn,b.xyz...,b.I,b.J,b.K,
        c.expn,c.xyz...,c.I,c.J,c.K,
        d.expn,d.xyz...,d.I,d.J,d.K)
end

fB(i,l1,l2,p,a,b,r,g) = binomial_prefactor(i,l1,l2,p-a,p-b)*B0(i,r,g)
B0(i,r,g) = fact_ratio2(i,r)*(4g)^(r-i)
#fact_ratio2(a,b) = factorial(a)/factorial(b)/factorial(a-2b)
fact_ratio2(a,b) = prod((b+1):a)/factorial(a-2b) # faster

function Bterm(i1,i2,r1,r2,u,l1,l2,l3,l4,
        Px,Ax,Bx,Qx,Cx,Dx,
        gamma1,gamma2,delta)
    # THO eq. 2.22
    return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)*
        (-1)^i2*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)*
        (-1)^u*fact_ratio2(i1+i2-2*(r1+r2),u)*
        (Qx-Px)^(i1+i2-2*(r1+r2)-2*u)/delta^(i1+i2-2*(r1+r2)-u)
end

function Barray(l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta)
    Imax = l1+l2+l3+l4+1
    B = OffsetArray(zeros(Float64,Imax),0:(Imax-1))
    for i1 in 0:(l1+l2)
        for i2 in 0:(l3+l4)
            for r1 in 0:div(i1,2)
                for r2 in 0:div(i2,2)
                    for u in 0:(div(i1+i2,2)-r1-r2)
                        I = i1+i2-2*(r1+r2)-u
                        B[I] += Bterm(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                      p,a,b,q,c,d,g1,g2,delta)
                    end
                end
            end
        end
    end
    return B
end

"""
    coulomb(a::CGBF,b::CGBF,c::CGBF,d::CGBF)

Compute the coulomb repulsion between four contracted Gaussian basis 
functions a,b,c,d.        
"""
coulomb(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb,a,b,c,d)

function all_twoe_ints(bfs,ERI=coulomb)
    n = length(bfs)
    totlen = div(n*(n+1)*(n*n+n+2),8)
    ints2e = zeros(Float64,totlen)
    Threads.@threads for (i,j,k,l) in collect(iiterator(n))
        ints2e[iindex(i,j,k,l)] = ERI(bfs[i],bfs[j],bfs[k],bfs[l])
    end
    return ints2e
end

function make2JmK(D::Array{Float64,2},Ints::Array{Float64,1})
    n = size(D,1)
    G = Array{Float64}(n,n)
    D1 = reshape(D,n*n)
    temp = Array{Float64}(n*n)
    for (i,j) in pairs(n)
        kl = 1
        for (k,l) in rpairs(n)
            temp[kl] = 2*Ints[iindex(i,j,k,l)]-Ints[iindex(i,k,j,l)]
            kl += 1
        end
        G[i,j] = G[j,i] = dot(D1,temp)
    end
    return G
end

dmat(U::Array{Float64,2},nocc::Int64) = U[:,1:nocc]*U[:,1:nocc]'
