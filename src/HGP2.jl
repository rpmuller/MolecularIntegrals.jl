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
function ssss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn,dxyz,mmax=0)
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
    return Kab*Kcd/sqrt(zeta+eta)*[Fgamma(m,T) for m in 0:mmax]   # HGP eq 12
end

function psss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn,dxyz,mmax=0)
    sarray = ssss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn, dxyz,mmax+1)
    # Recalculate a number of terms from ssss. When the code works, be more
    # judicious in what we recalculate.
    values = zeros(Float64,(3,mmax+1))
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    qxyz = gaussian_product_center(cexpn,cxyz,dexpn,dxyz)
    wxyz = gaussian_product_center(zeta,pxyz,eta,qxyz)
    for i in 1:3
        for m in 0:mmax
            values[i,m+1] = (pxyz[i]-axyz[i])*sarray[m+1] + (wxyz[i]-pxyz[i])*sarray[m+2]
        end
    end
    return values
end

function psps(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn,dxyz,mmax=0)
    sarray = ssss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn, dxyz,mmax+1)
    parray = psss(aexpn,axyz, bexpn,bxyz, cexpn,cxyz, dexpn, dxyz,mmax+1)
    # Recalculate a number of terms from ssss. When the code works, be more
    # judicious in what we recalculate.
    values = zeros(Float64,(3,3,mmax+1))
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    qxyz = gaussian_product_center(cexpn,cxyz,dexpn,dxyz)
    wxyz = gaussian_product_center(zeta,pxyz,eta,qxyz)
    for i in 1:3
        for j in 1:3
            for m in 0:mmax
                values[i,j,m+1] = (qxyz[j]-bxyz[j])*parray[i,m+1] + (wxyz[j]-qxyz[j])*parray[i,m+2] + 1/2/(zeta+eta)*sarray[m+2] 
            end
        end
    end
    return values
end

shell_indices = Dict(
    0 => [[0,0,0]], # 1
    1 => [[1,0,0],[0,1,0],[0,0,1]], # 3
    2 => [[2,0,0],[1,1,0],[1,0,1],[0,2,0],[0,1,1],[0,0,2]],
    3 => [[3,0,0],[2,1,0],[2,0,1],
            [1,2,0],[1,0,2],[1,1,1],
            [0,3,0],[0,2,1],[0,1,2],[0,0,3]], # 10
    4 => [[4,0,0],[3,1,0],[3,0,1],[2,2,0],[2,1,1],[2,0,2],
            [1,3,0],[1,2,1],[1,1,2],[1,0,3],
            [0,4,0],[0,3,1],[0,2,2],[0,1,3],[0,0,4]] # 15
)

"vrrindices - Generate indices for vrr in three steps:
- (0,0,0,0,0,0,m),  
- (ax,ay,az,0,0,0,m),  
- (ax,ay,az,cx,cy,cz,m) 
This version also adjusts the m terms based on a,c. This is the most efficient way to 
go through the recurrance relationships in practics."
function vrrindices(amax,cmax)
    indices = []
    mmax=amax+cmax
    # First generate (0,0,0, 0,0,0, m) 
    c=a=0
    for m in 0:(mmax-a-c) 
        push!(indices,(0,0,0,0,0,0,m))
    end

    # Now generate (ax,ay,az,0,0,0,m) 
    for a in 1:amax
        for (ax,ay,az) in shell_indices[a]
            for m in 0:(mmax-a-c)
                push!(indices,(ax,ay,az,0,0,0,m))
            end
        end
    end

    # Now build (ax,ay,az,cx,cy,cz,m)
    for a in 0:amax
        for (ax,ay,az) in shell_indices[a]
            for c in 1:cmax
                for (cx,cy,cz) in shell_indices[c]
                    for m in 0:(mmax-a-c)
                        push!(indices,(ax,ay,az,cx,cy,cz,m))
                    end
                end
            end
        end 
    end
    return indices
end

"prunem - Keep only the dictionary keys with m (last index) = 0"
prunem(d::Dict) = Dict(k[1:end-1] => v for (k,v) in d if k[end] == 0)

"vrr2 - iterative version of HGP vertical recurrance relations"
function vrr2(amax,cmax, aexpn,bexpn,cexpn,dexpn, axyz,bxyz,cxyz,dxyz)
    values = Dict()
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    qxyz = gaussian_product_center(cexpn,cxyz,dexpn,dxyz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    wxyz = gaussian_product_center(zeta,pxyz,eta,qxyz)
    rab2 = dist2(axyz-bxyz)
    rcd2 = dist2(cxyz-dxyz)
    rpq2 = dist2(pxyz-qxyz)
    T = zeta*eta/(zeta+eta)*rpq2
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    mmax=amax+cmax

    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (0,0,0, 0,0,0, m) using eq 12
    c=a=0
    for m in 0:(mmax-a-c) 
        values[(0,0,0, 0,0,0, m)] = Kab*Kcd/sqrt(zeta+eta)*Fgamma(m,T)
    end

    # Now generate (ax,ay,az,0,0,0,m) 
    # Eq 6a, with c=0 also:
    #   [a+1,0]m = (Pi-Ai)[a,0]m + (Wi-Pi)[a,0]m+1 
    #        + a_i/2zeta ([a-1,0]m - eta/zeta+eta[a-1,0]m+1)        # eq 6c
    for a in 1:amax
        for av in shell_indices[a]
            # TODO: find a nicer way to compute the index differences. Currently I compute the argmax of the (ax,ay,az) values,
            #   and use that to compute m = (ax,ay-1,az) or in whatever direction. But this is ugly code. I wrote the `unit` code
            #   to try to make this more idiomatic julia, but I'm not using it now.
            i = argmax(av)
            am = copy(av)
            am[i] -= 1
            am2 = copy(am)
            am2[i] -= 1
            for m in 0:(mmax-a-c)
                values[(av[1],av[2],av[3],0,0,0,m)] = (pxyz[i]-axyz[i])*values[(am[1],am[2],am[3],0,0,0,m)]+(wxyz[i]-pxyz[i])*values[(am[1],am[2],am[3],0,0,0,m+1)]
                if am2[i] >= 0
                    values[(av[1],av[2],av[3],0,0,0,m)] += am[i]/(2*zeta)*(values[(am2[1],am2[2],am2[3],0,0,0,m)]-eta/ze*values[(am2[1],am2[2],am2[3],0,0,0,m+1)])
                end
            end
        end
    end

    # Now build (ax,ay,az,cx,cy,cz,m)
    # The c-based version of 6a is:
    #   [a,c+1]m = (Qi-Bi)[a,c]m + (Wi-Qi)[a,c]m+1
    #       + c_i/2zeta ([a,c-1]m - eta/zeta+eta[a,c-1]m+1)         # eq 6b
    #       + a_i/2(zeta+eta)[a-1,c]m+1
    for a in 0:amax
        for av in shell_indices[a]
            i = argmax(av)
            am = copy(av)
            am[i] -= 1
            am2 = copy(am)
            am2[i] -= 1
            for c in 1:cmax
                for cv in shell_indices[c]
                    j = argmax(cv)
                    cm = copy(cv)
                    cm[i] -= 1
                    cm2 = copy(cm)
                    cm2[i] -= 1
                    for m in 0:(mmax-a-c)
                        values[(av[1],av[2],av[3],cv[1]cv[2],cv[3],m)] = (qxyz[j]-bxyz[j])*values[(av[1],av[2],av[3],cm[1],cm[2],cm[3],m)]
                            +(wxyz[j]-qxyz[j])*values[(av[1],av[2],av[3],cm[1],cm[2],cm[3],m+1)]
                            if cm2[j] >= 0
                                values[(av[1],av[2],av[3],0,0,0,m)] += am[i]/(2*zeta)*(values[(am2[1],am2[2],am2[3],0,0,0,m)]-eta/ze*values[(am2[1],am2[2],am2[3],0,0,0,m+1)])
                            end                                    
                    end
                end
            end
        end 
    end
    return prunem(values)
end

"unit(n,d) - create a n-dim unit vector in direction d"
function unit(n,d) 
    v = zeros(Int,n)
    v[d] = 1
    return v
end


