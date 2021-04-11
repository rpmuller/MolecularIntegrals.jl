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
    ze = zeta+eta
    pxyz = gaussian_product_center(aexpn,axyz,bexpn,bxyz)
    qxyz = gaussian_product_center(cexpn,cxyz,dexpn,dxyz)
    wxyz = gaussian_product_center(zeta,pxyz,eta,qxyz)
    for i in 1:3
        for j in 1:3
            for m in 0:mmax
                values[i,j,m+1] = (qxyz[j]-bxyz[j])*parray[i,m+1] + (wxyz[j]-qxyz[j])*parray[i,m+2] + 1/(2ze)*sarray[m+2] 
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
    for m in 0:mmax
        values[(0,0,0, 0,0,0, m)] = Kab*Kcd*Fgamma(m,T)/sqrt(ze)
    end

    # Now generate (ax,ay,az,0,0,0,m) 
    # Eq 6a, with c=0 also:
    #   [a+1,0]m = (Pi-Ai)[a,0]m + (Wi-Pi)[a,0]m+1 
    #        + a_i/2zeta ([a-1,0]m - eta/zeta+eta[a-1,0]m+1)        # eq 6b
    for a in 1:amax
        for ap in shell_indices[a]
            apx,apy,apz = ap
            i = argmax(ap) # Choose argmax(ap) as the direction to use for building new terms
            av,am = vdiffs(ap,i) #am = ap-1i, am2 = ap-2i in eq 6-7; i direction of change
            avx,avy,avz = av
            amx,amy,amz = am
            for m in 0:(mmax-a)
                values[(apx,apy,apz,0,0,0,m)] = (pxyz[i]-axyz[i])*values[(avx,avy,avz,0,0,0,m)]+(wxyz[i]-pxyz[i])*values[(avx,avy,avz,0,0,0,m+1)]
                if am[i] >= 0
                    values[(apx,apy,apz,0,0,0,m)] += av[i]/(2*zeta)*(values[(amx,amy,amz,0,0,0,m)]-eta/ze*values[(amx,amy,amz,0,0,0,m+1)])
                end
            end
        end
    end

    # Next build (0,0,0,cx,cy,cz,m)
    # The c-based version of 6a is:
    #   [0,c+1]m = (Qj-Bi)[0,c]m + (Wj-Qj)[0,c]m+1
    #       + c_j/2eta ([0,c-1]m - zeta/zeta+eta[0,c-1]m+1)         # eq 6c
    # 
    for c in 1:cmax
        for cp in shell_indices[c]
            cpx,cpy,cpz = cp
            j = argmax(cp)  # Choose argmax(cp) as the direction to use for building new terms
            cv,cm = vdiffs(cp,j)
            cvx,cvy,cvz = cv
            cmx,cmy,cmz = cm
            for m in 0:(mmax-c)
                values[(0,0,0,cpx,cpy,cpz,m)] = (qxyz[j]-bxyz[j])*values[(0,0,0,cvx,cvy,cvz,m)]+(wxyz[j]-qxyz[j])*values[(0,0,0,cvx,cvy,cvz,m+1)]
                if cm[j] >= 0
                    values[(0,0,0,cpx,cpy,cpz,m)] += cv[j]/(2*eta)*(values[(0,0,0,cmx,cmy,cmz,m)]-zeta/ze*values[(0,0,0,cmx,cmy,cmz,m+1)])
                end
            end
        end
    end

    # Now build (ax,ay,az,cx,cy,cz,m)
    # The c-based version of 6a is:
    #   [a,c+1]m = (Qj-Bi)[a,c]m + (Wj-Qj)[a,c]m+1
    #       + c_j/2eta ([a,c-1]m - zeta/zeta+eta[a,c-1]m+1)         # eq 6d
    #       + a_j/2(zeta+eta)[a-1,c]m+1
    for a in 1:amax
        for av in shell_indices[a]
            avx,avy,avz = av
            for c in 1:cmax
                for cp in shell_indices[c]
                    cpx,cpy,cpz = cp
                    j = argmax(cp)  # Choose argmax(cp) as the direction to use for building new terms
                    cv,cm = vdiffs(cp,j)
                    cvx,cvy,cvz = cv
                    cmx,cmy,cmz = cm

                    am,am2 = vdiffs(av,j) #am = av-1j, am2 = av-2j in eq 6-7
                    amx,amy,amz = am
                    #am2x,am2y,am2z = am2
                    for m in 0:(mmax-a-c)
                        values[(avx,avy,avz,cpx,cpy,cpz,m)] = (qxyz[j]-bxyz[j])*values[(avx,avy,avz,cvx,cvy,cvz,m)]+(wxyz[j]-qxyz[j])*values[(avx,avy,avz,cvx,cvy,cvz,m+1)]
                        if cm[j] >= 0
                            values[(avx,avy,avz,cpx,cpy,cpz,m)] += cv[j]/(2*eta)*(values[(avx,avy,avz,cmx,cmy,cmz,m)]-zeta/ze*values[(avx,avy,avz,cmx,cmy,cmz,m+1)])
                        end
                        if am[j] >= 0 
                            values[(avx,avy,avz,cpx,cpy,cpz,m)] += av[j]/(2*ze)*(values[(amx,amy,amz,cvx,cvy,cvz,m+1)])
                        end                             
                    end
                end
            end
        end 
    end
    return prunem(values)
end

"vdiffs(a) - Compute vector differences (ax,ay-1,az) (ax,ay-2,az) where y is the amax(ax,ay,az)
    return y,am,am2"
vdiffs(a,i) = a-unit(3,i),a-2*unit(3,i)

"unit(n,d) - create a n-dim unit vector in direction d"
function unit(n,d) 
    v = zeros(Int,n)
    v[d] = 1
    return v
end

