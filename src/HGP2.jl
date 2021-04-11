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
function vrr2(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    values = Dict()
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
    for ashell in 1:amax
        for ap in shell_indices[ashell]
            apx,apy,apz = ap
            i = argmax(ap) # Choose argmax(ap) as the direction to use for building new terms
            a,am = vdiffs(ap,i) #am = ap-1i, am2 = ap-2i in eq 6-7; i direction of change
            ax,ay,az = a
            amx,amy,amz = am
            for m in 0:(mmax-ashell)
                values[(apx,apy,apz,0,0,0,m)] = (P[i]-A[i])*values[(ax,ay,az,0,0,0,m)]+(W[i]-P[i])*values[(ax,ay,az,0,0,0,m+1)]
                if am[i] >= 0
                    values[(apx,apy,apz,0,0,0,m)] += a[i]/(2*zeta)*(values[(amx,amy,amz,0,0,0,m)]-eta/ze*values[(amx,amy,amz,0,0,0,m+1)])
                end
            end
        end
    end

    # Next build (0,0,0,cx,cy,cz,m)
    # The c-based version of 6a is:
    #   [0,c+1]m = (Qi-Bi)[0,c]m + (Wi-Qi)[0,c]m+1
    #       + ci/2eta ([0,c-1]m - zeta/zeta+eta[0,c-1]m+1)         # eq 6c
    # 
    for cshell in 1:cmax
        for cp in shell_indices[cshell]
            cpx,cpy,cpz = cp
            i = argmax(cp)  # Choose argmax(cp) as the direction to use for building new terms
            c,cm = vdiffs(cp,i)
            cx,cy,cz = c
            cmx,cmy,cmz = cm
            for m in 0:(mmax-cshell)
                values[(0,0,0,cpx,cpy,cpz,m)] = (Q[i]-C[i])*values[(0,0,0,cx,cy,cz,m)]+(W[i]-Q[i])*values[(0,0,0,cx,cy,cz,m+1)]
                if cm[i] >= 0
                    values[(0,0,0,cpx,cpy,cpz,m)] += c[i]/(2*eta)*(values[(0,0,0,cmx,cmy,cmz,m)]-zeta/ze*values[(0,0,0,cmx,cmy,cmz,m+1)])
                end
            end
        end
    end

    # Now build (ax,ay,az,cx,cy,cz,m)
    # The c-based version of 6a is:
    #   [a,c+1]m = (Qj-Bi)[a,c]m + (Wj-Qj)[a,c]m+1
    #       + c_j/2eta ([a,c-1]m - zeta/zeta+eta[a,c-1]m+1)         # eq 6d
    #       + a_j/2(zeta+eta)[a-1,c]m+1
    for ashell in 1:amax
        for a in shell_indices[ashell]
            ax,ay,az = a
            for cshell in 1:cmax
                for cp in shell_indices[cshell]
                    cpx,cpy,cpz = cp
                    j = argmax(cp)  # Choose argmax(cp) as the direction to use for building new terms
                    c,cm = vdiffs(cp,j)
                    cx,cy,cz = c
                    cmx,cmy,cmz = cm

                    am,am2 = vdiffs(a,j) #am = a-1j, am2 = a-2j in eq 6-7
                    amx,amy,amz = am
                    #am2x,am2y,am2z = am2
                    for m in 0:(mmax-ashell-cshell)
                        values[(ax,ay,az,cpx,cpy,cpz,m)] = (Q[j]-C[j])*values[(ax,ay,az,cx,cy,cz,m)]+(W[j]-Q[j])*values[(ax,ay,az,cx,cy,cz,m+1)]
                        if cm[j] >= 0
                            values[(ax,ay,az,cpx,cpy,cpz,m)] += c[j]/(2*eta)*(values[(ax,ay,az,cmx,cmy,cmz,m)]-zeta/ze*values[(ax,ay,az,cmx,cmy,cmz,m+1)])
                        end
                        if am[j] >= 0 
                            values[(ax,ay,az,cpx,cpy,cpz,m)] += a[j]/(2*ze)*(values[(amx,amy,amz,cx,cy,cz,m+1)])
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

