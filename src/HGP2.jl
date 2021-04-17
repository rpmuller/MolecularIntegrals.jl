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

# Note that we can prune more aggressively than this, since in practice
# we only need to return something like the VRRs from 
# [min(amax,bmax),0|min(cmax,dmax)0] to [amax+bmax,0|cmax+dmax,0] rather 
# than returning all terms starting from [00|00]. However, I don't know
# whether there is a savings here since we need to compute all of these
# anyway.

"prunem - Keep only the dictionary keys with m (last index) = 0"
prunem(d::Dict) = Dict(k[1:end-1] => v for (k,v) in d if k[end] == 0)

"vrr3 - compute and contract the vertical recurrance relations 
                                between shells ash,bsh,csh,dsh
"
function vrr3(ash::Shell,bsh::Shell,csh::Shell,dsh::Shell)
    amax,cmax = ash.L+bsh.L,csh.L,dsh.L
    values = OffsetArray(zeros(Float64,amax+1,amax+1,amax+1,cmax+1,cmax+1,cmax+1),
         0:amax,0:amax,0:amax, 0:cmax,0:cmax,0:cmax) 
    A,B,C,D = ash.xyz,bsh.xyz,csh.xyz,dsh.xyz
    for (aexpn,acoef) in zip(ash.expns,ash.coefs)
        for (bexpn,bcoef) in zip(bsh.expns,bsh.coefs)
            for (cexpn,ccoef) in zip(csh.expns,csh.coefs)
                for (dexpn,dcoef) in zip(dsh.expns,dsh.coefs)
                    values += acoef*bcoef*ccoef*dcoef*vrr3(amax,cmax, aexpn,bexpn,cexpn,dexpn,A,B,C,D)
                end
            end
        end
    end
    return values
end

"vrr3 - version of vrr2 that uses arrays instead of dicts"
function vrr3(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax
    values = OffsetArray(zeros(Float64,amax+1,amax+1,amax+1,cmax+1,cmax+1,cmax+1,mmax+1),
         0:amax,0:amax,0:amax, 0:cmax,0:cmax,0:cmax,0:mmax) 
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
    
    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (0,0,0, 0,0,0, m) using eq 12
    for m in 0:mmax
        values[0,0,0, 0,0,0, m] = Kab*Kcd*Fgamma(m,T)/sqrt(ze)
    end

    # Now generate (ax,ay,az,0,0,0,m) 
    # Eq 6a, with c=0 also:
    #   [a+1,0]m = (Pi-Ai)[a,0]m + (Wi-Pi)[a,0]m+1 
    #        + a_i/2zeta ([a-1,0]m - eta/zeta+eta[a-1,0]m+1)        # eq 6b
    for ashell in 1:amax
        for ap in shell_indices[ashell]
            apx,apy,apz = ap
            i = argmax(ap) # Choose argmax(ap) as the direction to use for building new terms
            a = vdiff(ap,i,-1)
            am = vdiff(ap,i,-2)
            ax,ay,az = a
            amx,amy,amz = am
            for m in 0:(mmax-ashell)
                values[apx,apy,apz,0,0,0,m] = (P[i]-A[i])*values[ax,ay,az,0,0,0,m]+(W[i]-P[i])*values[ax,ay,az,0,0,0,m+1]
                if am[i] >= 0
                    values[apx,apy,apz,0,0,0,m] += a[i]/(2*zeta)*(values[amx,amy,amz,0,0,0,m]-eta/ze*values[amx,amy,amz,0,0,0,m+1])
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
            c = vdiff(cp,i,-1)
            cm = vdiff(cp,i,-2)
            cx,cy,cz = c
            cmx,cmy,cmz = cm
            for m in 0:(mmax-cshell)
                values[0,0,0,cpx,cpy,cpz,m] = (Q[i]-C[i])*values[0,0,0,cx,cy,cz,m]+(W[i]-Q[i])*values[0,0,0,cx,cy,cz,m+1]
                if cm[i] >= 0
                    values[0,0,0,cpx,cpy,cpz,m] += c[i]/(2*eta)*(values[0,0,0,cmx,cmy,cmz,m]-zeta/ze*values[0,0,0,cmx,cmy,cmz,m+1])
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
                    c = vdiff(cp,j,-1)
                    cm = vdiff(cp,j,-2)
                    cx,cy,cz = c
                    cmx,cmy,cmz = cm

                    am = vdiff(a,j,-1)
                    amx,amy,amz = am
                    for m in 0:(mmax-ashell-cshell)
                        values[ax,ay,az,cpx,cpy,cpz,m] = (Q[j]-C[j])*values[ax,ay,az,cx,cy,cz,m]+(W[j]-Q[j])*values[ax,ay,az,cx,cy,cz,m+1]
                        if cm[j] >= 0
                            values[ax,ay,az,cpx,cpy,cpz,m] += c[j]/(2*eta)*(values[ax,ay,az,cmx,cmy,cmz,m]-zeta/ze*values[ax,ay,az,cmx,cmy,cmz,m+1])
                        end
                        if am[j] >= 0 
                            values[ax,ay,az,cpx,cpy,cpz,m] += a[j]/(2*ze)*(values[amx,amy,amz,cx,cy,cz,m+1])
                        end                             
                    end
                end
            end
        end 
    end
    return values[:,:,:,:,:,:,0]
end

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
            a = vdiff(ap,i,-1)
            am = vdiff(ap,i,-2)
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
            c = vdiff(cp,i,-1)
            cm = vdiff(cp,i,-2)
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
                    c = vdiff(cp,j,-1)
                    cm = vdiff(cp,j,-2)
                    cx,cy,cz = c
                    cmx,cmy,cmz = cm

                    am = vdiff(a,j,-1)
                    amx,amy,amz = am
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

"vdiff(a,i,n) - Move vector a by n unit vectors in the i direction"
vdiff(a,i,n) = a+n*unit(3,i)

"unit(n,d) - create a n-dim unit vector in direction d"
function unit(n,d) 
    v = zeros(Int,n)
    v[d] = 1
    return v
end

"""
hrr - compute two electron integrals using the horizontal recurrence relations

The HRRs are given in HGP eq 18: 
    (a(b+1),cd) = ((a+1)b,cd) + (Ai-Bi)(ab,cd)

Note that the HRRs apply to either contracted or primitive integrals.
"""
function hrr2(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

    # Get the relevant vrr terms. 
    # Interesting that the Gaussian exponents are simply a pass-through to vrr.
    vrrs = vrr2(ashell+bshell,cshell+dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D) 

    values = Dict()
    # Put vrr values into the hrr values dictionary
    for (ax,ay,az,cx,cy,cz) in keys(vrrs)
        values[(ax,ay,az,0,0,0,cx,cy,cz,0,0,0)] = vrrs[(ax,ay,az,cx,cy,cz)]
    end

    # First build (ab,c0) from (a0,c0)
    for bs in 1:bshell 
        for bp in shell_indices[bs]
            bpx,bpy,bpz = bp
            j = argmax(bp)
            bx,by,bz = vdiff(bp,j,-1)
            for as in ashell:(ashell+bshell-bs) # -1 guess
                for a in shell_indices[as]
                    ax,ay,az = a
                    apx,apy,apz = vdiff(a,j,1)
                    for cs in 0:(cshell+dshell)
                        for c in shell_indices[cs]
                            cx,cy,cz = c
                            values[(ax,ay,az,bpx,bpy,bpz,cx,cy,cz,0,0,0)] = values[(apx,apy,apz,bx,by,bz,cx,cy,cz,0,0,0)] +
                                (A[j]-B[j])*values[(ax,ay,az,bx,by,bz,cx,cy,cz,0,0,0)]
                        end
                    end
                end
            end
        end
    end
    # now build (ab,cd) from (ab,c0)
    for ds in 1:dshell
        for dp in shell_indices[ds]
            dpx,dpy,dpz = dp
            j = argmax(dp)
            dx,dy,dz = vdiff(dp,j,-1)
            for cs in cshell:(cshell+dshell-ds) 
                for c in shell_indices[cs]
                    cx,cy,cz = c
                    cpx,cpy,cpz = vdiff(c,j,1)
                    for a in shell_indices[ashell]
                        ax,ay,az = a
                        for b in shell_indices[bshell]
                            bx,by,bz = b
                            values[(ax,ay,az,bx,by,bz,cx,cy,cz,dpx,dpy,dpz)] = values[(ax,ay,az,bx,by,bz,cpx,cpy,cpz,dx,dy,dz)] +
                                (C[j]-D[j])*values[(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz)]
                        end
                    end
                end
            end
        end
    end
    # We could also just return the subset of values where (b) = bshell and (d) = dshell.
    # But we've already calculated them
    return values
end

"hrr3 - hrr using arrays rather than dicts"
function hrr3(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

    # Get the relevant vrr terms. 
    # Interesting that the Gaussian exponents are simply a pass-through to vrr.
    vrrs = vrr3(ashell+bshell,cshell+dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D) 

    alim = ashell+bshell
    blim = bshell
    clim = cshell+dshell
    dlim = dshell
    values = OffsetArray(zeros(Float64,alim+1,alim+1,alim+1,blim+1,blim+1,blim+1,clim+1,clim+1,clim+1,dlim+1,dlim+1,dlim+1),
                            0:alim,0:alim,0:alim,0:blim,0:blim,0:blim,0:clim,0:clim,0:clim,0:dlim,0:dlim,0:dlim)
    # Put vrr values into the hrr values dictionary
    values[:,:,:,0,0,0,:,:,:,0,0,0] = vrrs[:,:,:,:,:,:]

    # First build (ab,c0) from (a0,c0)
    for bs in 1:bshell 
        for bp in shell_indices[bs]
            bpx,bpy,bpz = bp
            j = argmax(bp)
            bx,by,bz = vdiff(bp,j,-1)
            for as in ashell:(ashell+bshell-bs) # -1 guess
                for a in shell_indices[as]
                    ax,ay,az = a
                    apx,apy,apz = vdiff(a,j,1)
                    for cs in 0:(cshell+dshell)
                        for c in shell_indices[cs]
                            cx,cy,cz = c
                            values[ax,ay,az,bpx,bpy,bpz,cx,cy,cz,0,0,0] = values[apx,apy,apz,bx,by,bz,cx,cy,cz,0,0,0] +
                                (A[j]-B[j])*values[ax,ay,az,bx,by,bz,cx,cy,cz,0,0,0]
                        end
                    end
                end
            end
        end
    end
    # now build (ab,cd) from (ab,c0)
    for ds in 1:dshell
        for dp in shell_indices[ds]
            dpx,dpy,dpz = dp
            j = argmax(dp)
            dx,dy,dz = vdiff(dp,j,-1)
            for cs in cshell:(cshell+dshell-ds) 
                for c in shell_indices[cs]
                    cx,cy,cz = c
                    cpx,cpy,cpz = vdiff(c,j,1)
                    for a in shell_indices[ashell]
                        ax,ay,az = a
                        for b in shell_indices[bshell]
                            bx,by,bz = b
                            values[ax,ay,az,bx,by,bz,cx,cy,cz,dpx,dpy,dpz] = values[ax,ay,az,bx,by,bz,cpx,cpy,cpz,dx,dy,dz] +
                                (C[j]-D[j])*values[ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz]
                        end
                    end
                end
            end
        end
    end
    # We could also just return the subset of values where (b) = bshell and (d) = dshell.
    # But we've already calculated them
    return values
end