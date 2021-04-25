# HGP - A (hopefully) fast, simple implementation of Head-Gordon, Pople's [^1]
# ERI recurrence relations.
export vrr_array,hrr_array,chrr,cvrr, hrr_dict, vrr_widearray
#
# We're going to break this into several steps:
#
# 1. Primitive shell generation [ab,cd]
# 
# Given a shell-level description of a basis set ([ss,ss], [sp,ss], [dd,df], 
# including l [ll,ss]) generate all ERIs for the primitive basis functions in that shell.
#
# This is currently done in the routine `vrr` to generate terms [l0,m0]
# 
# 2. Contraction to (ab,cd)
#
# This is currently done after the primitive `vrr` step, by the `cvrr` call over
# shell variables.
#
# 3. Completion of contracted integrals
#
# `hrr` works to turn either primitive [l0,m0] integrals to [ij,kl] integrals,
# or to turn contracted (l0,m0) integrals to (ij,kl) integrals.
#
# It is most efficient to use chrr on the contracted integrals output
# by the shell version of `cvrr`.
# 
# 4. Integral Array Generation
# 
# After all contracted integrals have been computed, populate an integral record 
# with the relevant terms for use in an electronic structure theory code.
#
# The Utils.jl `iiterator` function loops through these in the correct order.
# The Basis.jl `Basis` struct generates a list of contracted basis functions
# with other data to allow efficient construction of the integral record.
#
#
# 5. Future optimizations
# Gill's work on PRISM suggests [^2][^3] being more flexible about when the basis function
# contraction is performed can reduce the operations count, but this will be simple and 
# likely fast.
# 
# There are many other recurrence relations to consider (MD [ref], LRL [ref], 
# OS [ref], Rys [refs]) but this should be a template for those others. In particular,
# several other schemes for the VRR have been proposed. For now, we'll just stick with the
# HGP version, since that is a well-written and reasonably self-contained paper.
# 
# 6. References
# [^1]: A method for two-electron Gaussian integral and integral derivative
#       evaluation using recurrence relations. Martin Head-Gordon and John
#       A. Pople. JCP, 89 (9), 5777, 1988.
# [^2]: Molecular Integrals Over Gaussian Basis Functions. Peter M. W. Gill. Adv.
#       Q. Chem., 25, 141 (1994).
# [^3]: The Prism Algorithm for Two-Electron Integrals. Peter M. W. Gill and John
#       A. Pople. IJQC, 40, 753 (1991).

"cvrr - compute and contract the vertical recurrence relations 
 between shells ash,bsh,csh,dsh. 
"
function cvrr(ash::Shell,bsh::Shell,csh::Shell,dsh::Shell)
    amax,cmax = ash.L+bsh.L,csh.L+dsh.L
    cvrrs = zeros(Float64,nao(amax),nao(cmax))
    A,B,C,D = ash.xyz,bsh.xyz,csh.xyz,dsh.xyz
    for (aexpn,acoef) in zip(ash.expns,ash.coefs)
        for (bexpn,bcoef) in zip(bsh.expns,bsh.coefs)
            for (cexpn,ccoef) in zip(csh.expns,csh.coefs)
                for (dexpn,dcoef) in zip(dsh.expns,dsh.coefs)
                    cvrrs += acoef*bcoef*ccoef*dcoef*vrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn,A,B,C,D)
                end
            end
        end
    end
    return cvrrs
end

function all_twoe_ints_chrr(bfs)
    fetcher = eri_fetcher(bfs)
    nints = length(collect(MolecularIntegrals.iiterator(length(bfs))))
    ints = zeros(Float64,nints)
    for (ishell,jshell,kshell,lshell) in keys(fetcher)
        hrrs = chrr(bfs.shells[ishell],bfs.shells[jshell],bfs.shells[kshell],bfs.shells[lshell])
        for (ijkl,hi,hj,hk,hl) in fetcher[ishell,jshell,kshell,lshell] 
            ints[ijkl] = hrrs[hi,hj,hk,hl]
        end
    end
    return ints
end        

# Questions:
#  - For ethane/sto3g, why am I getting calls to shells (8,8,11,11), where all the
#       a.m. is on the c/d aos.

function chrr(ash::Shell,bsh::Shell,csh::Shell,dsh::Shell)
    ao2m,m2ao = ao_arrays()
    ashell,bshell,cshell,dshell = ash.L,bsh.L,csh.L,dsh.L
    A,B,C,D = ash.xyz,bsh.xyz,csh.xyz,dsh.xyz
    vrrs = cvrr(ash,bsh,csh,dsh) 
    hrrs = zeros(Float64,nao(ashell+bshell),nao(bshell),nao(cshell+dshell),nao(dshell))
    hrrs[:,1,:,1] = vrrs[:,:]

    # First build (ab,c0) from (a0,c0)
    for bs in 1:bshell 
        for bp in shell_indices[bs]
            bpindex = m2ao[bp]
            j = argmax(bp)
            b = vdiff(bp,j,-1)
            bindex = m2ao[b]
            for as in 0:(ashell+bshell-bs)
                for a in shell_indices[as]
                    aindex = m2ao[a]
                    ap = vdiff(a,j,1)
                    apindex = m2ao[ap]
                    for cs in 0:(cshell+dshell)
                        for c in shell_indices[cs]
                            cindex = m2ao[c]
                            hrrs[aindex,bpindex,cindex,1] = hrrs[apindex,bindex,cindex,1] + 
                                (A[j]-B[j])*hrrs[aindex,bindex,cindex,1]
                        end
                    end
                end
            end
        end
    end
    # now build (ab,cd) from (ab,c0)
    for ds in 1:dshell
        for dp in shell_indices[ds]
            dpindex = m2ao[dp]
            j = argmax(dp)
            d = vdiff(dp,j,-1)
            dindex = m2ao[d]
            for cs in 0:(cshell+dshell-ds) 
                for c in shell_indices[cs]
                    cindex = m2ao[c]
                    cp = vdiff(c,j,1)
                    cpindex = m2ao[cp]
                    for as in 0:ashell
                        for a in shell_indices[as]
                            aindex = m2ao[a]
                            for bs in 0:bshell
                                for b in shell_indices[bs]
                                    bindex = m2ao[b]
                                    hrrs[aindex,bindex,cindex,dpindex] = hrrs[aindex,bindex,cpindex,dindex] +
                                        (C[j]-D[j])*hrrs[aindex,bindex,cindex,dindex]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return hrrs#[1:nao(ashell),:,1:nao(cshell),:]
end

"""
vrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`amax` and `cmax` are related to the sum of the shell angular
momenta for the `a+b`, and `c+d` shells, respectively.
    
The function returns an `n`x`m` array, where `n` is the number
of aos in the `a+b` shell, and `m` is the number of aos in the
`c+d` shell.

The speeds of `vrr_array` and `vrr_widearray` are roughly equivalent,
and both interfaces are being retained for convenience.
"""
function vrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1
    ao2m,m2ao = ao_arrays()
    vrrs = zeros(Float64,nao(amax),nao(cmax),mmax)

    # Try to speed this up using a dispatch table to call
    # the hand-written routines. Doesn't appear to help.
    #=
    dispatch = Dict(
        (0,0) => vrr_ss,
        (0,1) => vrr_sp,
        (1,0) => vrr_ps,
        (1,1) => vrr_pp
    )
    if haskey(dispatch,(amax,cmax))
        return dispatch[(amax,cmax)](aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    end
    =#

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

    # First generate (1,1,m) using eq 12
    for m in 1:(mmax)
        vrrs[1,1, m] = Kab*Kcd*Fgamma(m-1,T)/sqrt(ze)
    end

    # Generate (A,1,m) 
    # Eq 6a, with c=0 also:
    #   [a+1,0]m = (Pi-Ai)[a,1]m + (Wi-Pi)[a,1]m+1 
    #        + a_i/2zeta ([a-1,0]m - eta/zeta+eta[a-1,0]m+1)        # eq 6b

    for ashell in 1:amax
        for ap in shell_indices[ashell]
            apindex = m2ao[ap]
            i = argmax(ap) # Choose argmax(ap) as the direction to use for building new terms
            a = vdiff(ap,i,-1)
            aindex = m2ao[a]
            am = vdiff(ap,i,-2)
            for m in 1:(mmax-ashell) 
                vrrs[apindex,1,m] = (P[i]-A[i])*vrrs[aindex,1,m]+(W[i]-P[i])*vrrs[aindex,1,m+1]
            end
            if am[i] >= 0
                amindex = m2ao[am]
                for m in 1:(mmax-ashell) 
                    vrrs[apindex,1,m] += a[i]/(2*zeta)*(vrrs[amindex,1,m]-eta/ze*vrrs[amindex,1,m+1])
                end
            end
        end
    end

    # Next build (1,C,m)
    # The c-based version of 6a is:
    #   [0,c+1]m = (Qi-Bi)[1,c]m + (Wi-Qi)[1,c]m+1
    #       + ci/2eta ([1,c-1]m - zeta/zeta+eta[1,c-1]m+1)         # eq 6c
    # 
    for cshell in 1:cmax
        for cp in shell_indices[cshell]
            cpindex = m2ao[cp]
            i = argmax(cp)  # Choose argmax(cp) as the direction to use for building new terms
            c = vdiff(cp,i,-1)
            cindex = m2ao[c]
            cm = vdiff(cp,i,-2)
           for m in 1:(mmax-cshell) 
                vrrs[1,cpindex,m] = (Q[i]-C[i])*vrrs[1,cindex,m]+(W[i]-Q[i])*vrrs[1,cindex,m+1]
            end
            if cm[i] >= 0
                cmindex = m2ao[cm]
                for m in 1:(mmax-cshell) 
                    vrrs[1,cpindex,m] += c[i]/(2*eta)*(vrrs[1,cmindex,m]-zeta/ze*vrrs[1,cmindex,m+1])
                end
            end
        end
    end

    # Now build (A,C,m)
    # The c-based version of 6a is:
    #   [a,c+1]m = (Qj-Bi)[a,c]m + (Wj-Qj)[a,c]m+1
    #       + c_j/2eta ([a,c-1]m - zeta/zeta+eta[a,c-1]m+1)         # eq 6d
    #       + a_j/2(zeta+eta)[a-1,c]m+1
    for ashell in 1:amax
        for a in shell_indices[ashell]
            aindex = m2ao[a]
            for cshell in 1:cmax
                for cp in shell_indices[cshell]
                    cpindex = m2ao[cp]
                    j = argmax(cp)  # Choose argmax(cp) as the direction to use for building new terms
                    c = vdiff(cp,j,-1)
                    cindex = m2ao[c]
                    cm = vdiff(cp,j,-2)
                    am = vdiff(a,j,-1)
                    for m in 1:(mmax-ashell-cshell) 
                        vrrs[aindex,cpindex,m] = (Q[j]-C[j])*vrrs[aindex,cindex,m]+(W[j]-Q[j])*vrrs[aindex,cindex,m+1]
                    end
                    if cm[j] >= 0
                        cmindex = m2ao[cm]
                        for m in 1:(mmax-ashell-cshell) 
                            vrrs[aindex,cpindex,m] += c[j]/(2*eta)*(vrrs[aindex,cmindex,m]-zeta/ze*vrrs[aindex,cmindex,m+1])
                        end
                    end
                    if am[j] >= 0 
                        amindex = m2ao[am]
                        for m in 1:(mmax-ashell-cshell) 
                            vrrs[aindex,cpindex,m] += a[j]/(2*ze)*(vrrs[amindex,cindex,m+1])
                        end
                    end                             
                end
            end
        end 
    end
    return vrrs[:,:,1]
end

"""
vrr_widearray(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`amax` and `cmax` are related to the sum of the shell angular
momenta for the `a+b`, and `c+d` shells, respectively.
    
The function returns a six-dimensional array over the possible
powers of the `a+b` and `c+d` shell functions.

The speeds of `vrr_array` and `vrr_widearray` are roughly equivalent,
and both interfaces are being retained for convenience.
"""
function vrr_widearray(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax
    vrrs = OffsetArray(zeros(Float64,amax+1,amax+1,amax+1,cmax+1,cmax+1,cmax+1,mmax+1),
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
        vrrs[0,0,0, 0,0,0, m] = Kab*Kcd*Fgamma(m,T)/sqrt(ze)
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
                vrrs[apx,apy,apz,0,0,0,m] = (P[i]-A[i])*vrrs[ax,ay,az,0,0,0,m]+(W[i]-P[i])*vrrs[ax,ay,az,0,0,0,m+1]
                if am[i] >= 0
                    vrrs[apx,apy,apz,0,0,0,m] += a[i]/(2*zeta)*(vrrs[amx,amy,amz,0,0,0,m]-eta/ze*vrrs[amx,amy,amz,0,0,0,m+1])
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
                vrrs[0,0,0,cpx,cpy,cpz,m] = (Q[i]-C[i])*vrrs[0,0,0,cx,cy,cz,m]+(W[i]-Q[i])*vrrs[0,0,0,cx,cy,cz,m+1]
                if cm[i] >= 0
                    vrrs[0,0,0,cpx,cpy,cpz,m] += c[i]/(2*eta)*(vrrs[0,0,0,cmx,cmy,cmz,m]-zeta/ze*vrrs[0,0,0,cmx,cmy,cmz,m+1])
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
                        vrrs[ax,ay,az,cpx,cpy,cpz,m] = (Q[j]-C[j])*vrrs[ax,ay,az,cx,cy,cz,m]+(W[j]-Q[j])*vrrs[ax,ay,az,cx,cy,cz,m+1]
                        if cm[j] >= 0
                            vrrs[ax,ay,az,cpx,cpy,cpz,m] += c[j]/(2*eta)*(vrrs[ax,ay,az,cmx,cmy,cmz,m]-zeta/ze*vrrs[ax,ay,az,cmx,cmy,cmz,m+1])
                        end
                        if am[j] >= 0 
                            vrrs[ax,ay,az,cpx,cpy,cpz,m] += a[j]/(2*ze)*(vrrs[amx,amy,amz,cx,cy,cz,m+1])
                        end                             
                    end
                end
            end
        end 
    end
    return vrrs[:,:,:,:,:,:,0]
end


"""
hrr_array(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's horizontal recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`ashell`, `bshell`, `cshell` and `dshell` are the shell angular
momenta for the `a`, `b`, `c`, and `d` shells, respectively.
 
The function returns a (k,l,m,n)-dimensional array, where the 
dimensions correspond to the number of aos in the a,b,c,d shells.        
"""
function hrr_array(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    ao2m,m2ao = ao_arrays()
    # Get the relevant vrr terms. 
    vrrs = vrr_array(ashell+bshell,cshell+dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D) 
    hrrs = zeros(Float64,nao(ashell+bshell),nao(bshell),nao(cshell+dshell),nao(dshell))
    hrrs[:,1,:,1] = vrrs[:,:] # This is why this version is fastest

    # First build (ab,c0) from (a0,c0)
    for bs in 1:bshell 
        for bp in shell_indices[bs]
            bpindex = m2ao[bp]
            j = argmax(bp)
            b = vdiff(bp,j,-1)
            bindex = m2ao[b]
            for as in 0:(ashell+bshell-bs)
                for a in shell_indices[as]
                    aindex = m2ao[a]
                    ap = vdiff(a,j,1)
                    apindex = m2ao[ap]
                    for cs in 0:(cshell+dshell)
                        for c in shell_indices[cs]
                            cindex = m2ao[c]
                            hrrs[aindex,bpindex,cindex,1] = hrrs[apindex,bindex,cindex,1] + 
                                (A[j]-B[j])*hrrs[aindex,bindex,cindex,1]
                        end
                    end
                end
            end
        end
    end
    # now build (ab,cd) from (ab,c0)
    for ds in 1:dshell
        for dp in shell_indices[ds]
            dpindex = m2ao[dp]
            j = argmax(dp)
            d = vdiff(dp,j,-1)
            dindex = m2ao[d]
            for cs in 0:(cshell+dshell-ds) 
                for c in shell_indices[cs]
                    cindex = m2ao[c]
                    cp = vdiff(c,j,1)
                    cpindex = m2ao[cp]
                    for as in 0:ashell
                        for a in shell_indices[as]
                            aindex = m2ao[a]
                            for bs in 0:bshell
                                for b in shell_indices[bs]
                                    bindex = m2ao[b]
                                    hrrs[aindex,bindex,cindex,dpindex] = hrrs[aindex,bindex,cpindex,dindex] +
                                        (C[j]-D[j])*hrrs[aindex,bindex,cindex,dindex]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return hrrs[1:nao(ashell),:,1:nao(cshell),:]
end

"""
hrr_dict(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's horizontal recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`ashell`, `bshell`, `cshell` and `dshell` are the shell angular
momenta for the `a`, `b`, `c`, and `d` shells, respectively.
 
The function returns a dictionary containing entries for the
relevant integrals. E.g., `hrrs[ax,ay,az,bx,by,bz,cx,cy,cz,dpx,dpy,dpz]`
contains the integral corresponding to the bfs with powers `ax,ay,az`,
`bx,by,bz`, `cx,cy,cz`, `dx,dy,dz`.

`hrr_dict` is slower than `hrr_array`, but is kept for convenience.        
"""
function hrr_dict(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    vrrs = vrr_widearray(ashell+bshell,cshell+dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D) 
    hrrs = Dict{NTuple{12,Int},Float64}() 

    for as in 0:(ashell+bshell)
        for (ax,ay,az) in shell_indices[as]
            for cs in 0:(cshell+dshell)
                for (cx,cy,cz) in shell_indices[cs]
                    hrrs[ax,ay,az,0,0,0,cx,cy,cz,0,0,0] = vrrs[ax,ay,az,cx,cy,cz]
                end
            end
        end
    end

    # First build (ab,c0) from (a0,c0)
    for bs in 1:bshell 
        for bp in shell_indices[bs]
            bpx,bpy,bpz = bp       
            j = argmax(bp)
            bx,by,bz = vdiff(bp,j,-1)
            for as in 0:(ashell+bshell-bs)
                for a in shell_indices[as]
                    ax,ay,az = a
                    apx,apy,apz = vdiff(a,j,1)
                    for cs in 0:(cshell+dshell)
                        for c in shell_indices[cs]
                            cx,cy,cz = c
                            hrrs[ax,ay,az,bpx,bpy,bpz,cx,cy,cz,0,0,0] = hrrs[apx,apy,apz,bx,by,bz,cx,cy,cz,0,0,0] + 
                                (A[j]-B[j])*hrrs[ax,ay,az,bx,by,bz,cx,cy,cz,0,0,0]
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
            for cs in 0:(cshell+dshell-ds) 
                for c in shell_indices[cs]
                    cx,cy,cz = c
                    cpx,cpy,cpz = vdiff(c,j,1)
                    for as in 0:ashell
                        for a in shell_indices[as]
                            ax,ay,az = a
                            for bs in 0:bshell
                                for b in shell_indices[bs]
                                    bx,by,bz = b
                                    hrrs[ax,ay,az,bx,by,bz,cx,cy,cz,dpx,dpy,dpz] = hrrs[ax,ay,az,bx,by,bz,cpx,cpy,cpz,dx,dy,dz] +
                                        (C[j]-D[j])*hrrs[ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return hrrs
end
