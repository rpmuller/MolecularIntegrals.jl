export vrr,hrr,chrr,cvrr,all_twoe_ints_chrr
# HGP - A (hopefully) fast, simple implementation of Head-Gordon, Pople's [^1]
# ERI recurrence relations.
#
# References
# [^1]: A method for two-electron Gaussian integral and integral derivative
#       evaluation using recurrence relations. Martin Head-Gordon and John
#       A. Pople. JCP, 89 (9), 5777, 1988.
# [^2]: Molecular Integrals Over Gaussian Basis Functions. Peter M. W. Gill. Adv.
#       Q. Chem., 25, 141 (1994).
# [^3]: The Prism Algorithm for Two-Electron Integrals. Peter M. W. Gill and John
#       A. Pople. IJQC, 40, 753 (1991).

"""
    all_twoe_ints_chrr(bfs)

Make multiple calls to `chrr` to form all required two-electron
integrals for a molecular basis set `bfs`.

Returns a 1d array of these integrals [ijkl] with i<j, k<l, and ij<kl. 
"""
function all_twoe_ints_chrr(bfs)
    nbf = length(bfs)
    nints = length(collect(MolecularIntegrals.iiterator(nbf)))
    fetcher = eri_fetcher(bfs)
    ints = zeros(Float64,nints)
    for (ishell,jshell,kshell,lshell) in keys(fetcher)
        hrrs = chrr(bfs.shells[ishell],bfs.shells[jshell],bfs.shells[kshell],bfs.shells[lshell])
        for (ijkl,hi,hj,hk,hl) in fetcher[ishell,jshell,kshell,lshell] 
            ints[ijkl] = hrrs[hi,hj,hk,hl]
        end
    end
    return ints
end        

"""
    cvrr (ash,bsh,csh,dsh)

Compute and contract the vertical recurrence relations 
between shells ash,bsh,csh,dsh. 

This routine takes advantage of the fact that the optimal time for
contracting integrals is after the VRR step.   

This function returns a nao[ash.L+bsh.L] x nao[csh.L+dsh.L]
with the contracted coefficients (a0|b0). It calls the
primitive function vrr multiple times using different
exponents and contraction coefficients.        
"""
function cvrr(ash::Shell,bsh::Shell,csh::Shell,dsh::Shell)
    amax,cmax = ash.L+bsh.L,csh.L+dsh.L
    cvrrs = zeros(Float64,nao[amax],nao[cmax])

    A,B,C,D = ash.xyz,bsh.xyz,csh.xyz,dsh.xyz
    for (aexpn,acoef) in zip(ash.expns,ash.coefs)
        for (bexpn,bcoef) in zip(bsh.expns,bsh.coefs)
            for (cexpn,ccoef) in zip(csh.expns,csh.coefs)
                for (dexpn,dcoef) in zip(dsh.expns,dsh.coefs)
                    cvrrs += acoef*bcoef*ccoef*dcoef*vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn,A,B,C,D)
                end
            end
        end
    end
    return cvrrs
end

"""
    chrr (ash,bsh,csh,dsh)

Compute the contracted horizontal recurrence relations 
between shells ash,bsh,csh,dsh. 

This routine makes calls `cvrr` to form the input VRRs,
and then applies the HRRs to form the full set of two
electron integrals for the four shells ash,bsh,csh,dsh.
    
The routine returns a 4d array where each dimensional
ranges over the number of orbitals in the shell.
"""
function chrr(ash::Shell,bsh::Shell,csh::Shell,dsh::Shell)
    ashell,bshell,cshell,dshell = ash.L,bsh.L,csh.L,dsh.L
    A,B,C,D = ash.xyz,bsh.xyz,csh.xyz,dsh.xyz
    vrrs = cvrr(ash,bsh,csh,dsh) 
    hrrs = zeros(Float64,nao[ashell+bshell],nao[bshell],nao[cshell+dshell],nao[dshell])
    hrrs[:,1,:,1] = vrrs[:,:] 

    # First build (ab,c0) from (a0,c0)
    for bplus in 2:nao[bshell]
        j = shift_direction[bplus]
        b = shift_index[bplus,j]
        bs = shell_number[b]
        for a in 1:nao[ashell+bshell-bs-1]
            aplus = shift_index_plus[a,j]
            for c in 1:nao[cshell+dshell] 
                hrrs[a,bplus,c,1] = hrrs[aplus,b,c,1] + (A[j]-B[j])*hrrs[a,b,c,1]
             end
         end
     end

     # now build (ab,cd) from (ab,c0)
     for dplus in 2:nao[dshell]
        j = shift_direction[dplus]
        d = shift_index[dplus,j]
        ds = shell_number[d]
        for c in 1:nao[cshell+dshell-ds-1] 
            cplus = shift_index_plus[c,j]
            for a in 1:nao[ashell]
                for b in 1:nao[bshell]
                    hrrs[a,b,c,dplus] = hrrs[a,b,cplus,d] +(C[j]-D[j])*hrrs[a,b,c,d]
                 end
             end
         end
     end
     return hrrs[1:nao[ashell],:,1:nao[cshell],:]
 end


"""
vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`amax` and `cmax` are related to the sum of the shell angular
momenta for the `a+b`, and `c+d` shells, respectively.
    
The function returns an `n`x`m` array, where `n` is the number
of aos in the `a+b` shell, and `m` is the number of aos in the
`c+d` shell.
"""
function vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1
    # Removing hand-generated code and retiming:
    #=
    if mmax == 1
        return vrr_ss(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    elseif cmax == 0
        if amax == 1
            return vrr_ps(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        elseif amax == 2
            return vrr_ds(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        end
    elseif cmax == 1
        if amax == 0
            return vrr_sp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        elseif amax == 1
            return vrr_pp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        elseif amax == 2
            return vrr_dp(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        end
    elseif cmax == 2
        if amax == 0
            return vrr_sd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        elseif amax == 1
            return vrr_pd(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        end
    end
    =#
    vrrs = zeros(Float64,nao[amax],nao[cmax],mmax)

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
    for m in 1:mmax
        vrrs[1,1,m] = Kab*Kcd/sqrt(ze)*Fgamma(m-1,T)
    end

    # Generate (A,1,m) 
    # Eq 6a, with c=0 also:
    #   [a+1,0]m = (Pi-Ai)[a,1]m + (Wi-Pi)[a,1]m+1 
    #        + a_i/2zeta ([a-1,0]m - eta/zeta+eta[a-1,0]m+1)        # eq 6b

    for aplus in 2:nao[amax]
        ashell = shell_number[aplus]
        i = shift_direction[aplus]
        a = shift_index[aplus,i]
        lim = mmax-ashell
        for m in 1:lim
            vrrs[aplus,1,m] = (P[i]-A[i])*vrrs[a,1,m] + (W[i]-P[i])*vrrs[a,1,m+1]
        end
        aminus = shift_index[a,i]
        if aminus > 0
            a_i = index_values(a,i)
            for m in 1:lim
                vrrs[aplus,1,m] += a_i/(2*zeta)*(vrrs[aminus,1,m]-eta/ze*vrrs[aminus,1,m+1])
            end
        end
    end

    # Now build (A,C,m)
    # The c-based version of 6a is:
    #   [a,c+1]m = (Qj-Bi)[a,c]m + (Wj-Qj)[a,c]m+1
    #       + c_j/2eta ([a,c-1]m - zeta/zeta+eta[a,c-1]m+1)         # eq 6d
    #       + a_j/2(zeta+eta)[a-1,c]m+1
    for cplus in 2:nao[cmax]
        cshell = shell_number[cplus]
        i = shift_direction[cplus]
        c = shift_index[cplus,i]
        for a in 1:nao[amax]
            ashell = shell_number[a]    
            lim = mmax-cshell-ashell
            for m in 1:lim
                vrrs[a,cplus,m] = (Q[i]-C[i])*vrrs[a,c,m]+(W[i]-Q[i])*vrrs[a,c,m+1]
            end
            cminus = shift_index[c,i]
            if cminus > 0
                c_i = index_values(c,i)
                for m in 1:lim
                    vrrs[a,cplus,m] += c_i/(2*eta)*(vrrs[a,cminus,m]-zeta/ze*vrrs[a,cminus,m+1])
                end
            end
            aminus = shift_index[a,i]
            if aminus > 0 
                a_i = index_values(a,i)
                for m in 1:lim
                    vrrs[a,cplus,m] += a_i/(2*ze)*(vrrs[aminus,c,m+1])
                end
            end                             
        end
    end
    return vrrs[:,:,1]
end

"""
hrr(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's horizontal recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`ashell`, `bshell`, `cshell` and `dshell` are the shell angular
momenta for the `a`, `b`, `c`, and `d` shells, respectively.
 
The function returns a (k,l,m,n)-dimensional array, where the 
dimensions correspond to the number of aos in the a,b,c,d shells.        
"""
function hrr(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    # Get the relevant vrr terms. 
     vrrs = vrr(ashell+bshell,cshell+dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D) 
     hrrs = zeros(Float64,nao[ashell+bshell],nao[bshell],nao[cshell+dshell],nao[dshell])
     hrrs[:,1,:,1] = vrrs[:,:] 
 
     # First build (ab,c0) from (a0,c0)
     for bplus in 2:nao[bshell]
        j = shift_direction[bplus]
        b = shift_index[bplus,j]
        bs = shell_number[b]
        for a in 1:nao[ashell+bshell-bs-1] 
            aplus = shift_index_plus[a,j]
            for c in 1:nao[cshell+dshell] 
                hrrs[a,bplus,c,1] = hrrs[aplus,b,c,1] + (A[j]-B[j])*hrrs[a,b,c,1]
             end
         end
     end

     # now build (ab,cd) from (ab,c0)
     for dplus in 2:nao[dshell]
        j = shift_direction[dplus]
        d = shift_index[dplus,j]
        ds = shell_number[d]
        for c in 1:nao[cshell+dshell-ds-1] 
            cplus = shift_index_plus[c,j]
            for a in 1:nao[ashell]
                for b in 1:nao[bshell]
                    hrrs[a,b,c,dplus] = hrrs[a,b,cplus,d] +(C[j]-D[j])*hrrs[a,b,c,d]
                 end
             end
         end
     end
     return hrrs[1:nao[ashell],:,1:nao[cshell],:]
 end
 
 """
 vrr_autogen(amax,cmax)
 
 Autogenerate the unbranched, unrolled vrr routine between shells
 `amax` and `cmax`.
 """
 function vrr_autogen(amax,cmax)
     mmax = amax+cmax+1
     sa = llabel[amax]
     sc = llabel[cmax]
     indent = "   "
     asize = nao[amax]
     csize = nao[cmax]
     lines = []
     for m in 1:mmax
         push!(lines,"$indent vrrs[1,1,$m] = Kab*Kcd*Fgamma($(m-1),T)/sqrt(ze)")
     end
     slines = join(lines,"\n")
     
     lines = []
     for aplus in 2:nao[amax]
         ashell = shell_number[aplus]
         i = shift_direction[aplus]
         a = shift_index[aplus,i]
         lim = mmax-ashell
         aminus = shift_index[a,i]
         for m in 1:lim
             line = ["$indent vrrs[$aplus,1,$m] = PA[$i]*vrrs[$a,1,$m] + WP[$i]*vrrs[$a,1,$(m+1)]"]
             if aminus > 0
                 a_i = index_values(a,i)
                 push!(line," +\n$indent$indent")
                 push!(line,"$a_i/2zeta*(vrrs[$aminus,1,$m]-eta/ze*vrrs[$aminus,1,$(m+1)])")
             end
             push!(lines,join(line,""))
         end
     end
 
     for cplus in 2:nao[cmax]
         cshell = shell_number[cplus]
         i = shift_direction[cplus]
         c = shift_index[cplus,i]
         for a in 1:nao[amax]
             ashell = shell_number[a]    
             lim = mmax-cshell-ashell
             cminus = shift_index[c,i]
             aminus = shift_index[a,i]
             for m in 1:lim
                 line = ["$indent vrrs[$a,$cplus,$m] = QC[$i]*vrrs[$a,$c,$m] + WQ[$i]*vrrs[$a,$c,$(m+1)]"]
                 if cminus > 0
                     c_i = index_values(c,i)
                     push!(line," +\n$indent$indent")
                     push!(line,"$c_i/2eta*(vrrs[$a,$cminus,$m]-zeta/ze*vrrs[$a,$cminus,$(m+1)])")
                 end
                 if aminus > 0 
                     a_i = index_values(a,i)
                     push!(line," +\n$indent$indent")
                     push!(line,"$a_i/2ze*vrrs[$aminus,$c,$(m+1)]")
                 end
                 push!(lines,join(line,""))
             end                             
         end
     end
     glines = join(lines,"\n")
  
     function_template = """
 "
 vrr_$sa$sc(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
 
 Use Head-Gordon/Pople's vertical recurrence relations to compute
 an array of two-electron integrals.
 
 `A`, `B`, `C`, `D` are the centers of four Gaussian functions.
 `aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
 
 This is an auto generated routine specific to when the A shell 
 has $sa-type angular momentum, and the C shell has $sc-type
 angular momentum.
 "
 function vrr_$sa$sc(aexpn,bexpn,cexpn,dexpn, A,B,C,D)
     mmax = $mmax
     vrrs = zeros(Float64,$asize,$csize,mmax)
     
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
 
 $slines
 
 $glines
 
     return vrrs[:,:,1]
 end
 """
 
     return function_template
 end
 
 #= Works in progress:
 function vrr_shell_recursive(amax,cmax,mmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    P = gaussian_product_center(aexpn,A,bexpn,B)
    Q = gaussian_product_center(cexpn,C,dexpn,D)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    W = gaussian_product_center(zeta,P,eta,Q)

    if amax+cmax == 0
        rab2 = dist2(A-B)
        rcd2 = dist2(C-D)
        rpq2 = dist2(P-Q)
        T = zeta*eta*rpq2/ze
        Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
        Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
        vrrs = zeros(Float64,1,1,mmax)
        for m in 1:mmax
            vrrs[1,1,m] = Kab*Kcd/sqrt(ze)*Fgamma(m-1,T)
        end
        return vrrs
    end

    if cmax > 0
        vrrsm = vrr_shell_recursive(amax,cmax-1,mmax+1, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        # (a,c+dir) = shift(dir)(a,c)
        # if exists(a-dir): (a,c+dir) += (a-dir,c)
        # if exists(c-dir): (a,c+dir) += (a,c-dir)


    return vrrs
 end

 function vrr_shell_dict(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax = amax+cmax+1
    vrrs = Dict{NTuple{Int,2},Matrix{Float64}}

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
    
    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (1,1,m) using eq 12
    vrrs[0,0] = zeros(Float64,1,1,mmax)
    for m in 1:mmax
        vrrs[0,0][1,1,m] = Kab*Kcd/sqrt(ze)*Fgamma(m-1,T)
    end

    for ashell in 1:amax
        na = length(shell_indices[ashell])
        nc = 1
        ma = mmax-ashell
        vrrs[1,0] = zeros(Float64,na,nc,ma)
        ic = nc
        for iap = 1:na
            ap = shell_indices[iap]
            dir = argmax(ap)
            ia = shell_minus(iap,dir)
            for m in 1:ma
                vrrs[ashell,1][iap,ic,m] = PA[dir]*vrrs[ashell-1,0][ia,ic,m]+ WP[dir]vrrs[ashell-1,0][ia,ic,m+1]
            end
            iam = shell_minus(ia,dir)
            if exists(iam)
                for m in 1:ma
                    vrrs[ashell,1][iap,ic,m] += index_values2(ia,ashell-1,dir)*(vrrs[ashell-2,0][iam,ic,m] -eta/ze*vrrs[ashell-2,0][iam,ic,m+1])
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
=#

