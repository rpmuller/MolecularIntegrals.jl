using LoopVectorization

export vrr,hrr,chrr,cvrr,all_twoe_ints_chrr, vrr_autogen
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
    mmax=amax+cmax+1
    vrrs = zeros(Float64,mmax,nao[amax],nao[cmax])

    A,B,C,D = ash.xyz,bsh.xyz,csh.xyz,dsh.xyz
    for (dexpn,dcoef) in zip(dsh.expns,dsh.coefs), (cexpn,ccoef) in zip(csh.expns,csh.coefs)
        for (bexpn,bcoef) in zip(bsh.expns,bsh.coefs), (aexpn,acoef) in zip(ash.expns,ash.coefs)
            #fill!(vrrs,0) # I don't know if this is really necessary
            vrr!(vrrs, amax,cmax, aexpn,bexpn,cexpn,dexpn,A,B,C,D)
            #vrr_turbo!(vrrs, amax,cmax, aexpn,bexpn,cexpn,dexpn,A,B,C,D)
            for j in 1:nao[cmax], i in 1:nao[amax]
                cvrrs[i,j] += (acoef*bcoef*ccoef*dcoef)*vrrs[1,i,j]
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
    for c in 1:nao[cshell+dshell], a in 1:nao[ashell+bshell]
        hrrs[a,1,c,1] = vrrs[a,c]
    end

    # First build (ab,c0) from (a0,c0)
    for c in 1:nao[cshell+dshell], bplus in 2:nao[bshell]
        j = shift_direction[bplus]
        b = shift_index[bplus,j]
        bs = shell_number[b]
        for a in 1:nao[ashell+bshell-bs-1]
            aplus = shift_index_plus[a,j]
            hrrs[a,bplus,c,1] = hrrs[aplus,b,c,1] + (A[j]-B[j])*hrrs[a,b,c,1]
        end
    end

     # now build (ab,cd) from (ab,c0)
    for dplus in 2:nao[dshell]
       j = shift_direction[dplus]
       d = shift_index[dplus,j]
       ds = shell_number[d]
       for c in 1:nao[cshell+dshell-ds-1] 
        cplus = shift_index_plus[c,j]
           for b in 1:nao[bshell],a in 1:nao[ashell]
                hrrs[a,b,c,dplus] = hrrs[a,b,cplus,d] +(C[j]-D[j])*hrrs[a,b,c,d]
            end
        end
    end
    return hrrs
end


"""
vrr!(vrrs, amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

This version passes in the vrrs array to reuse space if possible.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`amax` and `cmax` are related to the sum of the shell angular
momenta for the `a+b`, and `c+d` shells, respectively.
    
The function returns an `n`x`m` array, where `n` is the number
of aos in the `a+b` shell, and `m` is the number of aos in the
`c+d` shell.
"""
function vrr!(vrrs, amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1

    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    ooz,ooe,ooze = 1/zeta,1/eta,1/ze
    oortze = sqrt(ooze)
    P = (aexpn*A + bexpn*B)*ooz
    Q = (cexpn*C + dexpn*D)*ooe
    W = (zeta*P + eta*Q)*ooze
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2*ooze
    KabKcd_rtze = 2pi*pi*sqrt(pi)*ooze*oortze*exp(-aexpn*bexpn*rab2*ooz-cexpn*dexpn*rcd2*ooe)
    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q

    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (1,1,m) using eq 12
    Tcrit=20.0 # Most code uses a much higher Tcrit (117)
    #boys_array = boys_array_Fgamma(mmax,T)
    for m in 1:mmax
        vrrs[m,1,1] = KabKcd_rtze*Fgamma(m-1,T)
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
        aminus = shift_index[a,i]
        if aminus > 0
            a_i = 0.5*ooz*ao2m[a][i]
            for m in 1:lim
                vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1] + a_i*(vrrs[m,aminus,1]-eta*ooze*vrrs[m+1,aminus,1])
            end
        else
            for m in 1:lim
                vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1]
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
        cminus = shift_index[c,i]
        for a in 1:nao[amax]
            ashell = shell_number[a]    
            lim = mmax-cshell-ashell
            aminus = shift_index[a,i]
            if cminus > 0
                c_i = 0.5*ooe*ao2m[c][i]
                if aminus > 0
                    a_i = 0.5*ooze*ao2m[a][i]
                    for m in 1:lim
                        vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                            a_i*vrrs[m+1,aminus,c] +
                            c_i*(vrrs[m,a,cminus]-zeta*ooze*vrrs[m+1,a,cminus])
                    end
                else
                    for m in 1:lim
                        vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                            c_i*(vrrs[m,a,cminus]-zeta*ooze*vrrs[m+1,a,cminus])
                    end
                end
            elseif aminus > 0
                a_i = 0.5*ooze*ao2m[a][i]
                for m in 1:lim
                    vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                        a_i*vrrs[m+1,aminus,c]
                end
            else
                for m in 1:lim
                    vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c]
                end
            end
        end
    end
    return nothing
end

function vrr_turbo!(vrrs, amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1

    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    ooz,ooe,ooze = 1/zeta,1/eta,1/ze
    oortze = sqrt(ooze)
    P = (aexpn*A + bexpn*B)*ooz
    Q = (cexpn*C + dexpn*D)*ooe
    W = (zeta*P + eta*Q)*ooze
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2*ooze
    KabKcd_rtze = 2pi*pi*sqrt(pi)*ooze*oortze*exp(-aexpn*bexpn*rab2*ooz-cexpn*dexpn*rcd2*ooe)
    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q

    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (1,1,m) using eq 12
    Tcrit=20.0 # Most code uses a much higher Tcrit (117)
    #boys_array = boys_array_Fgamma(mmax,T)
    for m in 1:mmax
        vrrs[m,1,1] = KabKcd_rtze*Fgamma(m-1,T)
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
        aminus = shift_index[a,i]
        if aminus > 0
            a_i = 0.5*ooz*ao2m[a][i]
            @turbo for m in 1:lim
                vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1] + a_i*(vrrs[m,aminus,1]-eta*ooze*vrrs[m+1,aminus,1])
            end
        else
            @turbo for m in 1:lim
                vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1]
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
        cminus = shift_index[c,i]
        for a in 1:nao[amax]
            ashell = shell_number[a]    
            lim = mmax-cshell-ashell
            aminus = shift_index[a,i]
            if cminus > 0
                c_i = 0.5*ooe*ao2m[c][i]
                if aminus > 0
                    a_i = 0.5*ooze*ao2m[a][i]
                    @turbo for m in 1:lim
                        vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                            a_i*vrrs[m+1,aminus,c] +
                            c_i*(vrrs[m,a,cminus]-zeta*ooze*vrrs[m+1,a,cminus])
                    end
                else
                    @turbo for m in 1:lim
                        vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m,a,c] + 
                            c_i*(vrrs[m,a,cminus]-zeta*ooze*vrrs[m+1,a,cminus])
                    end
                end
            elseif aminus > 0
                a_i = 0.5*ooze*ao2m[a][i]
                @turbo for m in 1:lim
                    vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                        a_i*vrrs[m+1,aminus,c]
                end
            else
                @turbo for m in 1:lim
                    vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c]
                end
            end
        end
    end
    return nothing
end

function vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1
    vrrs = zeros(Float64,mmax,nao[amax],nao[cmax])
    vrr!(vrrs, amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    return vrrs
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
    for j in 1:nao[cshell+dshell], i in 1:nao[ashell+bshell]
        hrrs[i,1,j,1] = vrrs[1,i,j]
    end
 
     # First build (ab,c0) from (a0,c0)
    for c in 1:nao[cshell+dshell], bplus in 2:nao[bshell]
        j = shift_direction[bplus]
        b = shift_index[bplus,j]
        bs = shell_number[b]
        for a in 1:nao[ashell+bshell-bs-1] 
            aplus = shift_index_plus[a,j]
            hrrs[a,bplus,c,1] = hrrs[aplus,b,c,1] + (A[j]-B[j])*hrrs[a,b,c,1]
        end
     end

     # now build (ab,cd) from (ab,c0)
    for dplus in 2:nao[dshell]
        j = shift_direction[dplus]
        d = shift_index[dplus,j]
        ds = shell_number[d]
        for c in 1:nao[cshell+dshell-ds-1] 
            cplus = shift_index_plus[c,j]
            for b in 1:nao[bshell], a in 1:nao[ashell]
                hrrs[a,b,c,dplus] = hrrs[a,b,cplus,d] +(C[j]-D[j])*hrrs[a,b,c,d]
            end
        end
    end
    return hrrs
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
    indent = "    "
    asize = nao[amax]
    csize = nao[cmax]
    lines = []
    for m in 1:mmax
        push!(lines,"$(indent)vrrs[1,1,$m] = KabKcd_rtze*Fgamma($(m-1),T)")
    end
    slines = join(lines,"\n")

    palines = amax>0 ? "$(indent)PA = P-A\n$(indent)WP = W-P" : ""
    qclines = cmax>0 ? "$(indent)QC = Q-C\n$(indent)WQ = W-Q" : ""

    lines = []
    for aplus in 2:nao[amax]
        ashell = shell_number[aplus]
        i = shift_direction[aplus]
        a = shift_index[aplus,i]
        lim = mmax-ashell
        aminus = shift_index[a,i]
        for m in 1:lim
            line = ["$(indent)vrrs[$aplus,1,$m] = PA[$i]*vrrs[$a,1,$m] + WP[$i]*vrrs[$a,1,$(m+1)]"]
            if aminus > 0
                a_i = ao2m[a][i]
                #push!(line," +\n$indent$indent") # newline + indent
                push!(line," + ") # append to old line
                push!(line,"$a_i*0.5*ooz*(vrrs[$aminus,1,$m]-eta*ooze*vrrs[$aminus,1,$(m+1)])")
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
                line = ["$(indent)vrrs[$a,$cplus,$m] = QC[$i]*vrrs[$a,$c,$m] + WQ[$i]*vrrs[$a,$c,$(m+1)]"]
                if cminus > 0
                    c_i = ao2m[c][i]
                    #push!(line," +\n$indent$indent") # newline + indent
                    push!(line," + ") # append to old line
                    push!(line,"$c_i*0.5*ooe*(vrrs[$a,$cminus,$m]-zeta*ooze*vrrs[$a,$cminus,$(m+1)])")
                end
                if aminus > 0 
                    a_i = ao2m[a][i]
                    #push!(line," +\n$indent$indent") # newline + indent
                    push!(line," + ") # append to old line
                    push!(line,"$a_i*0.5*ooze*vrrs[$aminus,$c,$(m+1)]")
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
    vrrs = zeros(Float64,$asize,$csize,$mmax)
    
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    ooz,ooe,ooze = 1/zeta,1/eta,1/ze
    oortze = sqrt(ooze)
    P = (aexpn*A + bexpn*B)*ooz
    Q = (cexpn*C + dexpn*D)*ooe
    W = (zeta*P + eta*Q)*ooze
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2*ooze
    KabKcd_rtze = 2pi*pi*sqrt(pi)*ooze*exp(-aexpn*bexpn*rab2*ooz-cexpn*dexpn*rcd2*ooe)*oortze
$palines
$qclines
$slines
$glines
    return vrrs[:,:,1]
end
"""

    return function_template
end

