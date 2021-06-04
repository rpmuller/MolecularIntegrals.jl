### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ b7c63d26-0827-491f-908a-02c3c8bf5d0c
begin
	using BenchmarkTools, LinearAlgebra, OffsetArrays, SpecialFunctions, StaticArrays
	using Plots, Profile, ProfileSVG, PlutoUI
	using LoopVectorization
end

# ╔═╡ 7f1163a2-5a4b-4892-a5ed-cac70c41d3b7
md"# Standalone Quantum Chemistry Molecular Integrals"

# ╔═╡ e109b7eb-725d-499c-b136-9552fff9da38
PlutoUI.TableOfContents()

# ╔═╡ 927fecfa-bb26-11eb-3126-6393fad683fd
md"
The following file is the time-intensive part of the [MolecularIntegrals.jl]()
package. This is called thousands of times in a typical calculation, and every
bit of speed is important here.

I'd be grateful for any help people could give me on speeding things up. Thus far,
I have not had any luck speeding this up further using macros like `@simd`
or `@turbo`. Often one of these macros speeds up the timing of the `vrr!` routine, but slows down the overall code for [reasons I don't understand](https://discourse.julialang.org/t/turbo-speeds-routine-slows-down-everything-else/62163).

I've read through and checked everything in the [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/) page and @ChrisRackauckas' [7 Julia Gotchas and How to Handle Them](https://www.stochasticlifestyle.com/7-julia-gotchas-handle/) blog post.

I believe there are still performance gains to be had, because the head to head timing of my code against those C/C++ libraries shows that the Julia code is still 5-10 times slower.

The best algorithms in common use make use of recurrence relations to generate integrals for higher angular momentum basis functions in terms of the lower angular momentum integrals, which are ultimate worked on in terms of incomplete error functions. These methods involve a certain amount of irregular memory access.
"

# ╔═╡ 9282d572-29e6-4c48-8981-4280cf0de817
md"## vrr routine"

# ╔═╡ 94bcfcc1-5216-4a2d-a136-0c8d3dae675a
md"""
Here's the version that has @turbo turned on. It's faster for the tests here, but
for some reason slows down the full 2e integral code.
"""

# ╔═╡ d4aa62ea-7bd3-4094-ac77-1d50a97d3b6d
md"""
Here's a version that has @inbounds turned on. It has similar results to @turbo: ie faster for the single-routine behchmark, but slower on the full 2e integral code.
"""

# ╔═╡ 7fe4e843-7543-45ea-8765-ed04dfab157e
md"## Timing results
The standalone `vrr` routine is coming in at ~28 μs on my development box (old mac mini); the one with `@turbo` clocks at ~19 μs, the one with @inbounds is ~25 μs.
"

# ╔═╡ 1fb56ae0-654d-4288-9da7-1ac4869115a4
md"""## Profiling
Profiling is not working well from Pluto. For comparison, here's a recent screenshot of the full MolecularIntegrals.jl profile from vscode.
![MolecularIntegrals.jl profile](https://raw.githubusercontent.com/rpmuller/MolecularIntegrals.jl/master/profile-2021-05-22.png)
`vrr!` is toward the bottom. `Fgamma` takes up a big chunk, and much of the rest are standard Julia functions like `*` and `Base.math.exp` and `getindex`.
"""

# ╔═╡ 63626a39-616b-4264-932d-114e3e7aca30
ProfileSVG.view(maxdepth=65)

# ╔═╡ 92bf1847-1b11-4a2d-9231-ef629ccbf0f6
md"## Utility functions
"

# ╔═╡ d01acace-1209-425f-826a-5efadb27b947
@inline factorial2(n::Int64) = prod(n:-2:1); # double factorial !!

# ╔═╡ 9174abad-51b6-4dbe-b326-2ee25b56b8af
@inline dist2(dxyz) = dot(dxyz,dxyz);

# ╔═╡ 0b75d079-bdb6-48c2-bb7f-d3d257bc17f6
@inline dist2(xyz1,xyz2) = dist2(xyz1-xyz2);

# ╔═╡ 7111c63d-bbc7-4343-8941-0b87228a9694
"gammainc - return the lower incomplete gamma function"
@inline gammainc(a,x) = gamma(a)*gamma_inc(a,x)[1];

# ╔═╡ b4dde3a6-7ab4-4f4f-8270-ff32d8d3e12b
"boys_array_gamma - Compute the Boys Fm(T) function.
This is currently one of the rate determining steps of the program, but
it's a *separate* problem than the one I want to explore in this notebook."
function boys_array_gamma(mmax,T,SMALL=1e-18)
    T = max(T,SMALL) # needs underflow protection because of inverse
    boys_array = zeros(Float64,mmax)
    ooT = 1/T
    denom = sqrt(ooT)
    boys_array[1] = 0.5*denom*gammainc(0.5,T) 
    for m in 2:mmax
        denom *= ooT
        # Could speed this up more by expressing gamma(m) in terms of gamma(m±1)
        boys_array[m] = 0.5*denom*gammainc(m-0.5,T) 
    end
    return boys_array
end;

# ╔═╡ ce86e563-3dee-4c0c-9cc7-a3efbfc3e326
md"## Global data
I define global data, which I understand is discouraged,
but which are all defined `const`. Timing increases
when I inline these in `vrr`."

# ╔═╡ 828390fb-fe8c-4a8d-8882-642ea78e3555
const shell_indices = Dict(
	# map from a shell l-value to the Cartesian version of m-values that 
	# are the powers of the Cartesian Gaussian basis functions. 
	# Uses `StaticArrays` for some speed gains."
    0 => [MVector(0,0,0)], # 1
    1 => [MVector(1,0,0),MVector(0,1,0),MVector(0,0,1)], # 3
    2 => [MVector(2,0,0),MVector(1,1,0),MVector(1,0,1),MVector(0,2,0),
			MVector(0,1,1),MVector(0,0,2)],
    3 => [MVector(3,0,0),MVector(2,1,0),MVector(2,0,1),
            MVector(1,2,0),MVector(1,0,2),MVector(1,1,1),
            MVector(0,3,0),MVector(0,2,1),MVector(0,1,2),MVector(0,0,3)], # 10
    4 => [MVector(4,0,0),MVector(3,1,0),MVector(3,0,1),MVector(2,2,0),
			MVector(2,1,1),MVector(2,0,2),
            MVector(1,3,0),MVector(1,2,1),MVector(1,1,2),
			MVector(1,0,3),
            MVector(0,4,0),MVector(0,3,1),MVector(0,2,2),
			MVector(0,1,3),MVector(0,0,4)] # 15
);

# ╔═╡ 3b321740-b76a-4b81-896d-a328398e6477
"make_m2ao - Map between a sequential list of `(mx,my,mz)` values and ao indices"
function make_m2ao(lmax=4)
    m2ao = Dict{Vector{Int64}, Int64}()
    iao = 0
    for l in 0:lmax
        for ms in shell_indices[l]
            iao += 1
            m2ao[ms] = iao
        end
    end
    return m2ao
end;

# ╔═╡ f364fe16-17d3-42d0-9ad4-f562168a4a35
function make_ao2m(lmax=4)
    ao2m = MVector{3, Int64}[]
    for l in 0:lmax
        for ms in shell_indices[l]
            push!(ao2m,ms)
        end
    end
    return ao2m
end;

# ╔═╡ f6b8a5d1-18fd-41d3-bd2d-9eb85001af2c
"make_nao - Number of AOs for system with l shells"
make_nao(l) = sum(length(shell_indices[i]) for i in 0:l);

# ╔═╡ c6d84970-9d30-4214-a456-eb03db812aa9
const ao2m = make_ao2m();

# ╔═╡ ec367e77-8eb0-4e66-89a1-249841e1ecfe
function make_shift_direction()
    n = length(ao2m)
    shift_direction = zeros(Int,n)
    for a in 1:n
        shift_direction[a] = argmax(ao2m[a])
    end
    return shift_direction
end;

# ╔═╡ b9142292-e0e8-4bb7-98df-16ab5c23a8c6
const shift_direction = make_shift_direction();

# ╔═╡ 82f11fe4-e95d-4f3d-be43-e65b9c0b3058
function make_shell_number()
    n = length(ao2m)
    shell_number = zeros(Int,n)
    for a in 1:n
        shell_number[a] = sum(ao2m[a])
    end
    return shell_number
end;

# ╔═╡ 58071e7c-13d3-406d-9e93-4ca5e0906b6e
const shell_number = make_shell_number();

# ╔═╡ 474085ee-1a0f-458a-af9a-c746e988dff9
const nao = OffsetArray([make_nao(l) for l in 0:4],0:4);

# ╔═╡ 6a50233d-3948-4a0f-835e-171587e615fc
begin
	# Define variables for the timing functions
	amax = cmax = 4  # Max angular momentum of the Gaussian function
	aexp = bexp = cexp = dexp = 1.0 # exponent of the Gaussian function
	A = [0.0,0.0,0.0] # coordinates of the Gaussian function centers.
	B = [0.0,0.0,1.0]
	C = [0.0,1.0,0.0]
	D = [1.0,0.0,0.0]

	# Allocate main space
    vrrs = zeros(Float64,amax+cmax+1,nao[amax],nao[cmax])
end;

# ╔═╡ 21bf15f1-661e-4859-86d5-12ace608a7b9
const m2ao = make_m2ao();

# ╔═╡ c376bf6b-a153-4efa-9f4a-0e299377a170
function make_shift_index()
    n = length(ao2m)
    shift_index = zeros(Int,n,3)
    for a in 1:n
        m = ao2m[a]
        for i in 1:3
            if m[i] == 0
                shift_index[a,i] = 0
            else
                mm = copy(m)
                mm[i] -= 1
                am = m2ao[mm]
                shift_index[a,i] = am
            end
        end
    end
    return shift_index
end;

# ╔═╡ 21fafb54-d7c7-4dfc-83ad-7142a618ab91
const shift_index = make_shift_index();

# ╔═╡ 52972c5d-2d1c-487f-826e-a9d1196219f7
"vrr!(vrrs, amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`amax` and `cmax` are related to the sum of the shell angular
momenta for the `a+b`, and `c+d` shells, respectively.
    
The function returns an `n`x`m` array, where `n` is the number
of aos in the `a+b` shell, and `m` is the number of aos in the
`c+d` shell.

Equations refer to 
[A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations](https://doi.org/10.1063/1.455553). Martin Head-Gordon and John A. Pople. JCP, 89 (9), 5777, 1988."
function vrr!(vrrs,amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1
	
	# Prefactors arising from multiplying Gaussians together
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    ooz,ooe,ooze = 1/zeta,1/eta,1/ze # "ooz = one over zeta, etc."
    oortze = sqrt(ooze)
    P = (aexpn*A + bexpn*B)*ooz
    Q = (cexpn*C + dexpn*D)*ooe
    W = (zeta*P + eta*Q)*ooze
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2*ooze
    KabKcd_rtze = 2pi*pi*sqrt(pi)*ooze*oortze*exp(
		-aexpn*bexpn*rab2*ooz-cexpn*dexpn*rcd2*ooe)
    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q

    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (1,1,m) using eq 12 using an array of calls to
	# the incomplete gamma function. This is ~40% of the time here,
	# but can be replaced with table lookup or interpolation. The
	# simple version is included here.
    boys_array = boys_array_gamma(mmax,T)
    for m in 1:mmax
        vrrs[m,1,1] = KabKcd_rtze*boys_array[m]
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
                vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1] +
					a_i*(vrrs[m,aminus,1]-eta*ooze*vrrs[m+1,aminus,1])
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
end;

# ╔═╡ b3f833a8-2195-4547-a6d8-4157116b2789
@btime vrr!($vrrs, $amax,$cmax, $aexp,$bexp,$cexp,$dexp, $A,$B,$C,$D)

# ╔═╡ fe06d209-8826-4410-b83e-ccf6e28136bf
@profile (for _ in 1:500; vrr!(vrrs, amax,cmax, aexp,bexp,cexp,dexp, A,B,C,D);end)

# ╔═╡ 7ebda10b-54e4-4def-bb72-831dc9637183
"vrr_turbo!(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`amax` and `cmax` are related to the sum of the shell angular
momenta for the `a+b`, and `c+d` shells, respectively.
    
The function returns an `n`x`m` array, where `n` is the number
of aos in the `a+b` shell, and `m` is the number of aos in the
`c+d` shell.

Equations refer to 
[A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations](https://doi.org/10.1063/1.455553). Martin Head-Gordon and John A. Pople. JCP, 89 (9), 5777, 1988."
function vrr_turbo!(vrrs,amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1
	
	# Prefactors arising from multiplying Gaussians together
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    ooz,ooe,ooze = 1/zeta,1/eta,1/ze # "ooz = one over zeta, etc."
    oortze = sqrt(ooze)
    P = (aexpn*A + bexpn*B)*ooz
    Q = (cexpn*C + dexpn*D)*ooe
    W = (zeta*P + eta*Q)*ooze
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2*ooze
    KabKcd_rtze = 2pi*pi*sqrt(pi)*ooze*oortze*exp(
		-aexpn*bexpn*rab2*ooz-cexpn*dexpn*rcd2*ooe)
    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q

    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (1,1,m) using eq 12 using an array of calls to
	# the incomplete gamma function. This is ~40% of the time here,
	# but can be replaced with table lookup or interpolation. The
	# simple version is included here.
    boys_array = boys_array_gamma(mmax,T)
    @turbo for m in 1:mmax
        vrrs[m,1,1] = KabKcd_rtze*boys_array[m]
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
                vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1] +
					a_i*(vrrs[m,aminus,1]-eta*ooze*vrrs[m+1,aminus,1])
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
                        vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
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
end;

# ╔═╡ d07a14a5-6a5f-4889-bdd3-89e8b134ffdf
@btime vrr_turbo!($vrrs, $amax,$cmax, $aexp,$bexp,$cexp,$dexp, $A,$B,$C,$D)

# ╔═╡ 581b59ec-4d43-4601-bd9a-990a50261603
"vrr_inbounds!(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)

Use Head-Gordon/Pople's vertical recurrence relations to compute
an array of two-electron integrals.

`A`, `B`, `C`, `D` are the centers of four Gaussian functions.
`aexpn`, `bexpn`, `cexpn`, `dexpn` are their exponents.
`amax` and `cmax` are related to the sum of the shell angular
momenta for the `a+b`, and `c+d` shells, respectively.
    
The function returns an `n`x`m` array, where `n` is the number
of aos in the `a+b` shell, and `m` is the number of aos in the
`c+d` shell.

Equations refer to 
[A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations](https://doi.org/10.1063/1.455553). Martin Head-Gordon and John A. Pople. JCP, 89 (9), 5777, 1988."
function vrr_inbounds!(vrrs,amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1
	
	# Prefactors arising from multiplying Gaussians together
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    ze = zeta+eta
    ooz,ooe,ooze = 1/zeta,1/eta,1/ze # "ooz = one over zeta, etc."
    oortze = sqrt(ooze)
    P = (aexpn*A + bexpn*B)*ooz
    Q = (cexpn*C + dexpn*D)*ooe
    W = (zeta*P + eta*Q)*ooze
    rab2 = dist2(A-B)
    rcd2 = dist2(C-D)
    rpq2 = dist2(P-Q)
    T = zeta*eta*rpq2*ooze
    KabKcd_rtze = 2pi*pi*sqrt(pi)*ooze*oortze*exp(
		-aexpn*bexpn*rab2*ooz-cexpn*dexpn*rcd2*ooe)
    PA = P-A
    WP = W-P
    QC = Q-C
    WQ = W-Q

    # HGP equation 6, with b=d=0:
    #   [a+1,c]m = (Pi-Ai)[a,c]m + (Wi-Pi)[a,c]m+1 
    #        + a_i/2zeta ([a-1,c]m - eta/zeta+eta[a-1,c]m+1)        # eq 6a
    #        + ci/2(zeta+eta)[a,c-1]m+1

    # First generate (1,1,m) using eq 12 using an array of calls to
	# the incomplete gamma function. This is ~40% of the time here,
	# but can be replaced with table lookup or interpolation. The
	# simple version is included here.
    boys_array = boys_array_gamma(mmax,T)
    for m in 1:mmax
        @inbounds vrrs[m,1,1] = KabKcd_rtze*boys_array[m]
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
                @inbounds vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1] +
					a_i*(vrrs[m,aminus,1]-eta*ooze*vrrs[m+1,aminus,1])
            end
        else
            for m in 1:lim
                @inbounds vrrs[m,aplus,1] = PA[i]*vrrs[m,a,1] + WP[i]*vrrs[m+1,a,1]
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
                        @inbounds vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                            a_i*vrrs[m+1,aminus,c] +
                            c_i*(vrrs[m,a,cminus]-zeta*ooze*vrrs[m+1,a,cminus])
                    end
                else
                    for m in 1:lim
                        @inbounds vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                            c_i*(vrrs[m,a,cminus]-zeta*ooze*vrrs[m+1,a,cminus])
                    end
                end
            elseif aminus > 0
                a_i = 0.5*ooze*ao2m[a][i]
                for m in 1:lim
                    @inbounds vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c] + 
                        a_i*vrrs[m+1,aminus,c]
                end
            else
                for m in 1:lim
                    @inbounds vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m+1,a,c]
                end
            end
        end
    end
    return nothing
end;

# ╔═╡ ea3639bd-89c5-4b8f-b2d0-066c677824af
@btime vrr_inbounds!($vrrs, $amax,$cmax, $aexp,$bexp,$cexp,$dexp, $A,$B,$C,$D)

# ╔═╡ Cell order:
# ╟─7f1163a2-5a4b-4892-a5ed-cac70c41d3b7
# ╠═b7c63d26-0827-491f-908a-02c3c8bf5d0c
# ╠═e109b7eb-725d-499c-b136-9552fff9da38
# ╟─927fecfa-bb26-11eb-3126-6393fad683fd
# ╟─9282d572-29e6-4c48-8981-4280cf0de817
# ╠═52972c5d-2d1c-487f-826e-a9d1196219f7
# ╟─94bcfcc1-5216-4a2d-a136-0c8d3dae675a
# ╠═7ebda10b-54e4-4def-bb72-831dc9637183
# ╠═d4aa62ea-7bd3-4094-ac77-1d50a97d3b6d
# ╠═581b59ec-4d43-4601-bd9a-990a50261603
# ╠═7fe4e843-7543-45ea-8765-ed04dfab157e
# ╠═6a50233d-3948-4a0f-835e-171587e615fc
# ╠═b3f833a8-2195-4547-a6d8-4157116b2789
# ╠═d07a14a5-6a5f-4889-bdd3-89e8b134ffdf
# ╠═ea3639bd-89c5-4b8f-b2d0-066c677824af
# ╟─1fb56ae0-654d-4288-9da7-1ac4869115a4
# ╠═fe06d209-8826-4410-b83e-ccf6e28136bf
# ╠═63626a39-616b-4264-932d-114e3e7aca30
# ╟─92bf1847-1b11-4a2d-9231-ef629ccbf0f6
# ╠═d01acace-1209-425f-826a-5efadb27b947
# ╠═9174abad-51b6-4dbe-b326-2ee25b56b8af
# ╠═0b75d079-bdb6-48c2-bb7f-d3d257bc17f6
# ╠═b4dde3a6-7ab4-4f4f-8270-ff32d8d3e12b
# ╠═7111c63d-bbc7-4343-8941-0b87228a9694
# ╠═3b321740-b76a-4b81-896d-a328398e6477
# ╠═f364fe16-17d3-42d0-9ad4-f562168a4a35
# ╠═f6b8a5d1-18fd-41d3-bd2d-9eb85001af2c
# ╠═c376bf6b-a153-4efa-9f4a-0e299377a170
# ╠═ec367e77-8eb0-4e66-89a1-249841e1ecfe
# ╠═82f11fe4-e95d-4f3d-be43-e65b9c0b3058
# ╟─ce86e563-3dee-4c0c-9cc7-a3efbfc3e326
# ╠═828390fb-fe8c-4a8d-8882-642ea78e3555
# ╠═58071e7c-13d3-406d-9e93-4ca5e0906b6e
# ╠═b9142292-e0e8-4bb7-98df-16ab5c23a8c6
# ╠═21fafb54-d7c7-4dfc-83ad-7142a618ab91
# ╠═c6d84970-9d30-4214-a456-eb03db812aa9
# ╠═474085ee-1a0f-458a-af9a-c746e988dff9
# ╠═21bf15f1-661e-4859-86d5-12ace608a7b9
