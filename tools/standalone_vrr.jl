### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 57b56aa8-d641-406d-9632-a9073ba69c29
import Pkg; Pkg.add("ProfileSVG")

# ╔═╡ b7c63d26-0827-491f-908a-02c3c8bf5d0c
begin
	using BenchmarkTools, LinearAlgebra, OffsetArrays, SpecialFunctions, StaticArrays
	using Plots, Profile, ProfileSVG
end

# ╔═╡ 927fecfa-bb26-11eb-3126-6393fad683fd
md"# Standalone Quantum Chemistry Molecular Integrals
The following file is the time-intensive part of the [MolecularIntegrals.jl]()
package, to make it easier to get help from other people about improving 
the performance of the code.
"

# ╔═╡ 281b3f52-d3fc-44f7-89bc-ec62b480ce32
md"
`vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)`

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
[A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations](https://doi.org/10.1063/1.455553). Martin Head-Gordon and John A. Pople. JCP, 89 (9), 5777, 1988.
"

# ╔═╡ 92bf1847-1b11-4a2d-9231-ef629ccbf0f6
md"Utility functions"

# ╔═╡ d01acace-1209-425f-826a-5efadb27b947
factorial2(n::Int64) = prod(n:-2:1); # double factorial !!

# ╔═╡ bc891c80-ac5c-484f-afbe-cf2752d99251
dist2(dx,dy,dz) = dx*dx+dy*dy+dz*dz; # TODO: use hypot()^2?

# ╔═╡ 9174abad-51b6-4dbe-b326-2ee25b56b8af
dist2(dxyz) = dot(dxyz,dxyz);

# ╔═╡ 0b75d079-bdb6-48c2-bb7f-d3d257bc17f6
dist2(xyz1,xyz2) = dist2(xyz1-xyz2);

# ╔═╡ 17a919da-822c-4d96-9753-46b7f241ef72
md"boys\_array\_gamma - Compute the Boys Fm(T) function"

# ╔═╡ 079d906b-3172-46be-b515-a9679682a42e
md"`gammainc` - return the lower incomplete gamma function"

# ╔═╡ 7111c63d-bbc7-4343-8941-0b87228a9694
@inline gammainc(a,x) = gamma(a)*gamma_inc(a,x)[1];

# ╔═╡ b4dde3a6-7ab4-4f4f-8270-ff32d8d3e12b
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

# ╔═╡ 8539ea6c-3257-4728-948f-9b354739a849
md"`shell_indices` map from a shell l-value to the Cartesian version of m-values that are the powers of the Cartesian Gaussian basis functions."

# ╔═╡ 828390fb-fe8c-4a8d-8882-642ea78e3555
const shell_indices = Dict(
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

# ╔═╡ cc71e7e8-1d0e-4c42-9312-2ee048738557
md"`make_m2ao` - Map between a sequential list of `(mx,my,mz)` values and ao indices"

# ╔═╡ 3b321740-b76a-4b81-896d-a328398e6477
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

# ╔═╡ 9430e278-c67b-43b7-b471-eca65a8c7a6c
const m2ao = make_m2ao();

# ╔═╡ 8f9d377b-416d-4567-a035-d3504fe25c89
const ao2m = make_ao2m();

# ╔═╡ f6b8a5d1-18fd-41d3-bd2d-9eb85001af2c
"make_nao - Number of AOs for system with l shells"
make_nao(l) = sum(length(shell_indices[i]) for i in 0:l);

# ╔═╡ ef214991-956b-40f4-93df-68c95090e965
const nao = OffsetArray([make_nao(l) for l in 0:4],0:4);

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

# ╔═╡ e0c579d2-61d7-4263-8e34-9d538abe927a
const shift_index = make_shift_index();

# ╔═╡ ec367e77-8eb0-4e66-89a1-249841e1ecfe
function make_shift_direction()
    n = length(ao2m)
    shift_direction = zeros(Int,n)
    for a in 1:n
        shift_direction[a] = argmax(ao2m[a])
    end
    return shift_direction
end;

# ╔═╡ 4b1e88e0-f78b-460a-ba2c-2677472a5dec
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

# ╔═╡ 4ca93908-c10b-4b0f-87c0-59338d14212f
const shell_number = make_shell_number();

# ╔═╡ 7ebda10b-54e4-4def-bb72-831dc9637183
function vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
    mmax=amax+cmax+1
    vrrs = zeros(Float64,mmax,nao[amax],nao[cmax])

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

    # First generate (1,1,m) using eq 12
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
                        vrrs[m,a,cplus] = QC[i]*vrrs[m,a,c]+WQ[i]*vrrs[m,a,c] + 
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
    return vrrs[1,:,:] # We only need the higher m values to compute intermediate terms, so don't return them
end;

# ╔═╡ 6a50233d-3948-4a0f-835e-171587e615fc
begin
	amax = cmax = 4
	aexp = bexp = cexp = dexp = 1.0
	A = [0.0,0.0,0.0]
	B = [0.0,0.0,1.0]
	C = [0.0,1.0,0.0]
	D = [1.0,0.0,0.0]
end

# ╔═╡ b3f833a8-2195-4547-a6d8-4157116b2789
@btime vrr(amax,cmax, aexp,bexp,cexp,dexp, A,B,C,D)

# ╔═╡ fe06d209-8826-4410-b83e-ccf6e28136bf
@profview [vrr(amax,cmax, aexp,bexp,cexp,dexp, A,B,C,D) for _ in 1:200]

# ╔═╡ Cell order:
# ╟─57b56aa8-d641-406d-9632-a9073ba69c29
# ╠═b7c63d26-0827-491f-908a-02c3c8bf5d0c
# ╟─927fecfa-bb26-11eb-3126-6393fad683fd
# ╟─281b3f52-d3fc-44f7-89bc-ec62b480ce32
# ╠═7ebda10b-54e4-4def-bb72-831dc9637183
# ╟─92bf1847-1b11-4a2d-9231-ef629ccbf0f6
# ╠═d01acace-1209-425f-826a-5efadb27b947
# ╠═bc891c80-ac5c-484f-afbe-cf2752d99251
# ╠═9174abad-51b6-4dbe-b326-2ee25b56b8af
# ╠═0b75d079-bdb6-48c2-bb7f-d3d257bc17f6
# ╟─17a919da-822c-4d96-9753-46b7f241ef72
# ╠═b4dde3a6-7ab4-4f4f-8270-ff32d8d3e12b
# ╠═079d906b-3172-46be-b515-a9679682a42e
# ╠═7111c63d-bbc7-4343-8941-0b87228a9694
# ╟─8539ea6c-3257-4728-948f-9b354739a849
# ╠═828390fb-fe8c-4a8d-8882-642ea78e3555
# ╠═cc71e7e8-1d0e-4c42-9312-2ee048738557
# ╠═3b321740-b76a-4b81-896d-a328398e6477
# ╠═f364fe16-17d3-42d0-9ad4-f562168a4a35
# ╠═9430e278-c67b-43b7-b471-eca65a8c7a6c
# ╠═8f9d377b-416d-4567-a035-d3504fe25c89
# ╠═f6b8a5d1-18fd-41d3-bd2d-9eb85001af2c
# ╠═ef214991-956b-40f4-93df-68c95090e965
# ╠═c376bf6b-a153-4efa-9f4a-0e299377a170
# ╠═e0c579d2-61d7-4263-8e34-9d538abe927a
# ╠═ec367e77-8eb0-4e66-89a1-249841e1ecfe
# ╠═4b1e88e0-f78b-460a-ba2c-2677472a5dec
# ╠═82f11fe4-e95d-4f3d-be43-e65b9c0b3058
# ╠═4ca93908-c10b-4b0f-87c0-59338d14212f
# ╠═6a50233d-3948-4a0f-835e-171587e615fc
# ╠═b3f833a8-2195-4547-a6d8-4157116b2789
# ╠═fe06d209-8826-4410-b83e-ccf6e28136bf
