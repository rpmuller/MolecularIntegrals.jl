using OffsetArrays
using StaticArrays
using DataStructures
export pgbf, cgbf, contract, addbf!, PGBF, CGBF, build_basis,eri_fetcher, Shell, Basis
export m2ao, shell_indices,nao

"""
    PGBF(expn,xyz,I,J,K,norm)

Create a primitive Gaussian basis function 
    g(x,y,z) = norm * (x-x0)^I (y-y0)^J (z-z0)^K exp(-expn*r^2)
The function parameters `xyz` correspond to [x0,y0,z0].
"""
mutable struct PGBF
    expn::Float64
    xyz::SVector{3,Float64}
    I::Int64
    J::Int64
    K::Int64
    norm::Float64
end

"""
    pgbf(expn,x=0,y=0,z=0,I=0,J=0,K=0,norm=1)

Helper function to create a normalized PGBF with some optional
defaults set.    
"""
function pgbf(expn,x=0.0,y=0.0,z=0.0,I=0,J=0,K=0,norm=1)
    p = PGBF(expn,SA[x,y,z],I,J,K,norm)
    normalize!(p)
    return p
end

function amplitude(bf::PGBF,x,y,z)
    dx,dy,dz = dxyz = bf.xyz-SA[x,y,z]
    r2 = dist2(dxyz)
    return bf.norm*(dx^bf.I)*(dy^bf.J)*(dz^bf.K)*exp(-bf.expn*r2)
end
(bf::PGBF)(x,y,z) = amplitude(bf::PGBF,x,y,z)

function normalize!(pbf::PGBF)
    pbf.norm /= sqrt(overlap(pbf,pbf))
end

"""
    CGBF(xyz,I,J,K,norm,[pbgfs],[coefs])

Create a contracted Gaussian basis function made of 
the functions in [pgbfs] with coefficients [coefs].
Also track the origin `xyz` and powers `I,J,K`.
"""
mutable struct CGBF
    xyz::SVector{3,Float64}
    I::Int64
    J::Int64
    K::Int64
    norm::Float64
    pgbfs::Vector{PGBF}
    coefs::Vector{Float64}
end

"""
    cgbf(expn,x=0,y=0,z=0,I=0,J=0,K=0,norm=1)

Helper function to create a CGBF with optional defaults.
"""
cgbf(x=0.0,y=0.0,z=0.0,I=0,J=0,K=0) = CGBF(SA[x,y,z],I,J,K,1.0,[],[])

amplitude(bf::CGBF,x,y,z) = bf.norm*sum(c*amplitude(pbf,x,y,z) for (c,pbf) in primitives(bf))
(bf::CGBF)(x,y,z) = amplitude(bf::CGBF,x,y,z)

function normalize!(bf::CGBF)
    bf.norm /= sqrt(overlap(bf,bf))
end

primitives(a::CGBF) = zip(a.coefs,a.pgbfs)

function addbf!(cbf::CGBF,expn,coef)
    Base.push!(cbf.pgbfs,pgbf(expn,cbf.xyz...,cbf.I,cbf.J,cbf.K))
    Base.push!(cbf.coefs,coef)
    normalize!(cbf)
end

contract(f,a::CGBF,b::CGBF) = a.norm*b.norm*sum(ca*cb*f(abf,bbf) for (ca,abf) in primitives(a) for (cb,bbf) in primitives(b))
function contract(f,a::CGBF,b::CGBF,c::CGBF,d::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            for (cc,cbf) in primitives(c)
                for (cd,dbf) in primitives(d)
                    s += ca*cb*cc*cd*f(abf,bbf,cbf,dbf)
                end
            end
        end
    end
    return a.norm*b.norm*c.norm*d.norm*s
end

function shells(atoms::Vector{Atom},name="sto3g")
    data = basis_data[name]
    shs = Shell[]
    for atom in atoms
        for (sym,primlist) in data[atom.atno]
            L = lvalue[sym]
            expns = [expn for (expn,coef) in primlist]
            coefs = [coef for (expn,coef) in primlist]
            bfs = CGBF[]
            for (I,J,K) in shell_indices[L]
                cbf = cgbf(atom.xyz...,I,J,K)
                push!(bfs,cbf)
                for (expn,coef) in zip(expns,coefs)
                    addbf!(cbf,expn,coef)
                end
            end
            push!(shs,Shell(atom.xyz,lvalue[sym],expns,coefs,bfs))
        end
    end
    return shs
end

function build_basis(atoms::Vector{Atom},name="sto3g")
    shs = shells(atoms,name)
    bfs = CGBF[]
    ishs = Int[]
    mshs = Int[]
    for (ish,sh) in enumerate(shs)
        for (msh,(I,J,K)) in enumerate(shell_indices[sh.L])
            cbf = cgbf(sh.xyz...,I,J,K)
            push!(bfs,cbf)
            push!(ishs,ish)
            push!(mshs,msh)
            for (expn,coef) in zip(sh.expns,sh.coefs)
                addbf!(cbf,expn,coef)
            end
        end
    end
    return Basis(bfs,shs,ishs,mshs)
end

"""
    Shell(xyz,L,expns,coeffs)

Structure for a basis function shell, containing multiple
CGBFs of different angular momenta.        
"""
mutable struct Shell
    xyz::SVector{3,Float64}
    L::Int
    expns::Vector{Float64}
    coefs::Vector{Float64}
    cgbfs::Vector{CGBF}
end
nbf(sh::Shell) = length(shell_indices(sh.L))
function pgbfs(sh::Shell)
    x,y,z = sh.xyz
    pgbfs = PGBF[]
    for (I,J,K) in shell_indices[sh.L]
        for (expn,coef) in zip(sh.expns,sh.coefs)
            push!(pgbfs,pgbf(expn,x,y,z,I,J,K))
        end
    end
    return pgbfs
end


"""
    Basis(cgbfs,shells,ishell,mshell)

Structure to hold a basis set, with info about shells and other data
"""
mutable struct Basis # subset of AbstractVector{CGBF}?
    cgbfs::Vector{CGBF}
    shells::Vector{Shell}
    ishell::Vector{Int64} # Which shell corresponds to bf i
    mshell::Vector{Int64} # Which m-value (px,dxy, etc.) corresponds to bf i
end
Base.size(b::Basis) = (length(b.cgbfs),)
Base.length(b::Basis) = length(b.cgbfs)
Base.iterate(b::Basis,i::Int) = iterate(b.cgbfs,i)
Base.iterate(b::Basis) = iterate(b.cgbfs)
Base.getindex(b::Basis, i::Int) = b.cgbfs[i]
nbf(b::Basis) = length(b.cgbfs)
nshells(b::Basis) = length(b.shells)

"""
  eri_fetcher

Compute all of the required ijkl terms that go into an ERI
record, and compute the required calls to hrr, and how to 
unpack the results into the record.

eri_fetcher returns a dictionary such that the integral structure
may be formed via:
``` 
    fetcher[ishell,jshell,kshell,lshell] = (ijkl,hi,hj,hk,hl)
    hrrs = hrr(ishell,jshell,kshell,lshell)
    ints[ijkl] = hrrs[hi,hj,hk,hl]
```    
"""
function eri_fetcher(bfs::Basis)
    fetcher = DefaultOrderedDict{NTuple{4, Int}, Vector{NTuple{5, Int}}}(Vector{NTuple{5, Int}}[])
    for (index,ijkl) in enumerate(iiterator(length(bfs)))
        i,j,k,l = ijkl
        li,mi = bfs.ishell[i],bfs.mshell[i]
        lj,mj = bfs.ishell[j],bfs.mshell[j]
        lk,mk = bfs.ishell[k],bfs.mshell[k]
        ll,ml = bfs.ishell[l],bfs.mshell[l]
        push!(fetcher[li,lj,lk,ll],(index,mi,mj,mk,ml))
    end
    return fetcher
end

function eri_fetcher(shs::Vector{Shell})
    fetcher = DefaultOrderedDict{NTuple{4, Int}, Vector{NTuple{5, Int}}}(Vector{NTuple{5, Int}}[])
    nsh = length(shs)
    ijk2lm = make_ijk2lm()

    i = 0
    index = 0
    for (ish,shi) in enumerate(shs)
        shi = shs[ish]
        for ibf in shi.cgbfs
            i += 1
            mi = ijk2lm[ibf.I,ibf.J,ibf.K][2]

            j = 0
            for (jsh,shj) in enumerate(take(shs,ish))
                for jbf in shj.cgbfs
                    j += 1
                    if j > i break end
                    mj = ijk2lm[jbf.I,jbf.J,jbf.K][2]

                    k = 0
                    for (ksh,shk) in enumerate(shs)
                        for kbf in shk.cgbfs
                            k += 1
                            mk = ijk2lm[kbf.I,kbf.J,kbf.K][2]

                            l = 0
                            for (lsh,shl) in enumerate(take(shs,ksh))
                                @inbounds for lbf in shl.cgbfs
                                    l += 1
                                    if l > k break end
                                    if k*l > i*j break end
                                    ml = ijk2lm[lbf.I,lbf.J,lbf.K][2]
                                    index += 1
                                    push!(fetcher[ish,jsh,ksh,lsh],(index,mi,mj,mk,ml))
                                end # lbf
                            end # lsh

                        end # kbf
                    end # ksh

                end # jbf
            end # jsh

        end # ibf
    end # ish
       
    return fetcher
end

# shell_indices map from a shell l-value to the Cartesian version of m-values that are the 
#   powers of the Cartesian Gaussian basis functions.
# 
# If desired, we can also invert shell_indices to map IJK triplets to l,m pairs:
# IJK2lm = Dict(IJK =>(l,m) for l in 0:4 for (m,IJK) in enumerate(shell_indices[l]))
const shell_indices = Dict(
    0 => [MVector(0,0,0)], # 1
    1 => [MVector(1,0,0),MVector(0,1,0),MVector(0,0,1)], # 3
    2 => [MVector(2,0,0),MVector(1,1,0),MVector(1,0,1),MVector(0,2,0),MVector(0,1,1),MVector(0,0,2)],
    3 => [MVector(3,0,0),MVector(2,1,0),MVector(2,0,1),
            MVector(1,2,0),MVector(1,0,2),MVector(1,1,1),
            MVector(0,3,0),MVector(0,2,1),MVector(0,1,2),MVector(0,0,3)], # 10
    4 => [MVector(4,0,0),MVector(3,1,0),MVector(3,0,1),MVector(2,2,0),MVector(2,1,1),MVector(2,0,2),
            MVector(1,3,0),MVector(1,2,1),MVector(1,1,2),MVector(1,0,3),
            MVector(0,4,0),MVector(0,3,1),MVector(0,2,2),MVector(0,1,3),MVector(0,0,4)] # 15
)

llabel = Dict(0=>"s",1=>"p",2=>"d",3=>"f",4=>"g",5=>"h")
lvalue = merge(Dict((v,k) for (k,v) in llabel),Dict((uppercase(v),k) for (k,v) in llabel))

bflabel(bf) = llabel[bf.I+bf.J+bf.K]*bfpow("x",bf.I)*bfpow("y",bf.J)*bfpow("z",bf.K)
function bfpow(s,j) 
	if j == 0
		return ""
	elseif j == 1
		return s
	end
	return "$s$j"
end	

"make_m2ao - Make map between a sequential list of (mx,my,mz) values and ao indices"
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
end
function make_ao2m(lmax=4)
    ao2m = MVector{3, Int64}[]
    for l in 0:lmax
        for ms in shell_indices[l]
            push!(ao2m,ms)
        end
    end
    return ao2m
end

const m2ao = make_m2ao()
const ao2m = make_ao2m()

"make_nao - Number of AOs for system with l shells"
make_nao(l) = sum(length(shell_indices[i]) for i in 0:l)
const nao = OffsetArray([make_nao(l) for l in 0:4],0:4)

# Data functions and tables over basis functions
# todo: keep a constant of maxao or something
# and use it to construct tables
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
end
const shift_index = make_shift_index()

function make_shift_direction()
    n = length(ao2m)
    shift_direction = zeros(Int,n)
    for a in 1:n
        shift_direction[a] = argmax(ao2m[a])
    end
    return shift_direction
end
const shift_direction = make_shift_direction()

function make_shell_number()
    n = length(ao2m)
    shell_number = zeros(Int,n)
    for a in 1:n
        shell_number[a] = sum(ao2m[a])
    end
    return shell_number
end
const shell_number = make_shell_number()

function make_shift_index_plus()
    n = length(ao2m)
    maxsh = shell_number[n]
    shift_index = zeros(Int,n,3)
    for a in 1:n
        m = ao2m[a]
        for i in 1:3
            if sum(m) == maxsh
                shift_index[a,i] = 0
            else
                mm = copy(m)
                mm[i] += 1
                am = m2ao[mm]
                shift_index[a,i] = am
            end
        end
    end
    return shift_index
end
const shift_index_plus = make_shift_index_plus()

function make_ijk2lm(maxl=4)
    ijk2lm = Dict{NTuple{3,Int},NTuple{2,Int}}()
    for l in 0:maxl
        for (m,(i,j,k)) in enumerate(shell_indices[l])
            ijk2lm[i,j,k] = (l,m)
        end
    end
    return ijk2lm
end
