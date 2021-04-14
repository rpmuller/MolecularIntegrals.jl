
export pgbf, cgbf, contract, addbf!, PGBF, CGBF, build_basis

mutable struct PGBF
    expn::Float64
    x::Float64
    y::Float64
    z::Float64
    I::Int64
    J::Int64
    K::Int64
    norm::Float64
end


function pgbf(expn,x=0,y=0,z=0,I=0,J=0,K=0,norm=1)
    p = PGBF(expn,x,y,z,I,J,K,norm)
    normalize!(p)
    return p
end

function amplitude(bf::PGBF,x,y,z)
    dx,dy,dz = x-bf.x,y-bf.y,z-bf.z
    r2 = dist2(dx,dy,dz)
    return bf.norm*(dx^bf.I)*(dy^bf.J)*(dz^bf.K)*exp(-bf.expn*r2)
end
(bf::PGBF)(x,y,z) = amplitude(bf::PGBF,x,y,z)

function normalize!(pbf::PGBF)
    pbf.norm /= sqrt(overlap(pbf,pbf))
end

mutable struct CGBF
    x::Float64
    y::Float64
    z::Float64
    I::Int64
    J::Int64
    K::Int64
    norm::Float64
    pgbfs::Vector{PGBF}
    coefs::Vector{Float64}
end

cgbf(x=0,y=0,z=0,I=0,J=0,K=0) = CGBF(x,y,z,I,J,K,1.0,[],[])

function amplitude(bf::CGBF,x,y,z)
    s = 0
    for (c,pbf) in primitives(bf)
        s += c*amplitude(pbf,x,y,z)
    end
    return bf.norm*s
end
(bf::CGBF)(x,y,z) = amplitude(bf::CGBF,x,y,z)

function normalize!(bf::CGBF)
    bf.norm /= sqrt(overlap(bf,bf))
end

primitives(a::CGBF) = zip(a.coefs,a.pgbfs)

function addbf!(cbf::CGBF,expn,coef)
    Base.push!(cbf.pgbfs,pgbf(expn,cbf.x,cbf.y,cbf.z,cbf.I,cbf.J,cbf.K))
    Base.push!(cbf.coefs,coef)
    normalize!(cbf)
end

function contract(f,a::CGBF,b::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            s += ca*cb*f(abf,bbf)
        end
    end
    return a.norm*b.norm*s
end

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


function build_basis(mol::Vector{Atom},name="sto3g")
    data = basis_data[name]
    bfs = []
    for atom in mol
        for btuple in data[atom.atno]
            sym,primlist = btuple
            for (I,J,K) in sym2power[sym]
                cbf = cgbf(atom.x,atom.y,atom.z,I,J,K)
                push!(bfs,cbf)
                for (expn,coef) in primlist
                    addbf!(cbf,expn,coef)
                end
            end
        end
    end
    return bfs
end

sym2power = Dict(
    'S' => [(0,0,0)],
    'P' => [(1,0,0),(0,1,0),(0,0,1)],
    'D' => [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
    )

# shell_indices map from a shell l-value to the Cartesian version of m-values that are the 
#   powers of the Cartesian Gaussian basis functions. These m-values can be generated by
#   the following function (commented out because it is unused):
# mvalues(l) = [[I,J,K] for K in 0:l for J in 0:l for I in 0:l if I+J+K == l]	
# 
# If desired, we can also invert shell_indices to map IJK triplets to l,m pairs:
# IJK2lm = Dict(IJK =>(l,m) for l in 0:4 for (m,IJK) in enumerate(shell_indices[l]))
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

# Basis function labeling routines. Work with CGBF or PBGF. Maybe assign a type union.
llabel = Dict(0=>"s",1=>"p",2=>"d",3=>"f",4=>"g",5=>"h")
bflabel(bf) = llabel[bf.I+bf.J+bf.K]*bfpow("x",bf.I)*bfpow("y",bf.J)*bfpow("z",bf.K)
function bfpow(s,j) 
	if j == 0
		return ""
	elseif j == 1
		return s
	end
	return "$s$j"
end	

