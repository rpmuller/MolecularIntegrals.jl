
export pgbf, cgbf, contract

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

function push!(cbf::CGBF,expn,coef)
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
