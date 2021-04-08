export Atom

mutable struct Atom
    atno::Int
    x::Float64
    y::Float64
    z::Float64
end

const BOHR = 0.52918
tobohr(x::Float64) = x/BOHR
function tobohr!(at::Atom)
    at.x /= BOHR
    at.y /= BOHR
    at.z /= BOHR
end
function tobohr!(mol::Vector{Atom})
    for at in mol
        tobohr!(at)
    end
end

nuclear_repulsion(a::Atom,b::Atom)= a.atno*b.atno/sqrt(dist2(a.x-b.x,a.y-b.y,a.z-b.z))
function nuclear_repulsion(mol::Vector{Atom})
    nr = 0
    for (i,j) in spairs(nat(mol))
        nr += nuclear_repulsion(mol[i],mol[j])
    end
    return nr
end

nel(mol::Vector{Atom}) = sum([at.atno for at in mol])
nat(mol::Vector{Atom}) = length(mol)

# Other molecule methods to implement
# nocc, nclosed, nopen, nup, ndown, stoich, mass,
# center_of_mass, center!

# Array of symbols, masses

