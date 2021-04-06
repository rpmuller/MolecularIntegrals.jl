export Atom, Molecule

mutable struct Atom
    atno::Int
    x::Float64
    y::Float64
    z::Float64
end

mutable struct Molecule
    atomlist::Array{Atom,1}
end

function push!(mol::Molecule,at::Atom)
    Base.push!(atomlist,at)
end

tobohr(x::Float64) = x/0.52918
function tobohr!(at::Atom)
    at.x /= 0.52918
    at.y /= 0.52918
    at.z /= 0.52918
end
function tobohr!(mol::Molecule)
    for at in mol.atomlist
        tobohr!(at)
    end
end

nuclear_repulsion(a::Atom,b::Atom)= a.atno*b.atno/sqrt(dist2(a.x-b.x,a.y-b.y,a.z-b.z))
function nuclear_repulsion(mol::Molecule)
    nr = 0
    for (i,j) in spairs(nat(mol))
        nr += nuclear_repulsion(mol.atomlist[i],mol.atomlist[j])
    end
    return nr
end

nel(mol::Molecule) = sum([at.atno for at in mol.atomlist])
nat(mol::Molecule) = length(mol.atomlist)

# Other molecule methods to implement
# nocc, nclosed, nopen, nup, ndown, stoich, mass,
# center_of_mass, center!

# Array of symbols, masses

