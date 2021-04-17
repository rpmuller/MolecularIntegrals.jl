using Printf

export Atom

mutable struct Atom
    atno::Int
    xyz::Vector{Float64}
end

function atom(sym::String,xyz::Vector{Float64},units=:Angstrom)
    atno = sym2no[sym]
    if units == :Angstrom xyz /= 0.52918 end
    return Atom(atno,xyz)
end
function atom(atno::Int,xyz::Vector{Float64},units=:Angstrom)
    if units == :Angstrom xyz /= 0.52918 end
    return Atom(atno,xyz)
end
function atom(atno::Int,x::Float64,y::Float64,z::Float64,units=:Angstrom)
    xyz = [x,y,z]
    if units == :Angstrom xyz /= 0.52918 end
    return Atom(atno,xyz)
end
function atom(sym::String,x::Float64,y::Float64,z::Float64,units=:Angstrom)
    atno = sym2no[sym]
    xyz = [x,y,z]
    if units == :Angstrom xyz /= 0.52918 end
    return Atom(atno,xyz)
end

nuclear_repulsion(a::Atom,b::Atom)= a.atno*b.atno/sqrt(dist2(a.xyz,b.xyz))
function nuclear_repulsion(atoms::Vector{Atom})
    nr = 0
    for (i,j) in spairs(nat(atoms))
        nr += nuclear_repulsion(atoms[i],atoms[j])
    end
    return nr
end

nel(atoms::Vector{Atom}) = sum([at.atno for at in atoms])
nat(atoms::Vector{Atom}) = length(atoms)

# Other molecule methods to implement
# nocc, nclosed, nopen, nup, ndown, stoich, mass,
# center_of_mass, center!

# Array of symbols, masses
symbol = [
    "H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl","Pb","Bi","Po","At","Rn"];
sym2no = Dict((sym => i) for (i,sym) in enumerate(symbol));
mass = [
    1.0008, 4.0026,
    6.941,9.0122,
    10.811,12.011,14.007,15.999,18.998,20.179,
    22.990,24.305,
    26.982,28.086,30.974,32.066,35.453,39.948,
    39.098, 40.078,
    44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845,
    58.9332, 58.6934, 63.546,65.39,
    69.723, 72.61, 74.9216, 78.96, 79.904, 83.80,
    85.4678, 87.62,
    88.90686, 91.224, 92.90638, 95.94, 98, 101.07,
    102.90550, 106.42, 107.8682, 112.411,
    114.818, 118.710, 121.760, 127.60, 126.90447, 131.29,
    132.90545, 137.327, 138.9055, 140.11, 140.90765, 144.24,
    145.0, 150.36, 151.964,
    157.25, 158.92534, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967,
    178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    200.59];  

    function read_xyz(lines)
        nat = parse(Int,strip(lines[1]))
        title = strip(lines[2])
        atoms = []
        for i in 1:nat
            words = split(lines[2+i])
            sym = words[1]
            atno = sym2no[sym]
            xyz = [parse(Float64,w) for w in words[2:4]]
            push!(atoms,Atom(atno,xyz))
        end
        return atoms
    end
    function write_xyz(atoms,fname,title="Written by MolecularIntegrals.jl")
        f = open(fname,"w")
        write(f,"$(length(atoms))\n")
        write(f,"$title\n")
        for atom in atoms
            @printf(f,"%-10s %10.6f %10.6f %10.6f\n",
                symbol[atom.atno],atom.xyz[1],atom.xyz[2],atom.xyz[3])
        end
        close(f)
    end  
    calcmass(atoms) = sum(mass[atom.atno] for atom in atoms)
    com(atoms) = sum(mass[atom.atno]*[atom.x,atom.y,atom.z] for atom in atoms)/calcmass(atoms)
    function translate!(atoms,xyz)
        for atom in atoms
            atom.xyz += xyz
        end
        return nothing
    end 
    # using Bio3DView, able to view("mol.xyz") for a three.js view of molecule.