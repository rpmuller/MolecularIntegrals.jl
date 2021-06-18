using MolecularIntegrals, Profile
bfs = build_basis(ethane,"cc-pvdz") #"6-31G")
@profview MolecularIntegrals.all_twoe_ints_rys(bfs)