using MolecularIntegrals, Profile
bfs = build_basis(ethane,"6-31G")
@profview MolecularIntegrals.all_twoe_ints(bfs)