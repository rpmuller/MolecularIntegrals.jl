using MolecularIntegrals

bfs = build_basis(ethane,"6-31G")
@time MolecularIntegrals.all_twoe_ints(bfs)
#@time MolecularIntegrals.all_twoe_ints(bfs,MolecularIntegrals.coulomb_hgp)
