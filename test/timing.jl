using MolecularIntegrals, BenchmarkTools


bfs = build_basis(ethane,"6-31G")
@btime MolecularIntegrals.all_twoe_ints(bfs);
#@time MolecularIntegrals.all_twoe_ints(bfs,MolecularIntegrals.coulomb_hgp);
