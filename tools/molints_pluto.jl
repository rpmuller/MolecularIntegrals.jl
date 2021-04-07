### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ b7a5404c-9570-11eb-0562-47e5f0701cd0
using MolecularIntegrals, Profile, ProfileSVG

# ╔═╡ f49f47dd-cad0-49e8-9dd4-04ae57a8fb28
md"""
# MolecularIntegrals.jl

## TODO/Tasks:
- Profile code
- Benchmark timings
- Replace as many hand-written routines with library functions
  - dist2 with hypot^2?
  - [X] gammainc, gser, gcf in nuclear attraction ints with Special Functions gamma_inc
  - Similar replacement of beta functions?
- Treat constants more intelligently, and make them `const`
- Move ProfileSVG, Test, Plots dependencies to tools and tests subdirectories
"""

# ╔═╡ 6f3be0b3-856d-45d8-a870-187f09587195
md"## Profiling and Timing"

# ╔═╡ 3b737810-5ac7-4c9e-9307-6e04763ca3d1
bfs = build_basis(ethane,"6-31G")

# ╔═╡ 11bb6095-c88e-41b1-a569-bc705005db32
length(bfs.bfs)

# ╔═╡ 5d3b0917-97d6-4588-8d0c-921288051428
MolecularIntegrals.all_twoe_ints(bfs)

# ╔═╡ 731fee55-4688-467b-b703-5e0a9d94c9e3
MolecularIntegrals.all_twoe_ints(bfs,MolecularIntegrals.coulomb_hgp)

# ╔═╡ aa89c6af-49d9-41db-94a6-3c87a92a4298
md"The Huzinaga and HGP methods showed 7.1 and 44.1 seconds, respectively. The HGP time is very surprising. This was once relatively fast. I guess it shouldn't be too hard to find efficiencies!"

# ╔═╡ 135dae13-7bba-4936-9006-942d2798e158
@profview MolecularIntegrals.all_twoe_ints(bfs)

# ╔═╡ Cell order:
# ╠═b7a5404c-9570-11eb-0562-47e5f0701cd0
# ╟─f49f47dd-cad0-49e8-9dd4-04ae57a8fb28
# ╟─6f3be0b3-856d-45d8-a870-187f09587195
# ╠═3b737810-5ac7-4c9e-9307-6e04763ca3d1
# ╠═11bb6095-c88e-41b1-a569-bc705005db32
# ╠═5d3b0917-97d6-4588-8d0c-921288051428
# ╠═731fee55-4688-467b-b703-5e0a9d94c9e3
# ╠═aa89c6af-49d9-41db-94a6-3c87a92a4298
# ╠═135dae13-7bba-4936-9006-942d2798e158
