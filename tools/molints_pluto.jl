### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ cd0b0c44-4457-4ab6-ac00-3519ccdf1388
import Pkg; Pkg.add("Plots")

# ╔═╡ b7a5404c-9570-11eb-0562-47e5f0701cd0
using MolecularIntegrals, Profile, ProfileSVG

# ╔═╡ 39a6c055-e9e6-4736-afe5-e46d51c97139
using MolecularIntegrals:Fgamma

# ╔═╡ 1eb1804e-2894-40b7-86dc-98ec4db2b3d7
using Plots

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

# ╔═╡ b07aa4a1-4372-4bc7-9432-2c8aa29b6485
md"## Fgamma investigations"

# ╔═╡ 20dede5e-5b86-4176-8f4a-52609f73971b
begin
	xs = 0:0.01:2
	plot(xs,Fgamma.(1,xs),label="m=0")
	plot!(xs,Fgamma.(2,xs),label="m=1")
	plot!(xs,Fgamma.(3,xs),label="m=2")
end

# ╔═╡ 6f3be0b3-856d-45d8-a870-187f09587195
md"## Profiling and Timing"

# ╔═╡ 3b737810-5ac7-4c9e-9307-6e04763ca3d1
bfs = build_basis(ethane,"6-31G")

# ╔═╡ 11bb6095-c88e-41b1-a569-bc705005db32
length(bfs)

# ╔═╡ 5d3b0917-97d6-4588-8d0c-921288051428
MolecularIntegrals.all_twoe_ints(bfs)

# ╔═╡ aa89c6af-49d9-41db-94a6-3c87a92a4298
md"The Huzinaga and HGP methods showed 7.1 and 44.1 seconds, respectively. The HGP time is very surprising. This was once relatively fast. I guess it shouldn't be too hard to find efficiencies!"

# ╔═╡ 135dae13-7bba-4936-9006-942d2798e158
@profview MolecularIntegrals.all_twoe_ints(bfs)

# ╔═╡ Cell order:
# ╠═b7a5404c-9570-11eb-0562-47e5f0701cd0
# ╟─f49f47dd-cad0-49e8-9dd4-04ae57a8fb28
# ╠═b07aa4a1-4372-4bc7-9432-2c8aa29b6485
# ╠═39a6c055-e9e6-4736-afe5-e46d51c97139
# ╠═cd0b0c44-4457-4ab6-ac00-3519ccdf1388
# ╠═1eb1804e-2894-40b7-86dc-98ec4db2b3d7
# ╠═20dede5e-5b86-4176-8f4a-52609f73971b
# ╟─6f3be0b3-856d-45d8-a870-187f09587195
# ╠═3b737810-5ac7-4c9e-9307-6e04763ca3d1
# ╠═11bb6095-c88e-41b1-a569-bc705005db32
# ╠═5d3b0917-97d6-4588-8d0c-921288051428
# ╠═aa89c6af-49d9-41db-94a6-3c87a92a4298
# ╠═135dae13-7bba-4936-9006-942d2798e158
