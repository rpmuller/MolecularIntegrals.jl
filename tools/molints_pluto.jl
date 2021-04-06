### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 0633eb4a-5f22-4170-8689-687aaf4bbfea
import Pkg; Pkg.add("ProfileSVG")

# ╔═╡ b7a5404c-9570-11eb-0562-47e5f0701cd0
using MolecularIntegrals, Plots

# ╔═╡ 773d9fcf-6aab-4244-b129-adc284515178
using MolecularIntegrals: gammainc

# ╔═╡ 4e2368da-e0b7-4f58-a3d7-2972c3b010b4
using SpecialFunctions

# ╔═╡ 18fdb77c-1c00-40f0-998b-131e17e55d77
using Profile, ProfileSVG

# ╔═╡ 0f918147-5ead-4c4f-8186-8d626438d6eb
names(MolecularIntegrals)

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

# ╔═╡ 27d0ad39-1402-4dcd-8051-410e23a58be6
md"""
## Using SpecialFunctions
"""

# ╔═╡ 9a98d5f6-b812-4b91-b413-fa831dc2bf99
rng = 0.1:0.1:4

# ╔═╡ 1601f1db-2312-42a2-8855-7f4a06d690b7
begin
	plot(rng,[gammainc(0.5,x) for x in rng],label="gammainc(0+1/2,x)",
	title="PyQuante/NumRec gammainc")
	plot!(rng,[gammainc(1.5,x) for x in rng],label="gammainc(1+1/2,x)")
	plot!(rng,[gammainc(2.5,x) for x in rng],label="gammainc(2+1/2,x)")
	plot!(rng,[gammainc(3.5,x) for x in rng],label="gammainc(3+1/2,x)")
end

# ╔═╡ 2a4774a0-6853-42d2-a925-3254a524a2eb
begin
	plot(rng,[gamma(0,x) for x in rng],label="gamma(0,x)",
		title="Special functions gamma")
	plot!(rng,[gamma(1,x) for x in rng],label="gamma(1,x)")
	plot!(rng,[gamma(2,x) for x in rng],label="gamma(2,x)")
end

# ╔═╡ 9c92d2c4-3538-48c5-b9e8-7f066025acfd
begin
	plot(rng,[1.77*gamma_inc(0.5,x)[1] for x in rng],label="gamma(0,x)")
	plot!(rng,[gammainc(0.5,x) for x in rng],label="gammainc(0+1/2,x)")
end

# ╔═╡ 0d901600-25a5-4b30-95aa-5a7d43bf08d9
md"Okay, they match with the hand-tuned factor of 1.77, but where does this come from? From gamma(1/2)=1.772"

# ╔═╡ 23be8ec4-7f53-4554-a631-f7a303f8f7da
gamma(1/2)

# ╔═╡ 4906f825-a974-4e24-b8a9-c495effc609f
md"Therefore, we can replace the gammainc(a,x) from numrec with the full lower incomplete gamma function, which is computed via:"

# ╔═╡ b233210d-0815-487d-8715-14bcd3fdb0e4
lower_incomplete_gamma(a,x) = gamma(a)*gamma_inc(a,x)[1]

# ╔═╡ 6f3be0b3-856d-45d8-a870-187f09587195
md"## Profiling"

# ╔═╡ d73c0439-b3d2-406e-bebf-a557880eaa3f
begin
	s = pgbf(1.0)
	px = pgbf(1.0,0,0,0,1,0,0)
	
	c = cgbf(0.0,0.0,0.0)
	addbf!(c,1,1)
	
	c2 = cgbf(0,0,0)
	addbf!(c2,1,0.2)
	addbf!(c2,0.5,0.2)
end

# ╔═╡ 656ffed9-fa0e-4615-8610-fa25fd0a35e4
md"This worked momentarily, and now it just shows random reprs:"

# ╔═╡ 135dae13-7bba-4936-9006-942d2798e158
@profview for _ in 1:20; coulomb(s,s,px,px); end

# ╔═╡ be7b289c-a0fd-4a11-a684-40657ee83a2f
md"## Timing"

# ╔═╡ 2b97fb2c-6600-42ce-b91d-03e9b74aa504
md"### HGP terms"

# ╔═╡ 760d1d0b-15be-4a1e-973d-4ec998506d48
function test_vrr()
	ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
	aexpn=bexpn=cexpn=dexpn=1.0
	aI=aJ=aK=0
	cI=cJ=cK=0
	M=0

	for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in [
		(0.,0.,0., 0,0,0, 0,0,0, 4.37335456733),
		(0.,0.,0., 1,0,0, 1,0,0, 0.182223107579),
		(0.,0.,0., 0,1,0, 0,1,0, 0.182223107579),
		(0.,0.,0., 0,0,1, 0,0,1, 0.182223107579),

		(0.,0.,0., 2,0,0, 2,0,0,  0.223223306785),
		(0.,0.,0., 0,2,0, 0,2,0,  0.223223306785),
		(0.,0.,0., 0,0,2, 0,0,2,  0.223223306785),

		(1.,2.,3., 1,0,0, 1,0,0, -5.63387712455e-06),
		(1.,2.,3., 0,1,0, 0,1,0, -0.000116463120359),
		(1.,2.,3., 0,0,1, 0,0,1, -0.000301178525749),

		(1.,2.,3., 2,0,0, 2,0,0, 0.00022503308545040895),
		(1.,2.,3., 0,2,0, 0,2,0, 0.0006102470883881907),
		(1.,2.,3., 0,0,2, 0,0,2, 0.0013427831014563411),

		(0.,0.,0., 1,1,0, 1,1,0, 0.0136667330685),
		(0.,0.,0., 0,1,1, 0,1,1, 0.0136667330685),
		(0.,0.,0., 1,0,1, 1,0,1, 0.0136667330685),

		(3.,2.,1., 1,1,0, 1,1,0, 5.976771621486971e-5),
		(3.,2.,1., 0,1,1, 0,1,1, 1.5742904443905067e-6),
		(3.,2.,1., 1,0,1, 1,0,1, 4.00292848649699e-6)
	]

		val1 = MolecularIntegrals.vrr(
			aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
			cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
		val2 = MolecularIntegrals.vrr(
			cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
			aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
		#println("vrr results: $result $val1 $val2")
		#println("$ax,$ay,$az, $aI,$aJ,$aK, $cI,$cJ,$cK")
	end
end

# ╔═╡ 115dec55-8393-421d-af59-b6b829407fc8
@time test_vrr()

# ╔═╡ Cell order:
# ╠═0633eb4a-5f22-4170-8689-687aaf4bbfea
# ╠═b7a5404c-9570-11eb-0562-47e5f0701cd0
# ╠═0f918147-5ead-4c4f-8186-8d626438d6eb
# ╟─f49f47dd-cad0-49e8-9dd4-04ae57a8fb28
# ╟─27d0ad39-1402-4dcd-8051-410e23a58be6
# ╠═773d9fcf-6aab-4244-b129-adc284515178
# ╠═4e2368da-e0b7-4f58-a3d7-2972c3b010b4
# ╠═9a98d5f6-b812-4b91-b413-fa831dc2bf99
# ╠═1601f1db-2312-42a2-8855-7f4a06d690b7
# ╠═2a4774a0-6853-42d2-a925-3254a524a2eb
# ╠═9c92d2c4-3538-48c5-b9e8-7f066025acfd
# ╠═0d901600-25a5-4b30-95aa-5a7d43bf08d9
# ╠═23be8ec4-7f53-4554-a631-f7a303f8f7da
# ╠═4906f825-a974-4e24-b8a9-c495effc609f
# ╠═b233210d-0815-487d-8715-14bcd3fdb0e4
# ╟─6f3be0b3-856d-45d8-a870-187f09587195
# ╠═18fdb77c-1c00-40f0-998b-131e17e55d77
# ╠═d73c0439-b3d2-406e-bebf-a557880eaa3f
# ╟─656ffed9-fa0e-4615-8610-fa25fd0a35e4
# ╠═135dae13-7bba-4936-9006-942d2798e158
# ╟─be7b289c-a0fd-4a11-a684-40657ee83a2f
# ╟─2b97fb2c-6600-42ce-b91d-03e9b74aa504
# ╠═760d1d0b-15be-4a1e-973d-4ec998506d48
# ╠═115dec55-8393-421d-af59-b6b829407fc8
