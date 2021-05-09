### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 54e259b0-23ba-4605-8236-c8ff05889aba
using Plots

# ╔═╡ 9d4250e3-4b5d-4af7-86b7-57da81bda6d4
using BenchmarkTools

# ╔═╡ 7d206253-167e-42e0-9181-3cfe94a83b43
using QuadGK

# ╔═╡ 29c80c8a-afa1-11eb-34ad-43d1456634ee
md"# Approximating the Boys function for electron repulsion integrals
Rick Muller, Sandia National Labs

The Boys function is defined by $Fm(T)=\int_0^1 u^{2m}\exp(-Tu^2)du$ and 
is a key element 
of electron repulsion integrals over s-type functions. Because these integrals are
often used in recurrence relationships, evaluating the Boys function can be a major
time component of the runtime of electronic structure codes. As such, it's
worth figuring out how to approximate this function as efficiently as possible.

I'm going to walk through some approximations in this workbook to try to get
a feeling for the different approximation methods and how accurate and fast they
are.

"

# ╔═╡ 7dd1b7b0-bdc4-4ff0-99f5-3aa29598c4d8
md"## Nearly exact Fm function from libint"

# ╔═╡ cb4f06dd-2134-4ac6-b643-13b672638cb0
md"
The Boys function can be computed in a number of different ways, including using 
the (lower) incomplete $\gamma$ function. 

Here is the version of the Boys function taken from the libint reference function.
"

# ╔═╡ d5b77649-107c-40d2-812c-f61adbd65e88
"Reference Fm from libint"
function fm_ref(m,T,eps = 1e-10)
    denom = (m + 0.5)
    term = exp(-T) /2denom
    old_term = 0.0
    sum = term
    while true
        denom += 1
        old_term = term
        term = old_term * T / denom
        sum += term
        term > sum*eps || old_term < term || break
    end
    return sum
end

# ╔═╡ fe492349-52cd-456b-a44d-f68c545bde01
md"
## Asymptotic approximation for large `T`
One extremely simple approximation is the asymptotic form that holds when `T` is big:"

# ╔═╡ 542f4a7b-959b-4e4b-b531-8375e23902db
factorial2(n::Int64) = prod(n:-2:1)

# ╔═╡ 2a735726-bcc8-483e-bc68-68a9d245c4b7
Fasymp(m,T) = sqrt(pi/2)*factorial2(2m-1)/(2T)^(m+0.5)

# ╔═╡ 95418c51-07f6-4e5c-ae23-16061376bfd1
md"Plotting the functions shows that they are very close, even for relatively small values of `T`."

# ╔═╡ f83eaace-d16c-46ba-adcd-9545e326e54c
begin
	x = 0:50
	scatter(x,fm_ref.(0,x),label="F0(T)")
	scatter!(x,fm_ref.(1,x),label="F1(T)")
	scatter!(x,fm_ref.(2,x),label="F2(T)")
	plot!(x,Fasymp.(0,x),label="F0(T)")
	plot!(x,Fasymp.(1,x),label="F1(T)")
	plot!(x,Fasymp.(2,x),label="F2(T)")
end

# ╔═╡ 1ba537d0-9b98-4cce-b2b8-3dec55fe010b
md"We can see how close by plotting the error between the exact and the asymptotic values. Taking `T=10` as a cutoff value looks like it would be correct to chemical accuracy."

# ╔═╡ ef11111a-bbb1-4da4-8f1a-84a3794f6562
begin
	x1 = 0:30
	plot(x1,Fasymp.(0,x1)-fm_ref.(0,x1),label="Err0(T)", yaxis=:log)
	plot!(x1,Fasymp.(1,x1)-fm_ref.(1,x1),label="Err1(T)")
	plot!(x1,Fasymp.(2,x1)-fm_ref.(2,x1),label="Err2(T)")
end

# ╔═╡ 71cdeb54-3ed4-4766-b381-f73190a92b11
md"In the current (5/9/2021) version of the code I have the cutoff set at `Tcrit=20`,
and values lower than this lead to test cases failing. It's important to note that
most other codes have this value set very high (`Tcrit=117`), but I believe this is
also to insure that the downward recursion can be performed assuming that
$\exp(-T)$ is small."

# ╔═╡ 6821afd0-d2eb-4e84-9ef4-f7512b45387b
md"The asymptotic approximation should be a lot faster. Let's see how much. Here's the time of the nearly exact approximation:"

# ╔═╡ 64b18a29-083e-4e34-945f-c3f190c179fa
@btime fm_ref(1.0,10.0)

# ╔═╡ 42c6c2c7-6489-4c77-a803-261e8a619ae5
md"Can't see it in the notebook, but this is coming in at 200-500 nsec.

Let's time the asymptotic value:"

# ╔═╡ df1383d2-cc02-447c-840a-86d66b48cbab
@btime Fasymp(1,10)

# ╔═╡ 90ca75fe-a933-4cad-8b88-ace5358f3104
md"Comes in at 1-2 nsec. Much faster."

# ╔═╡ 971a389c-421f-479f-9942-f7fc71ae988b
md"
## Naive numerical approximation using `QuadGK`
Most of the fast methods of approximating the Boys function involve clever uses 
of numerical integration. To start out, let's just use a canned numerical integration 
scheme from the `QuadGK` Julia module to set a baseline:"

# ╔═╡ e6c91380-d452-4619-9eef-0812d7e1745b
function ffactory(m,T)
	fm(u) = u^2m*exp(-T*u*u)
end

# ╔═╡ 85a7951e-02cb-46da-8347-35b74f617567
md"We need to pass in a lower `rtol` so that we can integrate faster."

# ╔═╡ 2547855a-4be1-4f0a-bda6-7a415560a741
@btime quadgk(ffactory(1,10),0,1,rtol=1e-4)

# ╔═╡ 6e134f8a-177a-453a-aed6-df0efe9f6059
md"Using `rtol=1e-4` (which actually reports a error bound of 3.8e-7), we can integrate this in ~240 nsec. It's faster, but not dramatically so."

# ╔═╡ Cell order:
# ╟─29c80c8a-afa1-11eb-34ad-43d1456634ee
# ╠═54e259b0-23ba-4605-8236-c8ff05889aba
# ╠═9d4250e3-4b5d-4af7-86b7-57da81bda6d4
# ╠═7d206253-167e-42e0-9181-3cfe94a83b43
# ╟─7dd1b7b0-bdc4-4ff0-99f5-3aa29598c4d8
# ╟─cb4f06dd-2134-4ac6-b643-13b672638cb0
# ╠═d5b77649-107c-40d2-812c-f61adbd65e88
# ╟─fe492349-52cd-456b-a44d-f68c545bde01
# ╠═2a735726-bcc8-483e-bc68-68a9d245c4b7
# ╠═542f4a7b-959b-4e4b-b531-8375e23902db
# ╟─95418c51-07f6-4e5c-ae23-16061376bfd1
# ╠═f83eaace-d16c-46ba-adcd-9545e326e54c
# ╟─1ba537d0-9b98-4cce-b2b8-3dec55fe010b
# ╠═ef11111a-bbb1-4da4-8f1a-84a3794f6562
# ╟─71cdeb54-3ed4-4766-b381-f73190a92b11
# ╟─6821afd0-d2eb-4e84-9ef4-f7512b45387b
# ╟─64b18a29-083e-4e34-945f-c3f190c179fa
# ╟─42c6c2c7-6489-4c77-a803-261e8a619ae5
# ╠═df1383d2-cc02-447c-840a-86d66b48cbab
# ╟─90ca75fe-a933-4cad-8b88-ace5358f3104
# ╟─971a389c-421f-479f-9942-f7fc71ae988b
# ╠═e6c91380-d452-4619-9eef-0812d7e1745b
# ╟─85a7951e-02cb-46da-8347-35b74f617567
# ╠═2547855a-4be1-4f0a-bda6-7a415560a741
# ╟─6e134f8a-177a-453a-aed6-df0efe9f6059
