### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 1d2793ad-34cb-4704-b139-f55e9af06317
using Plots

# ╔═╡ 1196b181-459a-47e8-b57e-f320b0471a89
begin
	dates = ["4/6","4/24","4/25","4/26","4/27","4/28","4/30","5/1","5/7","5/8"]
	days = [1,18,19,20,21,22,24,25,31,32]
	seconds = [44.1,22,7.83,2.12,2.04,1.11,0.938,0.831,0.631,0.535]
end

# ╔═╡ 9b3d18a3-a176-4402-874c-45f15d1f6e14
plot(days,seconds,yaxis=:log,xlabel="Days",ylabel="Seconds",title="Benchmark speed for Ethane/6-31G Integrals",label=false)

# ╔═╡ Cell order:
# ╠═1d2793ad-34cb-4704-b139-f55e9af06317
# ╠═1196b181-459a-47e8-b57e-f320b0471a89
# ╠═9b3d18a3-a176-4402-874c-45f15d1f6e14
