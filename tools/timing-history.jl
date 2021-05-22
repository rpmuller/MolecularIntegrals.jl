### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 1d2793ad-34cb-4704-b139-f55e9af06317
using Plots

# ╔═╡ 3957e56c-2d81-4cb8-9825-963ae0dc49e4
using Dates

# ╔═╡ 1196b181-459a-47e8-b57e-f320b0471a89
begin
	dates = ["4/6/2021","4/24/2021","4/25/2021","4/26/2021","4/27/2021",
		"4/28/2021","4/30/2021","5/1/2021","5/7/2021","5/8/2021","5/22/2021"]
	days = [1,18,19,20,21,22,24,25,31,32,46]
	seconds = [44.1,22,7.83,2.12,2.04,1.11,0.938,0.831,0.631,0.535,0.366]
end

# ╔═╡ 9b3d18a3-a176-4402-874c-45f15d1f6e14
plot(days,seconds,xlabel="Days",ylabel="Seconds",title="Benchmark speed for Ethane/6-31G Integrals",label=false)

# ╔═╡ 59ebedbe-b0c9-439b-bec6-a0aa094e6c84
md"Here's how to do this using dates instead of days. You need to define a 
date format object to parse the dates properly. If you don't specify a year,
it sets it arbitrarily to 0001, so it's less work just to add years."

# ╔═╡ 0cd27fa0-2e58-4754-973a-66d758560dbc
df = DateFormat("m/d/y")

# ╔═╡ b0c6ce0c-5909-4d08-a60a-267f1bf7a2e5
Date("4/1",df);

# ╔═╡ 0e2403f7-02f9-48a7-8954-0534fc90bcae
dates2 = Date.(dates,df);

# ╔═╡ e840276d-787b-41e2-8dff-b4f522307a0f
plot(dates2,seconds,yaxis=:log,xlabel="Dates",ylabel="Seconds",title="Benchmark speed for Ethane/6-31G Integrals",label=false)

# ╔═╡ Cell order:
# ╠═1d2793ad-34cb-4704-b139-f55e9af06317
# ╠═3957e56c-2d81-4cb8-9825-963ae0dc49e4
# ╠═1196b181-459a-47e8-b57e-f320b0471a89
# ╠═9b3d18a3-a176-4402-874c-45f15d1f6e14
# ╟─59ebedbe-b0c9-439b-bec6-a0aa094e6c84
# ╠═0cd27fa0-2e58-4754-973a-66d758560dbc
# ╠═b0c6ce0c-5909-4d08-a60a-267f1bf7a2e5
# ╠═0e2403f7-02f9-48a7-8954-0534fc90bcae
# ╠═e840276d-787b-41e2-8dff-b4f522307a0f
