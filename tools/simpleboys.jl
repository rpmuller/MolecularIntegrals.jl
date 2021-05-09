### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 6b21238e-f8de-4c77-9c86-9dbd4691b6f5
begin
	# Imports
	using Plots
end

# ╔═╡ 4e16dbbe-b0fc-11eb-2349-ad7ccd0ff56b
md"""
# Simple questions and simple answers about the Boys function
Rick Muller

A large part of electronic structure theory lies in computing the repulsive
interaction between electrons, and a large part of electron repulsion 
is captured in the Boys function that computes the repulsion between
partially overlapping electronic wave functions.

The Boys function is given by [^boys]

$Fm(T)=\int_0^1 u^{2m}\exp(-Tu^2)du$

and was introduced by Boys in 196?. Since that time, a surprising amount
of work has gone into efficiently computing [^gjp] this function. The
best of these methods [^gjp] divides the integration over $u$ into multiple slices
and approximates the integrand in this region using Taylor or Chebyshev
polynomials. Efficient programs [^libint] to evaluate these terms have detailed
tables with nodes and weights for the interpolation and quadrature schemes.

One says "surprising" because we very much know what electron repulsion
looks like in crude terms: at long distances, the interaction should look like Coulombic $1/r$ repulsion, and at shorter distances the interaction rolls
over to a constant value with zero slope. This behavior was approximated by
Mataga and Nishimoto [^mn57] as

$F_{MN}(r) \sim \frac{1}{a+r}$

which we plot below.

"""

# ╔═╡ 8a736133-17f2-446d-9ab5-4bc46bdd556b
fmn(r,a=1) = 1/(a+r);

# ╔═╡ bf014265-210c-412f-9af4-5a654f557de9
begin
	rs = 0:100
	mnvals = 
	plot(rs,[fmn(r,1) for r in rs],label="fmn(r,a=1)",
		title="Behavior of the Mataga-Nishimoto approximation",
		xlabel="R (au)")
end

# ╔═╡ d8061039-367c-42ee-8a65-a32dca37fc67
md"""
We can use the following simple definition for the Boys function taken from
[^libint]:
"""

# ╔═╡ bc53663d-1f9f-43a6-8140-e8eeac1618bb
"Reference Fm from libint"
function fboys(m,T,eps = 1e-10)
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
end;

# ╔═╡ f6fa9869-b454-458d-95e3-818a9112cc52
md"""
The asymptotic behavior of the Boys function for large $T$ is known
to be

$F_\mathrm{asymp}(T,m) = \sqrt{\frac{\pi}{2}}(2m-1)!!/(2T)^{m+0.5}$

which leads to a slightly more complicated approximation form of

$F_\mathrm{approx}(T)=\sqrt{\frac{\pi}{2}}\frac{1}{a\sqrt{\pi/2}+\sqrt{2T}}$
"""

# ╔═╡ e7454913-e406-4b4a-856f-47103a7bf23d
fasymp(m,T) = sqrt(0.5*pi)*prod(2m-1:-2:1)*(2T)^(-m-0.5);

# ╔═╡ cdc2c6c9-e567-4ee6-a002-8fea660e35fa
function fapprox(T,a=1)
	pre = sqrt(0.5*pi)
	return pre/(a*pre + (2T)^(0.5))
end;

# ╔═╡ b5f42a47-eb6c-4b59-9827-a07b81618018
md"""
This approximation leads to good agreement with the Boys function.
"""

# ╔═╡ 105717a4-f481-479d-bfbb-dc35cc3c0997
begin
	plot(rs,[fboys(0,r) for r in rs],label="fboys(0,r)")
	plot!(rs,[fapprox(r,1) for r in rs],label="fapprox(r,a=1)")
	#title!("Parameter-free approximation of Boys function")
end

# ╔═╡ 0b513de7-cca5-45ca-a4fa-2f2604eb7da9
md"""
We note that although the behavior of `fapprox` plotted above is not exact, 
it is essentially parameter-free when `a=1`, and thus represents a very efficient starting point for further analysis.
"""

# ╔═╡ 4dd0b1e6-61f0-4daa-80f1-625dc793b33a
md"""
As an aside, we certainly don't mean to denigrate the excellent work that has
gone into the development of the approximations methods [^gjp] or their 
efficient implementation [^libint]. These are remarkably accurate and fast
methods that have been derived with great thought and applied repeatedly to
chemical calculations with excellent results. Our intent here is to merely
understand whether a more intuitive approach is possible to this admittedly
solved problem.
"""

# ╔═╡ fdd2691f-f26f-4c26-8cc6-0d8edaa5278a
md"""
Instead of

$F_\mathrm{approx}(T)=\sqrt{\frac{\pi}{2}}\frac{1}{a\sqrt{\pi/2}+\sqrt{2T}}$

we can expand the denominator in a term that will go to zero at large r:

$F_\mathrm{approx2}(T)=\sqrt{\frac{\pi}{2}}\frac{1}{a\exp(-bT)\sqrt{\pi/2}+\sqrt{2T}}$
"""

# ╔═╡ a0a09c9f-c6d2-422f-b32a-8e3f8131b3c3
function fapprox2(T,a=1,b=1)
	pre = sqrt(0.5*pi)
	return pre/(a*pre*exp(-b*T) + (2T)^(0.5))
end;

# ╔═╡ ae9658a3-9241-423d-9daa-d6644e8a206a
begin
	rs2=0:0.1:15
	plot(rs2,[fboys(0,r) for r in rs2],label="fboys(0,r)")
	plot!(rs2,[fapprox2(r,1,1) for r in rs2],label="fapprox2(r,a=1,b=1.5)")
	#title!("Parameter-free approximation of Boys function")
end

# ╔═╡ 26bfe04a-99bd-4089-a9fa-609abd30659b
md"""
This is accurate enough that it's worth just plotting the difference:
"""

# ╔═╡ 4c4c7b38-49a1-4615-a988-7acf5b715d71
begin
	b = 1
	plot(rs2,[fboys(0,r)-fapprox2(r,1,b) for r in rs2],label="Error(r)")
	err = round(maximum([abs(fboys(0,r)-fapprox2(r,1,b)) for r in rs2]),digits=3)
	title!("Max error = $err")
end

# ╔═╡ 8c8adf9b-5eff-4942-9a49-95a4aa1e2002
md"""
With very few parameters (essentially just estimating the `b=1.5` fit) it's possible to get a very good reproduction of the Boys function.

If rather than $exp(-bT)$ we assume there is some function of $T$ that multiplies
the decay of the constant term, what would that function look like?
"""

# ╔═╡ 9428670f-5cf7-4c9e-bb93-dd0747a34182
boysdenom(T) = log(1/fboys(0,T) - sqrt(2T)/sqrt(pi/2));

# ╔═╡ 1ff6a167-aa99-4a64-8621-ba9cee33cdcd
approxdenom(T) = -T-0.5;

# ╔═╡ 0dfe990a-afa9-4c33-bf5d-3a655f35a844
md"It turns out the function looks very much like $\exp(-T-0.5)$, apart from a little cusp."

# ╔═╡ 771d4eeb-72e6-4863-a114-ff203a42455b
begin
	plot(rs2,boysdenom.(rs2),label="Boys denom log decay")
	plot!(rs2,approxdenom.(rs2),label="approx log")
end

# ╔═╡ 997e6018-bb0e-4933-8b13-b904e96fc819
md"Plot this as `fapprox3`"

# ╔═╡ 894135da-1ef9-4170-aa49-1b6307e398eb
function fapprox3(T,a=1)
	pre = sqrt(0.5*pi)
	return pre/(a*pre*exp(-T-0.5) + (2T)^(0.5))
end;

# ╔═╡ 9dfb09bb-22d8-44f2-9c4b-98777e43de64
begin
	plot(rs2,[fboys(0,r) for r in rs2],label="fboys(0,r)")
	plot!(rs2,[fapprox3(r) for r in rs2],label="fapprox3(r)")
end

# ╔═╡ e530c663-afb5-46b8-a645-bd1a4eede764
md"Plot the difference close to the origin:"

# ╔═╡ 262cee5f-cdcb-4505-aee5-07077c815ffb
begin
	r3 = 0:0.1:5
	plot(r3,[fboys(0,r)-fapprox3(r) for r in r3],label="Error approx3")
end

# ╔═╡ 44107e8f-116c-475f-b6d4-2b5aecc2f1e8
md"## References
[^boys]: Boys 1961?
[^gjp]: Gill, Johnson, Pople 1991.
[^libint]: Libint
[^mn57]: K. Nishimoto, N. Mataga. *Z. Phys.*, 12, 335-338 (1957).
"

# ╔═╡ Cell order:
# ╟─6b21238e-f8de-4c77-9c86-9dbd4691b6f5
# ╟─4e16dbbe-b0fc-11eb-2349-ad7ccd0ff56b
# ╠═8a736133-17f2-446d-9ab5-4bc46bdd556b
# ╟─bf014265-210c-412f-9af4-5a654f557de9
# ╟─d8061039-367c-42ee-8a65-a32dca37fc67
# ╠═bc53663d-1f9f-43a6-8140-e8eeac1618bb
# ╟─f6fa9869-b454-458d-95e3-818a9112cc52
# ╠═e7454913-e406-4b4a-856f-47103a7bf23d
# ╠═cdc2c6c9-e567-4ee6-a002-8fea660e35fa
# ╟─b5f42a47-eb6c-4b59-9827-a07b81618018
# ╟─105717a4-f481-479d-bfbb-dc35cc3c0997
# ╟─0b513de7-cca5-45ca-a4fa-2f2604eb7da9
# ╟─4dd0b1e6-61f0-4daa-80f1-625dc793b33a
# ╟─fdd2691f-f26f-4c26-8cc6-0d8edaa5278a
# ╠═a0a09c9f-c6d2-422f-b32a-8e3f8131b3c3
# ╠═ae9658a3-9241-423d-9daa-d6644e8a206a
# ╟─26bfe04a-99bd-4089-a9fa-609abd30659b
# ╠═4c4c7b38-49a1-4615-a988-7acf5b715d71
# ╠═8c8adf9b-5eff-4942-9a49-95a4aa1e2002
# ╠═9428670f-5cf7-4c9e-bb93-dd0747a34182
# ╠═1ff6a167-aa99-4a64-8621-ba9cee33cdcd
# ╟─0dfe990a-afa9-4c33-bf5d-3a655f35a844
# ╠═771d4eeb-72e6-4863-a114-ff203a42455b
# ╟─997e6018-bb0e-4933-8b13-b904e96fc819
# ╠═894135da-1ef9-4170-aa49-1b6307e398eb
# ╠═9dfb09bb-22d8-44f2-9c4b-98777e43de64
# ╠═e530c663-afb5-46b8-a645-bd1a4eede764
# ╠═262cee5f-cdcb-4505-aee5-07077c815ffb
# ╠═44107e8f-116c-475f-b6d4-2b5aecc2f1e8
