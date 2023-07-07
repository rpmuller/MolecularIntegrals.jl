```julia
using Plots
using BenchmarkTools
using QuadGK
```
# Approximating the Boys function for electron repulsion integrals
Rick Muller, Sandia National Labs

The Boys function is defined by $Fm(T)=\int_0^1 u^{2m}\exp(-Tu^2)du$
and is a key element of electron repulsion integrals over s-type
functions. Because these integrals are often used in recurrence
relationships, evaluating the Boys function can be a major time
component of the runtime of electronic structure codes. As such, it's
worth figuring out how to approximate this function as efficiently as
possible.

I'm going to walk through some approximations in this workbook to try
to get a feeling for the different approximation methods and how
accurate and fast they are.


## Nearly exact Fm function from libint

The Boys function can be computed in a number of different ways,
including using the (lower) incomplete $\gamma$ function.

Here is the version of the Boys function taken from the libint
reference function.

```julia
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
```

## Asymptotic approximation for large `T`

One extremely simple approximation is the asymptotic form that holds
when `T` is big:

```julia
factorial2(n::Int64) = prod(n:-2:1)
Fasymp(m,T) = sqrt(pi/2)*factorial2(2m-1)/(2T)^(m+0.5)
```

Plotting the functions shows that they are very close, even for
relatively small values of `T`.


```julia
x = 0:50
scatter(x,fm_ref.(0,x),label="F0(T)")
scatter!(x,fm_ref.(1,x),label="F1(T)")
scatter!(x,fm_ref.(2,x),label="F2(T)")
plot!(x,Fasymp.(0,x),label="F0(T)")
plot!(x,Fasymp.(1,x),label="F1(T)")
plot!(x,Fasymp.(2,x),label="F2(T)")
```

We can see how close by plotting the error between the exact and the
asymptotic values. Taking `T=10` as a cutoff value looks like it would
be correct to chemical accuracy.

```julia
x1 = 0:30
plot(x1,Fasymp.(0,x1)-fm_ref.(0,x1),label="Err0(T)", yaxis=:log)
plot!(x1,Fasymp.(1,x1)-fm_ref.(1,x1),label="Err1(T)")
plot!(x1,Fasymp.(2,x1)-fm_ref.(2,x1),label="Err2(T)")
```

In the current (5/9/2021) version of the code I have the cutoff set
at `Tcrit=20`, and values lower than this lead to test cases failing.
It's important to note that most other codes have this value set very
high (`Tcrit=117`), but I believe this is also to insure that the
downward recursion can be performed assuming that $\exp(-T)$ is
small.

The asymptotic approximation should be a lot faster. Let's see how much. Here's the time of the nearly exact approximation:

```julia
@btime fm_ref(1.0,10.0)
```

Can't see it in the notebook, but this is coming in at 200-500 nsec.
Let's time the asymptotic value:

```julia
@btime Fasymp(1,10)
```

Comes in at 1-2 nsec. Much faster

## Numerical Integration: Plotting the integrand.

Most of the fast methods of approximating the Boys function involve
clever uses of numerical integration. To get a feeling for what kinds
of resolution such an integration needs, let's look at the integrand
that we will integrate:

```julia
function integrand(m,T)
	fm(u) = u^2m*exp(-T*u*u)
end
```

```julia
plot(integrand(0,0),0,1,label="F0(0)")
plot!(integrand(0,10),0,1,label="F0(10)")
plot!(integrand(1,0),0,1,label="F1(0)")
plot!(integrand(1,10),0,1,label="F1(10)")
plot!(integrand(10,0),0,1,label="F10(0)")
plot!(integrand(10,10),0,1,label="F10(10)")
```

The good news is that these are all smooth and should be
straightforward to integrate numerically. I can see why people thought
of Chebyshev interpolation for these functions.

## Naive numerical approximation using `QuadGK`

Let's just use a canned numerical integration scheme from the `QuadGK`
Julia module to set a baseline for the more sophisticated approaches.

We need to pass in a lower `rtol` so that we can integrate faster.

```julia
@btime quadgk(integrand(1,10),0,1,rtol=1e-4)
```

Using `rtol=1e-4` (which actually reports a error bound of 3.8e-7), we can integrate this in ~240 nsec. It's faster, but not dramatically so.

## Chebyshev interpolation and integration

The method that Gill, Johnson, and Pople [^GJP] use is to divide the
interval into widths of 2Δ, and to use a polynomial interpolation
technique to approximate the integral in this range.

Equation 39 in [^GJP] gives the following formula for Δ to generate an
approximation of ε:

$ \Delta = \left[\frac{2^n(n+1)!\varepsilon}{\max|f^{(n+1)}(X)|}\right]^{1/(n+1)}$

(I think the following is wrong. Upon further reading, I think the
exponent of $f$ represent the derivative, not a power.) Taking
$\max|f|=1$, we can implement this as:

```julia
delta(n,eps=1e-6) = (2^n*factorial(n+1)*eps)^(1/(n+1))
delta(3)
delta(7)
```

## Clenshaw-Curtis quadrature

[^GJP]'s method of integrating a Chebyshev expansion looks to be
equivalent to [Clenshaw-Curtis quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature)[^CCQ].

The Julia module `FastTransforms.jl`[^FTjl] provides nodes and weights
for Clenshaw-Curtis quadrature.

What if we just try to integrate the [0,1] interval in a single
integration step, instead of breaking the integral into Δ-sized
pieces?

First, we need to change the interval from [0,1] to [-1,1], via x=2u-1
⇒ u→(x+1)/2, du→dx/2. This means the integrand goes from

$Fm(T)=\int_0^1 u^{2m}\exp(-Tu^2)du$

to

$Fm(T)=\frac{1}{2}\int_{-1}^1 (x/2+1/2)^{2m}\exp(-T(x/2+1/2)^2)du$

This changes the integrand to:

```julia
function integrandx(m,T)
	fm(x) = (x/2+0.5)^2m*exp(-T*(x/2+0.5)^2)
end
```

Double-check that these still look right when plotted against [-1,1]:

```julia
plot(integrandx(0,0),-1,1,label="F0(0)")
plot!(integrandx(0,10),-1,1,label="F0(10)")
plot!(integrandx(1,0),-1,1,label="F1(0)")
plot!(integrandx(1,10),-1,1,label="F1(10)")
plot!(integrandx(10,0),-1,1,label="F10(0)")
plot!(integrandx(10,10),-1,1,label="F10(10)")
```

## References

[^GJP]: Two-Electron Repulsion Integrals Over Gaussian s Functions,
    IJQC, 40, 745, 1991.
[^CCQ]: [Wikipedia: Clenshaw-Curtis
    quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature).
[^FTjl]: [FastTransforms.jl](https://juliaapproximation.github.io/FastTransforms.jl/v0.2.0/)


