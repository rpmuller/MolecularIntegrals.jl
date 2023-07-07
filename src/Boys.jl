# Boys.jl contains different implementations and approximations to the Boys function
using SpecialFunctions

export farray, Fm, Fm_asymp, Fm_ref, farray_full, farray_recur

# Functions to make arrays out of multiple m values:
farray_full(mmax,T,fm=Fm) = fm.(0:mmax-1,T)

@inline function farray_recur(mmax,T,fm=Fm)
    Fms = zeros(Float64,mmax)
    Fms[mmax] = fm(mmax-1,T)
    emt = exp(-T)
    for m in mmax-1:-1:1
	    Fms[m] = (2*T*Fms[m+1]+emt)/(2*(m-1)+1)
    end
    return Fms
end

function farray(mmax,T,fm=Fm,Tmax=30)
    if T < Tmax return farray_recur(mmax,T,fm) end
    Fms = zeros(Float64,mmax)
    Fms[mmax] = Fm_asymp(mmax-1,T)
    for m in mmax-1:-1:1
	    Fms[m] = T*Fms[m+1]/(m-0.5)
    end
    return Fms		
end

function Fm(m,T,SMALL=1e-18)
    mhalf = m+0.5
    #T = max(T,SMALL) # Underflow; is there a better way to fix this?
    if T<SMALL return 1/(2m+1) end # This seems to work??
    return gamma(mhalf)*gamma_inc(mhalf,T,0)[1]/(2T^(mhalf))
end

Fm_asymp(m,T) = sqrt(pi/2)*factorial2(2m-1)/(2T)^(m+0.5)

# This is the ref function from libint
function Fm_ref(m,T,eps = 1e-10)
    denom = (m + 0.5)
    term = exp(-T) /2denom
    old_term = 0.0
    sum = term
    while term > sum*eps || old_term < term
        denom += 1
        old_term = term
        term = old_term * T / denom
        sum += term
    end
    return sum
end

# Old functions: hopefully no longer used

# This was an earlier attempt to evaluate the Fm values on a grid and
# interpolate between them
Fgamma(m,T) = Fm(m,T) # backwards compatibility in case something still calls this
# "Boys Fgamma function, using the lower incomplete gamma function."
# @inline function Fgamma(m,T,SMALL=1e-18,Tcrit=20.0) 
#     # Note, most programs use a much larger value for Tcrit (117)
#     mhalf = m+0.5
#     T = max(T,SMALL) # Evidently needs underflow protection
#     if T>Tcrit 
#         retval = sqrt(pi/2)*factorial2(2m-1)/(2T)^mhalf
#     else
#         #retval = 0.5*T^-mhalf*gammainc_nr(mhalf,T)
#         retval = fmspline(m,T)
#     end
#     return retval
# end
# using Interpolations
# function make_interpolator(mmax=10,Tmax=20.0)
#     Tgrid = 0:0.005:Tmax
#     mgrid = 0:mmax
#     fmvalues = [Fm(m,T) for m in mgrid, T in Tgrid]
#     itp = interpolate(fmvalues, BSpline(Cubic(Line(OnGrid()))))
#     sitp = scale(itp,mgrid,Tgrid)
#     return sitp
# end
# const fmspline = make_interpolator()
