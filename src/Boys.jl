# Boys.jl contains different implementations and approximations to the Boys function
using OffsetArrays

function boys_array_asymp(mmax,x)
    boys_array = OffsetArray(zeros(mmax+1),0:mmax)
    denom = 0.5/x # Only used for large x, so don't check for small x
    boys_array[1] = sqrt(0.5*pi*denom)
    for m in 2:mmax
        boys_array[m] = boys_array[m-1]*(2m-1)*denom
    end
    return boys_array
end

function boys_array_gamma(mmax,x,SMALL=1e-18)
    x = max(x,SMALL) # Evidently needs underflow protection
    boys_array = OffsetArray(zeros(mmax+1),0:mmax)
    oox = 1/x
    denom = sqrt(oox)
    boys_array[1] = 0.5*denom*gammainc(0.5,x)
    for m in 2:mmax
        denom *= oox
        # Can speed this up more by expressing gamma(m) in terms of gamma(m±1)
        boys_array[m] = fm_ref(m-1,x) #0.5*denom*gammainc(m-0.5,x) 
    end
    return boys_array
end

boys_array_trivial(mmax,x) = [Fgamma(m-1,x) for m in 1:mmax]

"Boys Fgamma function, using the lower incomplete gamma function."
function Fgamma(m,x,SMALL=1e-18,Tcrit=20.0) 
    # Note, most programs use a much larger value for Tcrit (117)
    mhalf = m+0.5
    x = max(x,SMALL) # Evidently needs underflow protection
    if x>Tcrit 
        retval = sqrt(pi/2)*factorial2(2m-1)/(2x)^mhalf
    else
        retval = 0.5*x^-mhalf*gammainc(mhalf,x)
    end
    return retval
end

"gammainc returns the lower incomplete gamma function"
gammainc(a,x) = gamma(a)*gamma_inc(a,x)[1]

# This is the ref function from libint
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


# These were taken from Numerical Recipes

function gammainc_nr(a::Float64,x::Float64)
    # This is the series version of gamma from pyquante. For reasons I don't get, it 
    # doesn't work around a=1. This works alright, but is only a stopgap solution
    # until Julia gets an incomplete gamma function programmed
    if abs(a-1) < 1e-3
        println("Warning: gammainc_series is known to have problems for a ~ 1")
    end
    if x < (a+1.0)
        #Use the series representation
        gam,gln = gser(a,x)
    else 
        #Use continued fractions
        gamc,gln = gcf(a,x)
        gam = 1-gamc
    end
    return exp(gln)*gam
end

function gser(a,x,ITMAX=100,EPS=3e-9)
    # Series representation of Gamma. NumRec sect 6.1.
    gln=loggamma(a)
    if x == 0
        return 0,gln
    end
    ap = a
    delt = s = 1/a
    for i in 1:ITMAX
        ap += 1
        delt *= (x/ap)
        s += delt
        if abs(delt) < abs(s)*EPS
            break
        end
    end
    return s*exp(-x+a*log(x)-gln),gln
end

function gcf(a::Float64,x::Float64,ITMAX::Int64=200,EPS::Float64=3e-9,FPMIN::Float64=1e-30)
    #Continued fraction representation of Gamma. NumRec sect 6.1"
    gln=loggamma(a)
    b=x+1.0-a
    c=1.0/FPMIN
    d=1.0/b
    h=d
    for i in 1:ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if abs(d) < FPMIN
            d=FPMIN
        end
        c=b+an/c
        if abs(c) < FPMIN
            c=FPMIN
        end
        d=1.0/d
        delt=d*c
        h=h*delt
        if abs(delt-1.) < EPS
            break
        end
    end
    gammcf = exp(-x+a*log(x)-gln)*h
    return gammcf,gln
end