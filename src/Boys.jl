# Boys.jl contains different implementations and approximations to the Boys function

# These were taken from Numerical Recipes

function gammainc(a::Float64,x::Float64)
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