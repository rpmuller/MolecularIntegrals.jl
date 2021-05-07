using MolecularIntegrals

open("src/HGPgen3.jl","w") do io
    for amax in 0:3
        for cmax in 0:(amax-1)
            f = MolecularIntegrals.vrr_autogen(amax,cmax)
            write(io,f)
            f = MolecularIntegrals.vrr_autogen(cmax,amax)
            write(io,f)
        end
        f = MolecularIntegrals.vrr_autogen(amax,amax)
        write(io,f)
    end
end

