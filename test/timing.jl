using MolecularIntegrals
using BenchmarkTools

function time_ethane()
    bfs = build_basis(ethane,"6-31G")
    @time MolecularIntegrals.all_twoe_ints(bfs);
    #@time MolecularIntegrals.all_twoe_ints(bfs,MolecularIntegrals.coulomb_hgp);
end

function vrr_timings()
    x=y=z=0.0
    xyz = [x,y,z]
    xyza = xyz + [0.1,0.05,0.025]
    xa,ya,za = xyza
    ex = 1.0
    A = B = xyz
    ax,ay,az = bx,by,bz = xyz
    C = D = xyza
    cx,cy,cz = dx,dy,dz = xyza
    amax = cmax = 4

    #print("# vrr_widearray($amax,$cmax) ")
    #@btime MolecularIntegrals.vrr_widearray($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    print("# vrr($amax,$cmax)      ")
    @btime MolecularIntegrals.vrr($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    #print("# vrr_aoloop($amax,$cmax) ")
    #@btime MolecularIntegrals.vrr_aoloop($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    nothing
end

# vrr2(4,4)            17.560 ms (166510 allocations: 6.46 MiB)
# vrr_widearray(4,4)   3.572 ms (104351 allocations: 3.85 MiB)
# vrr(4,4)        3.864 ms (81789 allocations: 1.70 MiB)

function hand_vrr_timings()
    x=y=z=0.0
    xyz = [x,y,z]
    xyza = xyz + [0.1,0.05,0.025]
    xa,ya,za = xyza
    ex = 1.0
    A = B = xyz
    ax,ay,az = bx,by,bz = xyz
    C = D = xyza
    cx,cy,cz = dx,dy,dz = xyza
    amax = cmax = 1
    #print("# vrr_widearray($amax,$cmax) ")
    #@btime MolecularIntegrals.vrr_widearray($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    print("# vrr($amax,$cmax)      ")
    @btime MolecularIntegrals.vrr($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    print("# vrr_pp               ")
    @btime MolecularIntegrals.vrr_pp($ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
end

# vrr_widearray(1,1)   23.402 μs (519 allocations: 20.45 KiB)
# vrr5(1,1)             33.262 μs (384 allocations: 15.94 KiB)
# vrr_pp                 1.586 μs (16 allocations: 2.66 KiB)


function hrr_timings()
    x=y=z=0.0
    xyz = [x,y,z]
    xyza = xyz + [0.1,0.05,0.025]
    xa,ya,za = xyza
    ex = 1
    A = B = xyz
    ax,ay,az = bx,by,bz = xyz
    C = D = xyza
    cx,cy,cz = dx,dy,dz = xyza
    aexpn=bexpn=cexpn=dexpn = ex
    ashell,bshell,cshell,dshell = (2,2,2,2)
    print("# hrr($ashell,$bshell,$cshell,$dshell) ")
    @btime MolecularIntegrals.hrr($ashell,$bshell,$cshell,$dshell, $aexpn,$bexpn,$cexpn,$dexpn, $A,$B,$C,$D);
    #print("# hrr_dict($ashell,$bshell,$cshell,$dshell)  ")
    #@btime MolecularIntegrals.hrr_dict($ashell,$bshell,$cshell,$dshell, $aexpn,$bexpn,$cexpn,$dexpn, $A,$B,$C,$D);
    nothing
end

# hrr2(2,2,2,2)   28.719 ms (282814 allocations: 12.60 MiB)
# hrr1(2,2,2,2)   20.398 ms (289908 allocations: 8.86 MiB)
# hrr5(2,2,2,2)   19.829 ms (257436 allocations: 6.44 MiB)

function ethane_timing()
    for bname in ["sto3g","6-31G"]
        bfs = build_basis(ethane,bname)
        println("Ethane: nbf=$(length(bfs))")
        fetcher = MolecularIntegrals.eri_fetcher(bfs)
        print("Huzinaga ")
        @btime MolecularIntegrals.all_twoe_ints($bfs);
        print("HGP      ")
        @btime MolecularIntegrals.all_twoe_ints_chrr($bfs);
    end
    return nothing
end

#vrr_timings()
#hand_vrr_timings()
#hrr_timings()
ethane_timing()