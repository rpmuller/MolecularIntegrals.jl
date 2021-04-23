using MolecularIntegrals
using BenchmarkTools

function time_ethane()
    bfs = build_basis(ethane,"6-31G")
    @time MolecularIntegrals.all_twoe_ints(bfs);
    #@time MolecularIntegrals.all_twoe_ints(bfs,MolecularIntegrals.coulomb_hgp);
end

function vrr_timings()
    println("VRR timings")
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

    print("# vrr2($amax,$cmax) ")
    @btime MolecularIntegrals.vrr2($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    print("# vrr1($amax,$cmax) ")
    @btime MolecularIntegrals.vrr1($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    print("# vrr5($amax,$cmax) ")
    @btime MolecularIntegrals.vrr($amax,$cmax, $ex,$ex,$ex,$ex, $xyz,$xyz,$xyza,$xyza);
    nothing
end

# vrr2(4,4)   17.560 ms (166510 allocations: 6.46 MiB)
# vrr1(4,4)   3.572 ms (104351 allocations: 3.85 MiB)
# vrr5(4,4)   3.864 ms (81789 allocations: 1.70 MiB)

function hrr_timings()
    println("HRR timing")
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
    print("# hrr2($ashell,$bshell,$cshell,$dshell) ")
    @btime MolecularIntegrals.hrr2($ashell,$bshell,$cshell,$dshell, $aexpn,$bexpn,$cexpn,$dexpn, $A,$B,$C,$D);
    print("# hrr1($ashell,$bshell,$cshell,$dshell) ")
    @btime MolecularIntegrals.hrr1($ashell,$bshell,$cshell,$dshell, $aexpn,$bexpn,$cexpn,$dexpn, $A,$B,$C,$D);
    print("# hrr5($ashell,$bshell,$cshell,$dshell) ")
    @btime MolecularIntegrals.hrr($ashell,$bshell,$cshell,$dshell, $aexpn,$bexpn,$cexpn,$dexpn, $A,$B,$C,$D);
    nothing
end

# hrr2(2,2,2,2)   28.719 ms (282814 allocations: 12.60 MiB)
# hrr1(2,2,2,2)   20.398 ms (289908 allocations: 8.86 MiB)
# hrr5(2,2,2,2)   19.829 ms (257436 allocations: 6.44 MiB)

vrr_timings()
hrr_timings()