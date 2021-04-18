using MolecularIntegrals

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
    ex = 1
    A = B = xyz
    ax,ay,az = bx,by,bz = xyz
    C = D = xyza
    cx,cy,cz = dx,dy,dz = xyza
    aexpn=bexpn=cexpn=dexpn = ex
    @time MolecularIntegrals.vrr2(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza);
    @time MolecularIntegrals.vrr(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza);
end
# Results for above case:
# vrr2 0.390068 seconds (791.39 k allocations: 45.704 MiB, 3.79% gc time, 99.51% compilation time)
# vrr  0.014459 seconds (6.29 k allocations: 289.375 KiB, 97.21% compilation time)

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
    #@time MolecularIntegrals.hrr_r(aexpn,ax,ay,az,aI,aJ,aK, bexpn,bx,by,bz,bI,bJ,bK, cexpn,cx,cy,cz,cI,cJ,cK, dexpn,dx,dy,dz,dI,dJ,dK);
    @time MolecularIntegrals.hrr2(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    @time MolecularIntegrals.hrr3(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    @time MolecularIntegrals.hrr(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
end
# Results for hrr, shells=2,1,2,1, xyz, xyza, ex as above:
# hrr2: 0.677294 seconds (1.03 M allocations: 58.055 MiB, 14.63% gc time, 98.78% compilation time)
# hrr3: 0.356590 seconds (835.99 k allocations: 45.983 MiB, 2.83% gc time, 98.99% compilation time)
# hrr:  0.003234 seconds (44.50 k allocations: 2.051 MiB)
# Results for hrr, shells = 2,2,2,2:
# hrr2 0.647392 seconds (1.26 M allocations: 68.603 MiB, 10.92% gc time, 94.79% compilation time)
# hrr3 0.425849 seconds (981.26 k allocations: 136.233 MiB, 4.71% gc time, 81.78% compilation time)
# hrr  0.017597 seconds (221.12 k allocations: 10.490 MiB)

hrr_timings()