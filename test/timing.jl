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
    @time MolecularIntegrals.vrr2(4,4, ex,ex,ex,ex, xyz,xyz,xyza,xyza);
    @time MolecularIntegrals.vrr1(4,4, ex,ex,ex,ex, xyz,xyz,xyza,xyza);
    @time MolecularIntegrals.vrr(4,4, ex,ex,ex,ex, xyz,xyz,xyza,xyza);
    nothing
end
# Results for 2,2 case:
# vrr2 0.001125 seconds (9.06 k allocations: 422.250 KiB)
# vrr1 0.000544 seconds (5.35 k allocations: 238.000 KiB)
# vrr  0.000525 seconds (3.94 k allocations: 152.641 KiB)

# Results for 4,4 case:
# vrr2 0.039209 seconds (169.34 k allocations: 7.163 MiB)
# vrr1 0.012553 seconds (107.18 k allocations: 4.554 MiB)
# vrr  0.009616 seconds (84.62 k allocations: 2.401 MiB)


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
    @time MolecularIntegrals.hrr2(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    @time MolecularIntegrals.hrr3(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    @time MolecularIntegrals.hrr1(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    @time MolecularIntegrals.hrr(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    nothing
end
# Results for hrr, shells = 2,2,2,2:
# hrr2 0.029898 seconds (290.39 k allocations: 13.413 MiB)
# hrr3 0.038984 seconds (182.71 k allocations: 93.745 MiB, 9.30% gc time)
# hrr1 0.075930 seconds (448.54 k allocations: 28.874 MiB)
# hrr  0.028568 seconds (265.20 k allocations: 7.264 MiB)

hrr_timings()