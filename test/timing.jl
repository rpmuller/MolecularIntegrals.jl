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
    # vrr3 turns out to be 2-3 times faster than vrr2.
    @time MolecularIntegrals.vrr2(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza);
    @time MolecularIntegrals.vrr3(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza);
end

vrr_timings()