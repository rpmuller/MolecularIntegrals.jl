using MolecularIntegrals

function time_ethane()
    bfs = build_basis(ethane,"6-31G")
    @time MolecularIntegrals.all_twoe_ints(bfs);
    #@time MolecularIntegrals.all_twoe_ints(bfs,MolecularIntegrals.coulomb_hgp);
end

function vrr_timings()
    A = B = xyz
    ax,ay,az = bx,by,bz = xyz
    C = D = xyza
    cx,cy,cz = dx,dy,dz = xyza
    aexpn=bexpn=cexpn=dexpn = ex
    val2 = MolecularIntegrals.vrr2(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza)
    val3 = MolecularIntegrals.vrr3(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza)
end