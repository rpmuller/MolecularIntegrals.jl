using MolecularIntegrals, Profile

function profints()
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
    #@profview MolecularIntegrals.hrr2(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    #@profview MolecularIntegrals.hrr3(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    #@profview MolecularIntegrals.hrr1(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
    @profview MolecularIntegrals.hrr(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D);
end

function profall()
    #bfs = build_basis(ethane,"cc-pvdz") #"6-31G")
    #@profview MolecularIntegrals.all_twoe_ints_chrr(bfs)
    bfs = build_basis(ethane,"6-31G") # use a smaller basis set for Rys for now
    @profview MolecularIntegrals.all_twoe_ints_rys(bfs)
end

profall()

