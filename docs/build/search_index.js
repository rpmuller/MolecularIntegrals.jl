var documenterSearchIndex = {"docs":
[{"location":"#MolecularIntegrals.jl","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"The goal of MolecularIntegrals.jl is to supply fast and hackable one- and two-electron integrals for computational chemistry calculations.  There are many excellent molecular integral packages available to Julia programmers, but few are written in Julia. This project will explore how fast we can make these integrals while maintaining a readable and hackable code base.","category":"page"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"MolecularIntegrals.jl strives to leverage the excellent work done by similar projects:","category":"page"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"libints and its Julia bindings Lints.jl\nPyscf, the libcint package.\nJuliaChem.jl's JERI bindings\nPyQuante's python and c/cython integrals, and the experimental Julia version","category":"page"},{"location":"#Basis-functions-and-other-structures","page":"MolecularIntegrals.jl","title":"Basis functions and other structures","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"PGBF\npgbf\nCGBF\ncgbf\nShell\nBasis","category":"page"},{"location":"#MolecularIntegrals.PGBF","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.PGBF","text":"PGBF(expn,xyz,I,J,K,norm)\n\nCreate a primitive Gaussian basis function      g(x,y,z) = norm * (x-x0)^I (y-y0)^J (z-z0)^K exp(-expn*r^2) The function parameters xyz correspond to [x0,y0,z0].\n\n\n\n\n\n","category":"type"},{"location":"#MolecularIntegrals.pgbf","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.pgbf","text":"pgbf(expn,x=0,y=0,z=0,I=0,J=0,K=0,norm=1)\n\nHelper function to create a normalized PGBF with some optional defaults set.    \n\n\n\n\n\n","category":"function"},{"location":"#MolecularIntegrals.CGBF","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.CGBF","text":"CGBF(xyz,I,J,K,norm,[pbgfs],[coefs])\n\nCreate a contracted Gaussian basis function made of  the functions in [pgbfs] with coefficients [coefs]. Also track the origin xyz and powers I,J,K.\n\n\n\n\n\n","category":"type"},{"location":"#MolecularIntegrals.cgbf","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.cgbf","text":"cgbf(expn,x=0,y=0,z=0,I=0,J=0,K=0,norm=1)\n\nHelper function to create a CGBF with optional defaults.\n\n\n\n\n\n","category":"function"},{"location":"#MolecularIntegrals.Shell","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.Shell","text":"Shell(xyz,L,expns,coeffs)\n\nStructure for a basis function shell, containing multiple CGBFs of different angular momenta.        \n\n\n\n\n\n","category":"type"},{"location":"#MolecularIntegrals.Basis","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.Basis","text":"Basis(cgbfs,shells,ishell,mshell)\n\nStructure to hold a basis set, with info about shells and other data\n\n\n\n\n\n","category":"type"},{"location":"#One-Electron-Integrals","page":"MolecularIntegrals.jl","title":"One Electron Integrals","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"overlap\nkinetic\nnuclear_attraction","category":"page"},{"location":"#MolecularIntegrals.overlap","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.overlap","text":"overlap(a::PGBF,b::PGBF)\n\nCompute the overlap between primitive Gaussian basis functions a and b.\n\n\n\n\n\noverlap(a::CGBF,b::CGBF)\n\nCompute the overlap between contracted Gaussian basis functions a and b.\n\n\n\n\n\noverlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)\n\nCompute the overlap between primitive Gaussian basis functions  defined by centers ax,ay,az, bx,by,bz, powers aI,aJ,aK bI,bJ,bK, and exponents aexpn,bexpn.\n\n\n\n\n\noverlap(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK)\n\nCompute the overlap between primitive Gaussian basis functions  defined by centers axyz, bxyz, powers aI,aJ,aK bI,bJ,bK, and exponents aexpn,bexpn.\n\n\n\n\n\n","category":"function"},{"location":"#MolecularIntegrals.kinetic","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.kinetic","text":"kinetic(a::PGBF,b::PGBF)\n\nCompute the kinetic energy between primitive Gaussian basis functions a and b.\n\n\n\n\n\nkinetic(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)\n\nCompute the kinetic energy between primitive Gaussian basis functions  defined by centers ax,ay,az, bx,by,bz, powers aI,aJ,aK bI,bJ,bK, and exponents aexpn,bexpn.\n\n\n\n\n\nkinetic(aexpn,axyz,aI,aJ,aK,bexpn,bxyz,bI,bJ,bK)\n\nCompute the kinetic energy between primitive Gaussian basis functions  defined by centers ax,yz, bxyz, powers aI,aJ,aK bI,bJ,bK, and exponents aexpn,bexpn.\n\n\n\n\n\nkinetic(a::CGBF,b::CGBF)\n\nCompute the kinetic energy between contracted Gaussian basis functions a and b.\n\n\n\n\n\n","category":"function"},{"location":"#MolecularIntegrals.nuclear_attraction","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.nuclear_attraction","text":"nuclear_attraction(a::PGBF,b::PGBF,cxyz)\n\nCompute the nuclear attraction energy between primitive Gaussian basis functions a and b and center cxyz.\n\n\n\n\n\nnuclear_attraction(a::PGBF,b::PGBF,c::Atom)\n\nCompute the nuclear attraction energy between primitive Gaussian basis functions a and b and atom c.\n\n\n\n\n\nnuclear_attraction(a::PGBF,b::PGBF,m::Vector{Atom})\n\nCompute the sum of nuclear attraction energy between primitive Gaussian basis  functions a and b and the vector of atoms in m.\n\n\n\n\n\nnuclear_attraction(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,cx,cy,cz)\n\nCompute the nuclear attraction energy between primitive Gaussian basis functions  defined by ax,ay,az, bx,by,bz, powers aI,aJ,aK, bI,bJ,bK, and exponents aexpn,bexpn and center cxyz.\n\n\n\n\n\nnuclear_attraction(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,cx,cy,cz)\n\nCompute the nuclear attraction energy between primitive Gaussian basis functions  defined by ax,ay,az, bx,by,bz, powers aI,aJ,aK, bI,bJ,bK, and exponents aexpn,bexpn and center cxyz.\n\n\n\n\n\nnuclear_attraction(a::CGBF,b::CGBF,cxyz)\n\nCompute the nuclear attraction energy between contracted Gaussian  basis functions a and b and center cxyz.\n\n\n\n\n\nnuclear_attraction(a::CGBF,b::CGBF,c::Atom)\n\nCompute the nuclear attraction energy between contracted Gaussian  basis functions a and b and atom c.\n\n\n\n\n\nnuclear_attraction(a::CGBF,b::CGBF,m::Vector{Atom})\n\nCompute the sum of nuclear attraction energy between contracted Gaussian  basis functions a and b and vector of atoms in m.\n\n\n\n\n\n","category":"function"},{"location":"#Two-Electron-Integrals","page":"MolecularIntegrals.jl","title":"Two Electron Integrals","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"MolecularIntegrals.jl supports a slower method based on Huzinaga's work, and a faster set of integrals based on Head-Gordon and Pople's work [HGP].","category":"page"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"The interface to the Huzinaga integrals is","category":"page"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"coulomb","category":"page"},{"location":"#MolecularIntegrals.coulomb","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.coulomb","text":"coulomb(aexpn,ax,ay,az,aI,aJ,aK,\n    bexpn,bx,by,bz,bI,bJ,bK,\n    cexpn,cx,cy,cz,cI,cJ,cK,\n    dexpn,dx,dy,dz,dI,dJ,dK)\n\nCompute the coulomb repulsion between four primitive Gaussian basis  functions defined by their exponents, xyz coordinates, and IJK powers.        \n\n\n\n\n\ncoulomb(a::PGBF,b::PGBF,c::PGBF,d::PGBF)\n\nCompute the coulomb repulsion between four primitive Gaussian basis  functions a,b,c,d.        \n\n\n\n\n\ncoulomb(a::CGBF,b::CGBF,c::CGBF,d::CGBF)\n\nCompute the coulomb repulsion between four contracted Gaussian basis  functions a,b,c,d.        \n\n\n\n\n\n","category":"function"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"These are useful for checking the values of other integrals, and have a  simple, mature interface.","category":"page"},{"location":"#Head-Gordon/Pople-Recurrence-Relations","page":"MolecularIntegrals.jl","title":"Head-Gordon/Pople Recurrence Relations","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"Integrals are computed using Head-Gordon and Pople's[HGP] recurrence relations using vertical (VRR) and horizontal (HRR) recurrence relations. In the notation of [HGP], vertical recurrence relations  construct integrals of the form [a0|c0] from kernels [00|00]^m via eq 6. Horizontal recurrence relations  construct integrals of the form [ab|cd] from kernels [a0|c0].","category":"page"},{"location":"#Vertical-recurrence-relations","page":"MolecularIntegrals.jl","title":"Vertical recurrence relations","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"MolecularIntegrals.jl returns VRRs either as a 2d array using vrr_array, or as  a 6d array using vrr_widearray:","category":"page"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"vrr_array\nvrr_widearray","category":"page"},{"location":"#MolecularIntegrals.vrr_array","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.vrr_array","text":"vrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)\n\nUse Head-Gordon/Pople's vertical recurrence relations to compute an array of two-electron integrals.\n\nA, B, C, D are the centers of four Gaussian functions. aexpn, bexpn, cexpn, dexpn are their exponents. amax and cmax are related to the sum of the shell angular momenta for the a+b, and c+d shells, respectively.\n\nThe function returns an nxm array, where n is the number of aos in the a+b shell, and m is the number of aos in the c+d shell.\n\nThe speeds of vrr_array and vrr_widearray are roughly equivalent, and both interfaces are being retained for convenience.\n\n\n\n\n\n","category":"function"},{"location":"#MolecularIntegrals.vrr_widearray","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.vrr_widearray","text":"vrr_widearray(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)\n\nUse Head-Gordon/Pople's vertical recurrence relations to compute an array of integrals.\n\nA, B, C, D are the centers of four Gaussian functions. aexpn, bexpn, cexpn, dexpn are their exponents. amax and cmax are related to the sum of the shell angular momenta for the a+b, and c+d shells, respectively.\n\nThe function returns a six-dimensional array over the possible powers of the a+b and c+d shell functions.\n\nThe speeds of vrr_array and vrr_widearray are roughly equivalent, and both interfaces are being retained for convenience.\n\n\n\n\n\n","category":"function"},{"location":"#Horizontal-recurrence-relations","page":"MolecularIntegrals.jl","title":"Horizontal recurrence relations","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"hrr_array\nhrr_dict","category":"page"},{"location":"#MolecularIntegrals.hrr_array","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.hrr_array","text":"hrr_array(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)\n\nUse Head-Gordon/Pople's horizontal recurrence relations to compute an array of two-electron integrals.\n\nA, B, C, D are the centers of four Gaussian functions. aexpn, bexpn, cexpn, dexpn are their exponents. ashell, bshell, cshell and dshell are the shell angular momenta for the a, b, c, and d shells, respectively.\n\nThe function returns a (k,l,m,n)-dimensional array, where the  dimensions correspond to the number of aos in the a,b,c,d shells.        \n\n\n\n\n\n","category":"function"},{"location":"#MolecularIntegrals.hrr_dict","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.hrr_dict","text":"hrr_dict(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)\n\nUse Head-Gordon/Pople's horizontal recurrence relations to compute an array of two-electron integrals.\n\nA, B, C, D are the centers of four Gaussian functions. aexpn, bexpn, cexpn, dexpn are their exponents. ashell, bshell, cshell and dshell are the shell angular momenta for the a, b, c, and d shells, respectively.\n\nThe function returns a dictionary containing entries for the relevant integrals. E.g., hrrs[ax,ay,az,bx,by,bz,cx,cy,cz,dpx,dpy,dpz] contains the integral corresponding to the bfs with powers ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz.\n\nhrr_dict is slower than hrr_array, but is kept for convenience.        \n\n\n\n\n\n","category":"function"},{"location":"#References","page":"MolecularIntegrals.jl","title":"References","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"[HGP]: A method for two-electron Gaussian integral and integral derivative   evaluation using recurrence relations.    Martin Head-Gordon and John A. Pople. JCP, 89 (9), 5777, 1988.","category":"page"}]
}