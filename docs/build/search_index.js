var documenterSearchIndex = {"docs":
[{"location":"#MolecularIntegrals.jl","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"The goal of MolecularIntegrals.jl is to supply fast and hackable one- and two-electron integrals for computational chemistry calculations.  There are many excellent molecular integral packages available to Julia programmers, but few are written in Julia. This project will explore how fast we can make these integrals while maintaining a readable and hackable code base.","category":"page"},{"location":"#One-Electron-Integrals","page":"MolecularIntegrals.jl","title":"One Electron Integrals","text":"","category":"section"},{"location":"#Two-Electron-Integrals","page":"MolecularIntegrals.jl","title":"Two Electron Integrals","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"MolecularIntegrals.jl supports a slower method based on Huzinaga's work, and a faster set of integrals based on Head-Gordon and Pople's work [HGP].","category":"page"},{"location":"#Huzinaga-Method","page":"MolecularIntegrals.jl","title":"Huzinaga Method","text":"","category":"section"},{"location":"#Head-Gordan/Pople-Recurrance-Relations","page":"MolecularIntegrals.jl","title":"Head-Gordan/Pople Recurrance Relations","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"vrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)","category":"page"},{"location":"#MolecularIntegrals.vrr_array-NTuple{10, Any}","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.vrr_array","text":"vrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)\n\nUse Head-Gordon/Pople's vertical recurrence relations to compute an array of two-electron integrals.\n\nA, B, C, D are the centers of four Gaussian functions. aexpn, bexpn, cexpn, dexpn are their exponents. amax and cmax are related to the sum of the shell angular momenta for the a+b, and c+d shells, respectively.\n\nThe function returns an nxm array, where n is the number of aos in the a+b shell, and m is the number of aos in the c+d shell.\n\nThe speeds of vrr_array and vrr_wide_array are roughly equivalent, and both interfaces are being retained for convenience.\n\n\n\n\n\n","category":"method"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"hrr_array(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D)","category":"page"},{"location":"#Other-Methods","page":"MolecularIntegrals.jl","title":"Other Methods","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"Rys quadrature","category":"page"},{"location":"#References","page":"MolecularIntegrals.jl","title":"References","text":"","category":"section"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"[HGP]: A method for two-electron Gaussian integral and integral derivative   evaluation using recurrence relations. Martin Head-Gordon and John   A. Pople. JCP, 89 (9), 5777, 1988.","category":"page"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"[G]: Molecular Integrals Over Gaussian Basis Functions. Peter M. W. Gill. Adv.   Q. Chem., 25, 141 (1994).","category":"page"},{"location":"","page":"MolecularIntegrals.jl","title":"MolecularIntegrals.jl","text":"[GP]: The Prism Algorithm for Two-Electron Integrals. Peter M. W. Gill and John   A. Pople. IJQC, 40, 753 (1991).","category":"page"}]
}
