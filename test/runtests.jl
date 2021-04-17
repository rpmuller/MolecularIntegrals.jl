using MolecularIntegrals, Test
using OffsetArrays

# Define functions to be used throughout
s = pgbf(1.0)
px = pgbf(1.0,0,0,0,1,0,0)

c = cgbf(0.0,0.0,0.0)
addbf!(c,1,1)

c2 = cgbf(0,0,0)
addbf!(c2,1,0.2)
addbf!(c2,0.5,0.2)

@testset verbose = true "MolecularIntegrals.jl" begin

    @testset "Utils" begin
        @test factorial2(6)==48
        @test collect(ipairs(3)) == [(1,1),(2,1),(2,2),(3,1),(3,2),(3,3)]
        @test collect(spairs(3))== [(2,1),(3,1),(3,2)]
        @test collect(rpairs(2)) == [(1,1),(1,2),(2,1),(2,2)]
        @test iindex(1,1,1,1) == 1
        @test iindex(1,1,1,2) == iindex(1,1,2,1) == iindex(1,2,1,1) == iindex(2,1,1,1) == 2
        @test iindex(1,1,2,2) == iindex(2,2,1,1) == 4
    end

    @testset "Basis" begin
        @test s(0,0,0) ≈ 0.71270547
        @test px(0,0,0) ≈ 0
        @test c(0,0,0) ≈ 0.71270547
        @test overlap(c2,c2) ≈ 1
    end

    @testset "OneInts" begin
        @test kinetic(1.0, 0.0,0.0,0.0, 0,0,0, 1.0, 0.0,0.0,0.0, 0,0,0) ≈ 2.9530518648229536
        @test kinetic(s,s) ≈ 1.5
        @test kinetic(c,c) ≈ 1.5

        @test nuclear_attraction(s,s,[0.0,0.0,0.0]) ≈ -1.59576912
        @test nuclear_attraction(c,c,[0.0,0.0,0.0]) ≈ -1.59576912
        @test nuclear_attraction(px,px,[0.0,0.0,0.0]) ≈ -1.06384608

        @test overlap(s,s) ≈ 1
        @test overlap(px,px) ≈ 1
        @test overlap(s,px) ≈ 0
        @test overlap(c,c) ≈ 1
    end

    @testset "Low level OneInts" begin
        @test MolecularIntegrals.overlap1d(0,0,0.0,0.0,1.0) == 1
        @test MolecularIntegrals.gaussian_product_center(s,s) == [0,0,0]
        @test MolecularIntegrals.binomial_prefactor(0,0,0,0.0,0.0) == 1
        
        @test MolecularIntegrals.Aterm(0,0,0,0,0,0,0,0,0) == 1.0
        @test MolecularIntegrals.MolecularIntegrals.Aarray(0,0,0,0,0,1) == OffsetArray([1.0],0:0)
        @test MolecularIntegrals.Aarray(0,1,1,1,1,1) == OffsetArray([1.0, -1.0],0:1)
        @test MolecularIntegrals.Aarray(1,1,1,1,1,1) == OffsetArray([1.5, -2.5, 1.0],0:2)
        @test MolecularIntegrals.Aterm(0,0,0,0,0,0,0,0,1) == 1.0
        @test MolecularIntegrals.Aterm(0,0,0,0,1,1,1,1,1) == 1.0
        @test MolecularIntegrals.Aterm(1,0,0,0,1,1,1,1,1) == -1.0
        @test MolecularIntegrals.Aterm(0,0,0,1,1,1,1,1,1) == 1.0
        @test MolecularIntegrals.Aterm(1,0,0,1,1,1,1,1,1) == -2.0
        @test MolecularIntegrals.Aterm(2,0,0,1,1,1,1,1,1) == 1.0
        @test MolecularIntegrals.Aterm(2,0,1,1,1,1,1,1,1) == -0.5
        @test MolecularIntegrals.Aterm(2,1,0,1,1,1,1,1,1) == 0.5
    end 

    @testset "ERI tests" begin
        @test coulomb(1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0) ≈ 4.37335458
        @test coulomb(s,s,s,s) ≈ 1.128379167
        @test coulomb(c,c,c,c) ≈ 1.128379167
        @test coulomb(1, 0,0,0, 0,0,0, 1, 0,0,1, 0,0,0, 1, 0,0,0, 0,0,0, 1, 0,0,1, 0,0,0) ≈ 1.6088672396
        @test coulomb(c,c,c2,c2) ≈ 1.0343247
        @test coulomb(s,s,px,px) ≈ 1.16599181
    end

    @testset "VRR tests" begin
        ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
        aexpn=bexpn=cexpn=dexpn=1.0
        aI=aJ=aK=0
        cI=cJ=cK=0
        M=0

        for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in [
            (0.,0.,0., 0,0,0, 0,0,0, 4.37335456733),
            (0.,0.,0., 1,0,0, 1,0,0, 0.182223107579),
            (0.,0.,0., 0,1,0, 0,1,0, 0.182223107579),
            (0.,0.,0., 0,0,1, 0,0,1, 0.182223107579),

            (0.,0.,0., 2,0,0, 2,0,0,  0.223223306785),
            (0.,0.,0., 0,2,0, 0,2,0,  0.223223306785),
            (0.,0.,0., 0,0,2, 0,0,2,  0.223223306785),

            (1.,2.,3., 1,0,0, 1,0,0, -5.63387712455e-06),
            (1.,2.,3., 0,1,0, 0,1,0, -0.000116463120359),
            (1.,2.,3., 0,0,1, 0,0,1, -0.000301178525749),

            (1.,2.,3., 2,0,0, 2,0,0, 0.00022503308545040895),
            (1.,2.,3., 0,2,0, 0,2,0, 0.0006102470883881907),
            (1.,2.,3., 0,0,2, 0,0,2, 0.0013427831014563411),

            (0.,0.,0., 1,1,0, 1,1,0, 0.0136667330685),
            (0.,0.,0., 0,1,1, 0,1,1, 0.0136667330685),
            (0.,0.,0., 1,0,1, 1,0,1, 0.0136667330685),

            (3.,2.,1., 1,1,0, 1,1,0, 5.976771621486971e-5),
            (3.,2.,1., 0,1,1, 0,1,1, 1.5742904443905067e-6),
            (3.,2.,1., 1,0,1, 1,0,1, 4.00292848649699e-6)
        ]

            val1 = MolecularIntegrals.vrr(
                aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
            val2 = MolecularIntegrals.vrr(
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
                aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
            @test isapprox(val1,val2)
            @test isapprox(val1,result)
        end
        
    end
    @testset "HRR tests" begin
        ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
        aexpn=bexpn=cexpn=dexpn=1.0
        aI=aJ=aK=0
        bI,bJ,bK = 1,0,1
        cI=cJ=cK=0
        dI,dJ,dK = 1,0,1

        for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in [
            (0.,0.,0., 0,0,0, 0,0,0, 0.0136667330685),
            (0.,0.,0., 1,0,0, 1,0,0, 0.00821630976139),
            (0.,0.,0., 0,1,0, 0,1,0, 0.00122024402397),
            (0.,0.,0., 0,0,1, 0,0,1, 0.00821630976139),

            (0.,0.,0., 2,0,0, 2,0,0,   0.0039759617781),
            (0.,0.,0., 0,2,0, 0,2,0,   0.000599953311785),
            (0.,0.,0., 0,0,2, 0,0,2,  0.0039759617781),

            (1.,2.,3., 1,0,0, 1,0,0, -1.1851316496333975e-6),
            (1.,2.,3., 0,1,0, 0,1,0,  -4.669991667384835e-6),
            (1.,2.,3., 0,0,1, 0,0,1, -3.474373852654044e-5),

            (1.,2.,3., 2,0,0, 2,0,0, 2.81002247462e-6),
            (1.,2.,3., 0,2,0, 0,2,0, 7.09856891538e-6),
            (1.,2.,3., 0,0,2, 0,0,2, 3.62153023224e-5),

            (0.,0.,0., 1,1,0, 1,1,0, 0.000599953311785),
            (0.,0.,0., 0,1,1, 0,1,1, 0.000599953311785),
            (0.,0.,0., 1,0,1, 1,0,1, 0.0116431617287),

            (3.,2.,1., 1,1,0, 1,1,0, 7.37307761485e-6),
            (3.,2.,1., 0,1,1, 0,1,1, 2.5333243119843164e-7),
            (3.,2.,1., 1,0,1, 1,0,1, 2.452115184675799e-6)
        ]
            val1 = MolecularIntegrals.hrr(
                aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
            val2 = MolecularIntegrals.hrr(
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK,
                aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
            @test val1 ≈ val2 ≈ result
            
        end
    end

    @testset "Molecular basis" begin
        bfs = build_basis(h2)
        @test length(bfs)==2
        l,r = bfs
        @test overlap(l,l) ≈ 1
        @test overlap(r,r) ≈ 1
        @test overlap(l,r) ≈ 0.6647387449282997
        @test kinetic(l,l) ≈ 0.76003188
        @test kinetic(r,r) ≈ 0.76003188
        @test kinetic(l,r) ≈ 0.24141861181119084
        @test coulomb(l,l,l,l) ≈ 0.7746059439196398
        @test coulomb(r,r,r,r) ≈ 0.7746059439196398
        @test coulomb(l,l,r,r) ≈ 0.5727937653511646
        @test coulomb(l,l,l,r) ≈ 0.4488373301593464
        @test coulomb(l,r,l,r) ≈ 0.3025451156654606

        @test MolecularIntegrals.lvalue["S"] == 0
    end

    @testset "HGP2 tests" begin
        x=y=z=0.0
        xyz = [x,y,z]
        xyza = xyz + [0.1,0.05,0.025]
        xa,ya,za = xyza
        ex = 1
        s0 = pgbf(ex,x,y,z)
        sa = pgbf(ex,xa,ya,za)
        px = pgbf(ex,x,y,z,1,0,0)
        

        # prunem function removes the m values, since at the end of vrr we only need the ones where m=0
        testd = Dict((1,2,3) => 4, (1,2,0)=> 3)
        @test MolecularIntegrals.prunem(testd) == Dict((1,2) => 3)

        @test MolecularIntegrals.unit(3,1) == [1,0,0]

        # Test the vrr2 code:
        val = MolecularIntegrals.vrr2(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza) 

        for (I1,J1,K1,I2,J2,K2) in [(1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0),
                                    (0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 1),
                                    (1, 0, 0, 1, 0, 0), (0, 1, 0, 0, 1, 0), (0, 0, 1, 0, 0, 1),
                                    (2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0), (0, 0, 2, 0, 0, 0),
                                    (0, 0, 0, 2, 0, 0), (0, 0, 0, 0, 2, 0), (0, 0, 0, 0, 0, 2),
                                    (2, 0, 0, 2, 0, 0), (0, 2, 0, 0, 2, 0), (0, 0, 2, 0, 0, 2)]
            @test val[(I1,J1,K1,I2,J2,K2)] == MolecularIntegrals.vrr(ex, x,y,z, I1,J1,K1, ex, x,y,z, ex, xa,ya,za, I2,J2,K2, ex, xa,ya,za,0)
        end

        # Test the vrr3 code:
        val3 = MolecularIntegrals.vrr2(2,2, ex,ex,ex,ex, xyz,xyz,xyza,xyza)
        for (I1,J1,K1,I2,J2,K2) in [(1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0),
                                    (0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 1),
                                    (1, 0, 0, 1, 0, 0), (0, 1, 0, 0, 1, 0), (0, 0, 1, 0, 0, 1),
                                    (2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0), (0, 0, 2, 0, 0, 0),
                                    (0, 0, 0, 2, 0, 0), (0, 0, 0, 0, 2, 0), (0, 0, 0, 0, 0, 2),
                                    (2, 0, 0, 2, 0, 0), (0, 2, 0, 0, 2, 0), (0, 0, 2, 0, 0, 2)]
            @test val[(I1,J1,K1,I2,J2,K2)] == val3[I1,J1,K1,I2,J2,K2]
        end

        # Using dense arrays for this is still pretty small: length(vrr3(d,d)) = 100
        #@show length(val3)

        # Test against coulomb() as well: The second test doesn't work:
        @test val[(0, 0, 0, 0, 0, 0)] ≈ MolecularIntegrals.coulomb(s0,s0,sa,sa)/s0.norm^2/sa.norm^2
        #@test val[(1, 0, 0, 0, 0, 0)] ≈ MolecularIntegrals.coulomb(px,s0,sa,sa)/px.norm/s0.norm/sa.norm^2

        # HRR tests:
        A = B = xyz
        ax,ay,az = bx,by,bz = xyz
        C = D = xyza
        cx,cy,cz = dx,dy,dz = xyza
        aexpn=bexpn=cexpn=dexpn = ex

        for (ashell,bshell,cshell,dshell) in [(0,0,0,0),(1,0,0,0),(1,0,0,1),(1,1,0,0),(1,1,0,1),
                (2,1,0,1),(2,1,2,0),(2,1,2,1),(3,1,0,0),(3,1,3,1),(0,2,0,0),(2,1,0,2)] 
            aI,aJ,aK = MolecularIntegrals.shell_indices[ashell][1]
            bI,bJ,bK = MolecularIntegrals.shell_indices[bshell][1]
            cI,cJ,cK = MolecularIntegrals.shell_indices[cshell][1]
            dI,dJ,dK = MolecularIntegrals.shell_indices[dshell][1]

            hrr2_vals = MolecularIntegrals.hrr2(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
            hrr_val = MolecularIntegrals.hrr(aexpn,ax,ay,az,aI,aJ,aK, bexpn,bx,by,bz,bI,bJ,bK, cexpn,cx,cy,cz,cI,cJ,cK, dexpn,dx,dy,dz,dI,dJ,dK)
            @test hrr2_vals[(aI,aJ,aK, bI,bJ,bK, cI,cJ,cK, dI,dJ,dK)] == hrr_val
        end

        ashell,bshell,cshell,dshell = (2,1,2,0)
        @show MolecularIntegrals.hrr2(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)
        @show MolecularIntegrals.hrr3(ashell,bshell,cshell,dshell, aexpn,bexpn,cexpn,dexpn, A,B,C,D)   
    
    end
  
end
