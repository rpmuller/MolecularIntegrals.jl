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
        @test MolecularIntegrals.nao[0] == 1
        @test MolecularIntegrals.nao[1] == 4
        @test MolecularIntegrals.nao[2] == 10

        @test MolecularIntegrals.m2ao[[0,0,0]] == 1
        @test MolecularIntegrals.m2ao[[0,0,1]] == 4

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

        # Move somewhere else?
        bfs = build_basis(h2)
        fetcher = MolecularIntegrals.eri_fetcher(bfs)
        @test_skip MolecularIntegrals.all_twoe_ints_chrr(bfs) ≈ MolecularIntegrals.all_twoe_ints(bfs)
        @test isapprox(MolecularIntegrals.all_twoe_ints(bfs),MolecularIntegrals.all_twoe_ints_rys(bfs),rtol=1e-7)
    
        # demonstrate that the normalization constants are not the same
        #  for all m-levels:
        #sh = Shell([0.0,0.0,0.0],2,[1.0],[1.0])
        #@show [pgbf.norm for pgbf in MolecularIntegrals.pgbfs(sh)]

    end

    @testset "OneInts" begin
        for (a,b,val) in [(0,0,1.0),(0,1,0.7468241328124271),(1,0,0.33333333333333)]
            @test MolecularIntegrals.Fgamma(a,b) ≈ val
        end

        @test kinetic(1.0, 0.0,0.0,0.0, 0,0,0, 1.0, 0.0,0.0,0.0, 0,0,0) ≈ 2.9530518648229536
        @test kinetic(s,s) ≈ 1.5
        @test kinetic(c,c) ≈ 1.5

        @test nuclear_attraction(s,s,[0.0,0.0,0.0]) ≈ -1.59576912
        @test nuclear_attraction(c,c,[0.0,0.0,0.0]) ≈ -1.59576912
        @test nuclear_attraction(px,px,[0.0,0.0,0.0]) ≈ -1.06384608
        @test nuclear_attraction(1,[0,0,0],0,0,0, 1,[0,0,0],0,0,0, [0,0,0]) ≈ -3.141592653589793

        @test overlap(s,s) ≈ 1
        @test overlap(px,px) ≈ 1
        @test overlap(s,px) ≈ 0
        @test overlap(c,c) ≈ 1

        @test overlap(1.0,[0.0,0.0,0.0],0,0,0, 1,[0.0,0.0,0.0],0,0,0) ≈ 1.9687012432

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
        # Tests from pyquante2:
        @test coulomb(1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0) ≈ 4.37335458
        @test coulomb(c,c,c,c) ≈ 1.128379167
        @test coulomb(s,s,px,px) ≈ 0.9403159725793302
        @test coulomb(s,s,s,px) == 0.0

        for (r,val) in [(0.0, 1.1283791670951362),
                        (1.0, 0.8427007900292186),
                        (2.0, 0.49766113257563993),
                        (3.0, 0.33332596983445223),
                        (4.0, 0.2499999961456855), # coulomb's law hereafter:
                        (5.0, 1/5), 
                        (6.0, 1/6),
                        (7.0, 1/7),
                        (8.0, 1/8),
                        (9.0, 1/9)]
            s3 = pgbf(1.0, 0.0,0.0,r)
            @test coulomb(s,s,s3,s3) ≈ val
        end
    end

    @testset "Rys tests" begin
        # Tests from pyquante2:
        @test coulomb_rys(c,c,c,c) ≈ 1.128379167
        @test coulomb_rys(s,s,px,px) ≈ 0.9403159725793302
        @test coulomb_rys(s,s,s,px) == 0.0

        for (r,val) in [(0.0, 1.1283791670951362),
                        (1.0, 0.8427007900292186),
                        (2.0, 0.49766113257563993),
                        (3.0, 0.33332596983445223),
                        (4.0, 0.2499999961456855), # coulomb's law hereafter:
                        (5.0, 1/5), 
                        (6.0, 1/6),
                        (7.0, 1/7),
                        (8.0, 1/8),
                        (9.0, 1/9)]
            s3 = pgbf(1.0, 0.0,0.0,r)
            @test coulomb_rys(s,s,s3,s3) ≈ val
        end
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
            vrrs = vrr(aI+aJ+aK, cI+cJ+cK, aexpn,bexpn,cexpn,dexpn, [ax,ay,az],[bx,by,bz],[cx,cy,cz],[dx,dy,dz])
            @test result ≈ vrrs[1,m2ao[[aI,aJ,aK]],m2ao[[cI,cJ,cK]]]
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
            ashell = aI+aJ+aK
            bshell = bI+bJ+bK
            cshell = cI+cJ+cK
            dshell = dI+dJ+dK
            hrrs = hrr(ashell,bshell,cshell,dshell,aexpn,bexpn,cexpn,dexpn,
                        [ax,ay,az],[bx,by,bz],[cx,cy,cz],[dx,dy,dz])
            #hrrds = MolecularIntegrals.hrr_dict(ashell,bshell,cshell,dshell,aexpn,bexpn,cexpn,dexpn,
            #            [ax,ay,az],[bx,by,bz],[cx,cy,cz],[dx,dy,dz])
            @test result ≈ hrrs[m2ao[[aI,aJ,aK]],m2ao[[bI,bJ,bK]],m2ao[[cI,cJ,cK]],m2ao[[dI,dJ,dK]]]
            #@test result ≈ hrrds[aI,aJ,aK,bI,bJ,bK,cI,cJ,cK,dI,dJ,dK]
        end
    end
  
end
