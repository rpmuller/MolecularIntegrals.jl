using MolecularIntegrals, Test

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
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @test s(0,0,0) ≈ 0.71270547
    @test px(0,0,0) ≈ 0
    c = cgbf(0.0,0.0,0.0)
    addbf!(c,1,1)
    @test c(0,0,0) ≈ 0.71270547
    c2 = cgbf(0,0,0)
    addbf!(c2,1,0.2)
    addbf!(c2,0.5,0.2)
    @test overlap(c2,c2) ≈ 1
end

@testset "One Electron Integrals" begin
    @test kinetic(1.0, 0.0,0.0,0.0, 0,0,0, 1.0, 0.0,0.0,0.0, 0,0,0) ≈ 2.9530518648229536
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @test kinetic(s,s) ≈ 1.5
    c = cgbf(0.0,0.0,0.0)
    addbf!(c,1,1)
    @test kinetic(c,c) ≈ 1.5
    @test nuclear_attraction(s,s,0.0,0.0,0.0) ≈ -1.59576912
    @test nuclear_attraction(c,c,0.0,0.0,0.0) ≈ -1.59576912

    @test overlap1d(0,0,0.0,0.0,1.0) == 1
    @test gaussian_product_center(s,s) == [0,0,0]
    @test overlap(s,s) ≈ 1
    @test overlap(px,px) ≈ 1
    @test overlap(s,px) ≈ 0
    @test binomial_prefactor(0,0,0,0.0,0.0) == 1
    @test overlap(c,c) ≈ 1
end
