using MolecularIntegrals, Test

@testset "Basis" begin
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @test s(0,0,0) ≈ 0.71270547
    @test px(0,0,0) ≈ 0
    c = cgbf(0.0,0.0,0.0)
    #push!(c,1,1)
    #@test c(0,0,0) ≈ 0.71270547
    #c2 = cgbf(0,0,0)
    #push!(c2,1,0.2)
    #push!(c2,0.5,0.2)
    #@test overlap(c2,c2) ≈ 1
end
