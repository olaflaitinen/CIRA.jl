using Test
using CIRA

@testset "CIRA benchmarks" begin
    models = benchmark_models()
    @test length(models) == 8
    for m in models
        r = analyze(m)
        @test r.rank >= 0
        @test r.rank <= m.n_params
        @test length(r.parameter_flags) == m.n_params
        @test length(r.practical_scores) == m.n_params
        @test r.n_identifiable_functions == r.rank
    end
    mm = analyze(michaelis_menten_model())
    @test mm.rank == 3
    @test mm.globally_identifiable
end
