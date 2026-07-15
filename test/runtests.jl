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
    @test !analyze(seir_model()).globally_identifiable
    @test !analyze(jak_stat_model()).globally_identifiable
    @test !analyze(pk_model()).globally_identifiable
    @test !analyze(lotka_volterra_model()).globally_identifiable
    @test !analyze(hiv_model()).globally_identifiable
    mm = analyze(michaelis_menten_model())
    @test mm.rank <= 2
    @test !mm.globally_identifiable
end
