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
        @test count(r.parameter_flags) <= r.rank
        @test r.globally_identifiable == (r.rank == m.n_params)
        @test r.reparameterization_dimension == (r.globally_identifiable ? 0 : r.rank)
    end

    results = Dict(m.name => analyze(m) for m in models)

    expected_rank = Dict(
        "SIR" => 2,
        "SEIR" => 3,
        "Goodwin oscillator" => 4,
        "JAK-STAT signalling" => 4,
        "Two-compartment PK (oral)" => 4,
        "Michaelis-Menten kinetics" => 2,
        "Lotka-Volterra" => 3,
        "HIV viral dynamics" => 5,
    )
    for (name, rk) in expected_rank
        @test results[name].rank == rk
    end

    @test results["SIR"].globally_identifiable
    @test results["SIR"].reparameterization_dimension == 0

    for name in ("SEIR", "Goodwin oscillator", "JAK-STAT signalling",
                 "Two-compartment PK (oral)", "Michaelis-Menten kinetics",
                 "Lotka-Volterra", "HIV viral dynamics")
        @test !results[name].globally_identifiable
        @test results[name].reparameterization_dimension == expected_rank[name]
    end

    mm = results["Michaelis-Menten kinetics"]
    @test mm.rank == 2
    @test !mm.globally_identifiable
    @test mm.reparameterization_dimension == 2

    goodwin = results["Goodwin oscillator"]
    @test goodwin.rank == 4
    @test !goodwin.globally_identifiable
    @test goodwin.reparameterization_dimension == 4
end
