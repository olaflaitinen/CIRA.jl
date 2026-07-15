using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CIRA

const TABLE1_NAMES = Dict(
    "SIR" => "SIR",
    "SEIR" => "SEIR",
    "Goodwin oscillator" => "Goodwin oscillator",
    "JAK-STAT signalling" => "JAK-STAT",
    "Two-compartment PK (oral)" => "PK (two-compartment)",
    "Michaelis-Menten kinetics" => "Michaelis-Menten",
    "Lotka-Volterra" => "Lotka-Volterra",
    "HIV viral dynamics" => "HIV",
)

function main()
    models = benchmark_models()
    header = rpad("Model", 28) * rpad("q", 5) * rpad("rank", 6) * rpad("global", 8) * rpad("reparam", 8)
    println(header)
    for m in models
        r = analyze(m)
        display_name = get(TABLE1_NAMES, r.model, r.model)
        line = rpad(display_name, 28) * rpad(string(r.n_params), 5) * rpad(string(r.rank), 6) * rpad(r.globally_identifiable ? "yes" : "no", 8) * rpad(string(r.reparameterization_dimension), 8)
        println(line)
    end
end

main()
