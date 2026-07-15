using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CIRA

function main()
    models = benchmark_models()
    header = rpad("Model", 28) * rpad("Discipline", 18) * rpad("n", 5) * rpad("q", 5) * rpad("rank", 6) * rpad("nif", 5) * rpad("global", 8) * rpad("reparam", 8)
    println(header)
    for m in models
        r = analyze(m)
        line = rpad(r.model, 28) * rpad(r.discipline, 18) * rpad(string(r.n_states), 5) * rpad(string(r.n_params), 5) * rpad(string(r.rank), 6) * rpad(string(r.n_identifiable_functions), 5) * rpad(r.globally_identifiable ? "yes" : "no", 8) * rpad(string(r.reparameterization_dimension), 8)
        println(line)
    end
end

main()
