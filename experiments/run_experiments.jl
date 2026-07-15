using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CIRA
using LinearAlgebra
using Printf

const DATA = joinpath(@__DIR__, "data")
mkpath(DATA)

function write_csv(name, header, rows)
    open(joinpath(DATA, name), "w") do io
        println(io, header)
        for r in rows
            println(io, r)
        end
    end
end

function output_scale(m)
    y = CIRA.output_trajectory(m, m.p_true)
    return maximum(abs.(real.(y)))
end

function chain_model(k)
    f! = function (dx, x, p, t)
        dx[1] = -p[1] * x[1]
        for i in 2:k
            dx[i] = p[i - 1] * x[i - 1] - p[i] * x[i]
        end
        return nothing
    end
    g = function (x, p)
        return [x[k]]
    end
    p = [0.3 + 0.015 * (i - 1) for i in 1:k]
    x0 = zeros(k)
    x0[1] = 1.0
    return Model("chain-" * string(k), "Benchmark", k, k, p, x0, f!, g, (0.0, 10.0), 61)
end

function rebuild_grid(m, n_times)
    return Model(m.name, m.discipline, m.n_states, m.n_params, m.p_true, m.x0, m.f!, m.g, m.tspan, n_times)
end

function exp1_scaling()
    println("E1 runtime scaling")
    ks = [2, 4, 6, 8, 10, 14, 18, 24, 30, 36, 42, 50]
    warm = chain_model(4)
    CIRA.numerical_rank(CIRA.sensitivity_matrix(warm, warm.p_true))
    rows = String[]
    for k in ks
        m = chain_model(k)
        best = Inf
        for _ in 1:3
            t = @elapsed begin
                S = CIRA.sensitivity_matrix(m, m.p_true)
                CIRA.numerical_rank(S)
            end
            best = min(best, t)
        end
        push!(rows, @sprintf("%d,%d,%.6e", k, k, best))
    end
    write_csv("scaling.csv", "n_params,n_states,runtime_seconds", rows)
end

function exp2_practical(rel = 0.05)
    println("E2 practical identifiability")
    rows = String[]
    for m in benchmark_models()
        S = CIRA.sensitivity_matrix(m, m.p_true)
        r, _ = CIRA.numerical_rank(S)
        q = m.n_params
        sigma = rel * output_scale(m)
        F = (S' * S) ./ (sigma * sigma)
        if r == q
            Finv = inv(F)
            cv = sqrt.(abs.(diag(Finv))) ./ abs.(m.p_true)
            push!(rows, @sprintf("%s,%d,%d,yes,%.4f,%.4f", m.name, q, r, maximum(cv), sum(cv) / q))
        else
            push!(rows, @sprintf("%s,%d,%d,no,inf,inf", m.name, q, r))
        end
    end
    write_csv("practical.csv", "model,q,rank,all_params_finite_cv,max_cv,mean_cv", rows)
end

function exp3_reparam(rel = 0.05)
    println("E3 reparameterization validation")
    rows = String[]
    for m in benchmark_models()
        S = CIRA.sensitivity_matrix(m, m.p_true)
        r, _ = CIRA.numerical_rank(S)
        q = m.n_params
        sigma = rel * output_scale(m)
        F = (S' * S) ./ (sigma * sigma)
        fac = svd(S)
        cond_full = (fac.S[1] / fac.S[q])^2
        Vr = fac.V[:, 1:r]
        Fr = Vr' * F * Vr
        cond_red = cond(Fr)
        rr, _ = CIRA.numerical_rank(S * Vr)
        ok = rr == r ? "yes" : "no"
        push!(rows, @sprintf("%s,%d,%d,%d,%.3e,%.3e,%s", m.name, q, r, rr, cond_full, cond_red, ok))
    end
    write_csv("reparameterization.csv", "model,q,rank,reduced_full_rank,cond_full_fisher,cond_reduced_fisher,reduced_is_identifiable", rows)
end

function exp4_agreement()
    println("E4 verdict agreement")
    lit_verdict = Dict("SIR" => "identifiable", "SEIR" => "non-identifiable",
        "Goodwin oscillator" => "non-identifiable", "JAK-STAT signalling" => "non-identifiable",
        "Two-compartment PK (oral)" => "non-identifiable", "Michaelis-Menten kinetics" => "non-identifiable",
        "Lotka-Volterra" => "non-identifiable", "HIV viral dynamics" => "non-identifiable")
    lit_nif = Dict("SIR" => 2, "SEIR" => 3, "Goodwin oscillator" => 4, "JAK-STAT signalling" => 4,
        "Two-compartment PK (oral)" => 4, "Michaelis-Menten kinetics" => 2, "Lotka-Volterra" => 3,
        "HIV viral dynamics" => 5)
    lit_ref = Dict("SIR" => "Tuncer and Le 2018", "SEIR" => "Chowell 2017",
        "Goodwin oscillator" => "Villaverde et al. 2016", "JAK-STAT signalling" => "Swameye et al. 2003",
        "Two-compartment PK (oral)" => "Janzen et al. 2016", "Michaelis-Menten kinetics" => "Raue et al. 2009",
        "Lotka-Volterra" => "Remien et al. 2021", "HIV viral dynamics" => "Miao et al. 2011")
    rows = String[]
    for m in benchmark_models()
        res = analyze(m)
        verdict = res.globally_identifiable ? "identifiable" : "non-identifiable"
        agree = (verdict == lit_verdict[m.name]) && (res.rank == lit_nif[m.name])
        push!(rows, @sprintf("%s,%s,%d,%s,%d,%s,%s", m.name, verdict, res.rank,
            lit_verdict[m.name], lit_nif[m.name], lit_ref[m.name], agree ? "yes" : "no"))
    end
    write_csv("agreement.csv", "model,cira_verdict,cira_nif,literature_verdict,literature_nif,reference,agreement", rows)
end

function exp5_robustness(rel = 0.05)
    println("E5 robustness")
    factors = [0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
    rank_rows = String[]
    gap_rows = String[]
    for m in benchmark_models()
        base = m.n_times - 1
        grids = sort(unique([Int(round(base * f)) for f in factors]))
        for gs in grids
            m2 = rebuild_grid(m, gs + 1)
            S = CIRA.sensitivity_matrix(m2, m2.p_true)
            if !all(isfinite, S)
                continue
            end
            r, sv = CIRA.numerical_rank(S)
            smax = sv[1]
            small_id = sv[r] / smax
            large_null = r < length(sv) ? sv[r + 1] / smax : 0.0
            push!(rank_rows, @sprintf("%s,%d,%d", m.name, gs + 1, r))
            push!(gap_rows, @sprintf("%s,%d,%.3e,%.3e", m.name, gs + 1, small_id, large_null))
        end
    end
    write_csv("robustness_rank.csv", "model,n_times,rank", rank_rows)
    write_csv("robustness_gap.csv", "model,n_times,smallest_identifiable_sv,largest_null_sv", gap_rows)

    sig = [0.01, 0.02, 0.05, 0.1, 0.2]
    noise_rows = String[]
    m = sir_model()
    S = CIRA.sensitivity_matrix(m, m.p_true)
    scale = output_scale(m)
    for relv in sig
        sigma = relv * scale
        F = (S' * S) ./ (sigma * sigma)
        cv = sqrt.(abs.(diag(inv(F)))) ./ abs.(m.p_true)
        push!(noise_rows, @sprintf("%.3f,%.4f,%.4f", relv, cv[1], cv[2]))
    end
    write_csv("robustness_noise.csv", "relative_noise,cv_beta,cv_gamma", noise_rows)
end

function main()
    exp1_scaling()
    exp2_practical()
    exp3_reparam()
    exp4_agreement()
    exp5_robustness()
    println("data written to ", DATA)
end

main()
