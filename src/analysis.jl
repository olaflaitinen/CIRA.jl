struct IdentifiabilityResult
    model::String
    discipline::String
    n_states::Int
    n_params::Int
    rank::Int
    n_identifiable_functions::Int
    globally_identifiable::Bool
    parameter_flags::BitVector
    practical_scores::Vector{Float64}
    reparameterization_dimension::Int
end

function analyze(model::Model; sigma::Float64 = 0.05, seed::Int = 20260715)
    Random.seed!(seed)
    S = sensitivity_matrix(model, model.p_true)
    r, _ = numerical_rank(S)
    flags = parameter_identifiability(S, r)
    scores = practical_scores(S, sigma)
    global_id = r == model.n_params
    reparam = global_id ? 0 : r
    return IdentifiabilityResult(model.name, model.discipline, model.n_states,
        model.n_params, r, r, global_id, flags, scores, reparam)
end
