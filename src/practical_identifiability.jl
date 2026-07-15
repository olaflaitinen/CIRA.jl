function fisher_information(S::AbstractMatrix{Float64}, sigma::Float64)
    return (transpose(S) * S) ./ (sigma * sigma)
end

function practical_scores(S::AbstractMatrix{Float64}, sigma::Float64)
    q = size(S, 2)
    scores = fill(Inf, q)
    r, _ = numerical_rank(S)
    if r < q
        return scores
    end
    F = fisher_information(S, sigma)
    Finv = inv(F)
    for j in 1:q
        d = Finv[j, j]
        scores[j] = d > 0.0 ? sqrt(d) : Inf
    end
    return scores
end
