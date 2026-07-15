function numerical_rank(S::AbstractMatrix{Float64})
    if size(S, 1) == 0 || size(S, 2) == 0
        return 0, Float64[]
    end
    factorization = svd(S)
    svals = factorization.S
    if length(svals) == 0 || svals[1] == 0.0
        return 0, svals
    end
    tol = maximum(size(S)) * eps(Float64) * svals[1]
    r = count(s -> s > tol, svals)
    return r, svals
end

function parameter_identifiability(S::AbstractMatrix{Float64}, full_rank::Int)
    q = size(S, 2)
    flags = falses(q)
    for j in 1:q
        cols = [c for c in 1:q if c != j]
        if length(cols) == 0
            flags[j] = full_rank > 0
        else
            rj, _ = numerical_rank(S[:, cols])
            flags[j] = rj < full_rank
        end
    end
    return flags
end
