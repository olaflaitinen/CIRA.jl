function sensitivity_matrix(model::Model, p::AbstractVector{<:Real})
    rows = length(output_trajectory(model, p))
    q = model.n_params
    S = Matrix{Float64}(undef, rows, q)
    for j in 1:q
        h = 1.0e-30
        pp = complex.(p)
        pp[j] += im * h
        yp = output_trajectory(model, pp)
        S[:, j] = imag.(yp) ./ h
    end
    return S
end
