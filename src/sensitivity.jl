function sensitivity_matrix(model::Model, p::Vector{Float64})
    base = output_trajectory(model, p)
    rows = length(base)
    q = model.n_params
    S = Matrix{Float64}(undef, rows, q)
    for j in 1:q
        h = 1.0e-6 * (1.0 + abs(p[j]))
        pp = copy(p)
        pm = copy(p)
        pp[j] += h
        pm[j] -= h
        yp = output_trajectory(model, pp)
        ym = output_trajectory(model, pm)
        S[:, j] = (yp .- ym) ./ (2.0 * h)
    end
    return S
end
