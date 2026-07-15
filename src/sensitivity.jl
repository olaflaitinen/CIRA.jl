function sensitivity_matrix(model::Model, p::Vector{Float64})
    base = output_trajectory(model, p)
    rows = length(base)
    q = model.n_params
    S = Matrix{Float64}(undef, rows, q)
    h = 1.0e-20
    for j in 1:q
        pc = ComplexF64.(p)
        pc[j] = pc[j] + im * h
        yc = output_trajectory(model, pc)
        S[:, j] = imag.(yc) ./ h
    end
    return S
end
