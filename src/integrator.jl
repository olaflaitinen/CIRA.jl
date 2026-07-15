function rk4_step(f!, x, p, t, dt)
    k1 = similar(x)
    f!(k1, x, p, t)
    k2 = similar(x)
    f!(k2, x .+ 0.5 .* dt .* k1, p, t + 0.5 * dt)
    k3 = similar(x)
    f!(k3, x .+ 0.5 .* dt .* k2, p, t + 0.5 * dt)
    k4 = similar(x)
    f!(k4, x .+ dt .* k3, p, t + dt)
    return x .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

function integrate(model::Model, p::Vector{Float64})
    t0, t1 = model.tspan
    n = model.n_times
    dt = (t1 - t0) / (n - 1)
    x = copy(model.x0)
    xs = Vector{Vector{Float64}}(undef, n)
    ts = Vector{Float64}(undef, n)
    xs[1] = copy(x)
    ts[1] = t0
    for i in 2:n
        t = t0 + (i - 2) * dt
        x = rk4_step(model.f!, x, p, t, dt)
        xs[i] = copy(x)
        ts[i] = t0 + (i - 1) * dt
    end
    return ts, xs
end

function output_trajectory(model::Model, p::Vector{Float64})
    ts, xs = integrate(model, p)
    g0 = model.g(xs[1], p)
    m = length(g0)
    vals = Vector{Float64}(undef, length(ts) * m)
    idx = 1
    for i in eachindex(ts)
        gi = model.g(xs[i], p)
        for j in 1:m
            vals[idx] = gi[j]
            idx += 1
        end
    end
    return vals
end
