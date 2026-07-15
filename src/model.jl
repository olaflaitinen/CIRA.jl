struct Model
    name::String
    discipline::String
    n_states::Int
    n_params::Int
    p_true::Vector{Float64}
    x0::Vector{Float64}
    f!::Function
    g::Function
    tspan::Tuple{Float64,Float64}
    n_times::Int
end
