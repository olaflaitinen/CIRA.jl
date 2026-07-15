function sir_model()
    f! = function (dx, x, p, t)
        beta = p[1]
        gamma = p[2]
        S = x[1]
        I = x[2]
        dx[1] = -beta * S * I
        dx[2] = beta * S * I - gamma * I
        dx[3] = beta * S * I
        return nothing
    end
    g = function (x, p)
        return [x[3]]
    end
    return Model("SIR", "Epidemiology", 3, 2, [0.9, 0.3], [0.99, 0.01, 0.0], f!, g, (0.0, 20.0), 41)
end

function seir_model()
    f! = function (dx, x, p, t)
        beta = p[1]
        sigma = p[2]
        gamma = p[3]
        rho = p[4]
        tau = beta * rho
        S = x[1]
        E = x[2]
        I = x[3]
        dx[1] = -tau * S * I
        dx[2] = tau * S * I - sigma * E
        dx[3] = sigma * E - gamma * I
        dx[4] = tau * S * I
        return nothing
    end
    g = function (x, p)
        return [x[4]]
    end
    return Model("SEIR", "Epidemiology", 4, 4, [0.8, 0.5, 0.3, 0.6], [0.99, 0.0, 0.01, 0.0], f!, g, (0.0, 20.0), 41)
end

function goodwin_model()
    f! = function (dx, x, p, t)
        a = p[1]
        b = p[2]
        c = p[3]
        d = p[4]
        e = p[5]
        k = p[6]
        X = x[1]
        Y = x[2]
        Z = x[3]
        dx[1] = a / (1.0 + Z^4) - b * X
        dx[2] = (c * d) * X - 0.4 * Y
        dx[3] = (e * k) * Y - 0.35 * Z
        return nothing
    end
    g = function (x, p)
        return [x[1]]
    end
    return Model("Goodwin oscillator", "Systems biology", 3, 6, [1.0, 0.3, 0.5, 0.4, 0.6, 0.35], [0.3, 0.2, 0.1], f!, g, (0.0, 20.0), 41)
end

function jak_stat_model()
    f! = function (dx, x, p, t)
        k1 = p[1]
        k2 = p[2]
        x1 = x[1]
        x2 = x[2]
        x3 = x[3]
        x4 = x[4]
        dx[1] = -k1 * x1
        dx[2] = k1 * x1 - k2 * x2
        dx[3] = k2 * x2 - 0.4 * x3
        dx[4] = 0.4 * x3 - 0.3 * x4
        return nothing
    end
    g = function (x, p)
        s1 = p[3]
        a1 = p[4]
        s2 = p[5]
        a2 = p[6]
        return [(s1 * a1) * (x[2] + 2.0 * x[3]), (s2 * a2) * (x[1] + x[3])]
    end
    return Model("JAK-STAT signalling", "Systems biology", 4, 6, [0.8, 0.5, 1.2, 0.9, 1.1, 0.7], [1.0, 0.0, 0.0, 0.0], f!, g, (0.0, 15.0), 31)
end

function pk_model()
    f! = function (dx, x, p, t)
        ka = p[1]
        ke = p[2]
        k12 = p[3]
        Ad = x[1]
        Ac = x[2]
        Ap = x[3]
        dx[1] = -ka * Ad
        dx[2] = ka * Ad - ke * Ac - k12 * Ac + 0.4 * Ap
        dx[3] = k12 * Ac - 0.4 * Ap
        return nothing
    end
    g = function (x, p)
        F = p[4]
        V = p[5]
        return [(F / V) * x[2]]
    end
    return Model("Two-compartment PK (oral)", "Pharmacology", 3, 5, [1.2, 0.4, 0.5, 0.8, 3.0], [1.0, 0.0, 0.0], f!, g, (0.0, 24.0), 49)
end

function michaelis_menten_model()
    f! = function (dx, x, p, t)
        kcat = p[1]
        Km = p[2]
        E0 = p[3]
        S = x[1]
        v = (kcat * E0) * S / (Km + S)
        dx[1] = -v
        dx[2] = v
        return nothing
    end
    g = function (x, p)
        return [x[2]]
    end
    return Model("Michaelis-Menten kinetics", "Biochemistry", 2, 3, [3.0, 0.5, 2.0], [1.0, 0.0], f!, g, (0.0, 10.0), 41)
end

function lotka_volterra_model()
    f! = function (dx, x, p, t)
        alpha = p[1]
        beta = p[2]
        delta = p[3]
        gamma = p[4]
        rate = beta * delta
        X = x[1]
        Y = x[2]
        dx[1] = alpha * X - rate * X * Y
        dx[2] = rate * X * Y - gamma * Y
        return nothing
    end
    g = function (x, p)
        return [x[1]]
    end
    return Model("Lotka-Volterra", "Ecology", 2, 4, [1.0, 0.6, 0.4, 0.5], [1.0, 0.5], f!, g, (0.0, 15.0), 31)
end

function hiv_model()
    f! = function (dx, x, p, t)
        lambda = p[1]
        d = p[2]
        beta = p[3]
        delta = p[4]
        prod = p[5]
        c = p[6]
        T = x[1]
        I = x[2]
        V = x[3]
        dx[1] = lambda - d * T - 0.001 * T * V
        dx[2] = beta * T * V - delta * I
        dx[3] = prod * I - c * V
        return nothing
    end
    g = function (x, p)
        return [x[3]]
    end
    return Model("HIV viral dynamics", "Virology", 3, 6, [10.0, 0.1, 0.002, 0.5, 50.0, 3.0], [100.0, 0.0, 1.0], f!, g, (0.0, 10.0), 41)
end

function benchmark_models()
    return [sir_model(), seir_model(), goodwin_model(), jak_stat_model(),
        pk_model(), michaelis_menten_model(), lotka_volterra_model(), hiv_model()]
end

function capability_matrix()
    tools = ["DAISY", "GenSSI 2.0", "SIAN", "StructuralIdentifiability.jl",
        "STRIKE-GOLDD", "COMBOS", "Profile likelihood", "CIRA"]
    capabilities = ["Global structural verdict", "Identifiable combinations",
        "Automatic minimal reparameterization", "Practical identifiability",
        "Probabilistic correctness guarantee"]
    matrix = [
        "Yes" "Yes" "Yes" "Yes" "Yes" "Yes" "No" "Yes";
        "Partial" "Partial" "No" "Yes" "Partial" "Yes" "No" "Yes";
        "No" "No" "No" "No" "Partial" "Partial" "No" "Yes";
        "No" "No" "No" "No" "No" "No" "Yes" "Yes";
        "No" "No" "Yes" "Yes" "No" "No" "No" "Yes"
    ]
    return capabilities, tools, matrix
end
