function sir_model()
    f! = function (dx, x, p, t)
        beta = p[1]
        gamma = p[2]
        S = x[1]
        I = x[2]
        dx[1] = -beta * S * I
        dx[2] = beta * S * I - gamma * I
        dx[3] = gamma * I
        return nothing
    end
    g = function (x, p)
        return [x[2]]
    end
    return Model("SIR", "Epidemiology", 3, 2, [0.9, 0.3], [0.99, 0.01, 0.0], f!, g, (0.0, 20.0), 41)
end

function seir_model()
    f! = function (dx, x, p, t)
        beta = p[1]
        v = p[2]
        psi = p[3]
        gamma = p[4]
        S = x[1]
        E = x[2]
        I = x[3]
        Q = x[5]
        dx[1] = -beta * S * I
        dx[2] = beta * S * I - v * E
        dx[3] = v * E - psi * I - (1.0 - psi) * gamma * I
        dx[4] = gamma * Q + (1.0 - psi) * gamma * I
        dx[5] = -gamma * Q + psi * I
        return nothing
    end
    g = function (x, p)
        return [x[5]]
    end
    return Model("SEIR_1_io", "Epidemiology", 5, 4, [0.9, 0.4, 0.3, 0.6], [0.99, 0.0, 0.01, 0.0, 0.0], f!, g, (0.0, 20.0), 41)
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
        dx[2] = c * X - d * Y
        dx[3] = e * Y - k * Z
        return nothing
    end
    g = function (x, p)
        return [x[1]]
    end
    return Model("Goodwin oscillator", "Systems biology", 3, 6, [1.0, 0.3, 0.8, 0.3, 0.8, 0.3], [0.1, 0.1, 0.1], f!, g, (0.0, 20.0), 41)
end

function jak_stat_model()
    f! = function (dx, x, p, t)
        u = 1.0
        t1 = p[1]
        t2 = p[2]
        t3 = p[3]
        t4 = p[4]
        t5 = p[5]
        t6 = p[6]
        t7 = p[7]
        t8 = p[8]
        t9 = p[9]
        t10 = p[10]
        t11 = p[11]
        t12 = p[12]
        t13 = p[13]
        t14 = p[14]
        t15 = p[15]
        t16 = p[16]
        t17 = p[17]
        t18 = p[18]
        t19 = p[19]
        t20 = p[20]
        t21 = p[21]
        t22 = p[22]
        x1 = x[1]
        x2 = x[2]
        x3 = x[3]
        x4 = x[4]
        x5 = x[5]
        x6 = x[6]
        x7 = x[7]
        x8 = x[8]
        x9 = x[9]
        x10 = x[10]
        dx[1] = t6 * x2 - t5 * x1 - 2.0 * t1 * u * x1
        dx[2] = -t6 * x2 + t5 * x1
        dx[3] = x6 * x3 * t2 - 3.0 * x3 * t2 + 2.0 * t1 * u * x1
        dx[4] = -t3 * x4 - x6 * x3 * t2 + 3.0 * x3 * t2
        dx[5] = t3 * x4 - x5 * t4
        dx[6] = (-x6 * x3 * x10 * t7 * t13 - x6 * x3 * t7 - 92.0 * x6 * x10 * x1 * t8 * t13^2 - 92.0 * x6 * x10 * t8 * t13 - 92.0 * x6 * x1 * t8 * t13 - x6 * x1 * t7 * t13 * x4 - 92.0 * x6 * t8 - x6 * t7 * x4 + 276.0 * x10 * x1 * t8 * t13^2 + 276.0 * x10 * t8 * t13 + 276.0 * x1 * t8 * t13 + 276.0 * t8) / (x10 * x1 * t13^2 + x10 * t13 + x1 * t13 + 1.0)
        dx[7] = -92.0 * x7 * t10 + x7 * x6 * t9 - 3.0 * x7 * t9 + 15180.0 * t10
        dx[8] = -x7 * t11 + 165.0 * t11
        dx[9] = -2.0 * x9 * u * t12
        dx[10] = (-x8 * t16 * x10 + x8 * t14 - t16 * x10 * t15) / (x8 + t15)
        return nothing
    end
    g = function (x, p)
        t11 = p[11]
        t17 = p[17]
        t18 = p[18]
        t19 = p[19]
        t20 = p[20]
        t21 = p[21]
        t22 = p[22]
        return [
            x[3] + x[1] + x[4],
            -x[9] * t18 + x[5] * t18 + t18 * x[3] + t18 * x[4] + (1.0 / 3.0) * t18,
            t19 * x[5] + t19 * x[4],
            -t20 * x[6] + 3.0 * t20,
            x[8] * t21,
            (x[8] * t22 * t17) / t11,
            x[10],
            -x[7] + 165.0,
        ]
    end
    return Model("JAK-STAT 1", "Systems biology", 10, 22, [
        0.2, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2,
    ], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 165.0, 1.0], f!, g, (0.0, 15.0), 31)
end

function pk_model()
    f! = function (dx, x, p, t)
        ka = p[1]
        kc = p[2]
        a1 = p[3]
        a2 = p[4]
        b1 = p[5]
        b2 = p[6]
        n = p[7]
        x0 = x[1]
        x1 = x[2]
        x2 = x[3]
        x3 = x[4]
        denom = ka * kc + ka * x0 + kc * x2
        dx[1] = (-ka * n * x0 - ka * kc * a1 * x0 + ka * kc * a1 * x1 - ka * a1 * x0^2 + ka * a1 * x0 * x1 - kc * a1 * x0 * x2 + kc * a1 * x1 * x2) / denom
        dx[2] = a2 * x0 - a2 * x1
        dx[3] = (ka * kc * b1 * x3 - ka * kc * b1 * x2 + ka * b1 * x0 * x3 - ka * b1 * x0 * x2 - n * kc * x2 + kc * b1 * x3 * x2 - kc * b1 * x2^2) / denom
        dx[4] = -b2 * x3 + b2 * x2
        return nothing
    end
    g = function (x, p)
        return [x[1]]
    end
    return Model("Pharm", "Pharmacology", 4, 7, [1.0, 1.0, 0.5, 0.3, 0.2, 0.4, 10.0], [5.0, 0.0, 0.0, 0.0], f!, g, (0.0, 24.0), 49)
end

function michaelis_menten_model()
    f! = function (dx, x, p, t)
        kcat = p[1]
        Km = p[2]
        E0 = p[3]
        S = x[1]
        v = kcat * E0 * S / (Km + S)
        dx[1] = -v
        dx[2] = v
        return nothing
    end
    g = function (x, p)
        return [x[2]]
    end
    return Model("Michaelis-Menten kinetics", "Biochemistry", 2, 3, [2.0, 1.0, 1.5], [5.0, 0.0], f!, g, (0.0, 10.0), 41)
end

function lotka_volterra_model()
    f! = function (dx, x, p, t)
        a = p[1]
        b = p[2]
        c = p[3]
        d = p[4]
        x1 = x[1]
        x2 = x[2]
        dx[1] = (a + b) * x1 - c * x1 * x2
        dx[2] = -a * b * x2 + d * x1 * x2
        return nothing
    end
    g = function (x, p)
        return [x[1]]
    end
    return Model("Modified LV for testing", "Ecology", 2, 4, [1.0, 2.0, 0.4, 0.6], [1.0, 0.5], f!, g, (0.0, 15.0), 31)
end

function hiv_model()
    f! = function (dx, x, p, t)
        lm = p[1]
        d = p[2]
        beta = p[3]
        a = p[4]
        k = p[5]
        u = p[6]
        b = p[7]
        c = p[8]
        q = p[9]
        h = p[10]
        x_state = x[1]
        y_state = x[2]
        v_state = x[3]
        w_state = x[4]
        z_state = x[5]
        dx[1] = lm - x_state * d - x_state * v_state * beta
        dx[2] = x_state * v_state * beta - a * y_state
        dx[3] = k * y_state - v_state * u
        dx[4] = -b * w_state + c * w_state * x_state * y_state - c * w_state * q * y_state
        dx[5] = c * w_state * q * y_state - h * z_state
        return nothing
    end
    g = function (x, p)
        return [x[4], x[5]]
    end
    return Model("HIV", "Virology", 5, 10, [10.0, 0.1, 0.05, 0.5, 5.0, 3.0, 0.7, 0.8, 0.9, 1.1], [100.0, 0.0, 0.001, 0.0, 0.0], f!, g, (0.0, 10.0), 41)
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
