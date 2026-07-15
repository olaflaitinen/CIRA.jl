using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CIRA

function main()
    caps, tools, matrix = capability_matrix()
    header = rpad("Capability", 40)
    for t in tools
        header = header * " | " * rpad(t, 28)
    end
    println(header)
    for i in eachindex(caps)
        line = rpad(caps[i], 40)
        for j in eachindex(tools)
            line = line * " | " * rpad(matrix[i, j], 28)
        end
        println(line)
    end
end

main()
