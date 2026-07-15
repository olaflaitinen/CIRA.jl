# CIRA.jl

CIRA (Certified Identifiability and Reparameterization Analysis) is a Julia package for analyzing parameter identifiability of mechanistic ordinary differential equation models. This release is a dependency light reference implementation that relies only on the Julia standard library.

## Features

- Local structural identifiability through the numerical rank of the output sensitivity matrix.
- Per parameter identifiability flags derived from column drop tests.
- Practical identifiability scoring based on the Fisher information matrix.
- Eight benchmark models spanning epidemiology, systems biology, pharmacology, biochemistry, ecology and virology.
- Scripts that reproduce the benchmark table and the capability comparison table.

## Requirements

- Julia 1.9 or newer. No external packages are required.

## Installation

```
git clone https://github.com/olaflaitinen/CIRA.jl.git
cd CIRA.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Quick start

```
using CIRA

for model in benchmark_models()
    result = analyze(model)
    println(result.model, " rank=", result.rank, " global=", result.globally_identifiable)
end
```

## Reproducing the tables

```
julia --project=. scripts/reproduce_table1.jl
julia --project=. scripts/reproduce_table2.jl
```

## Running the tests

```
julia --project=. -e "using Pkg; Pkg.test()"
```

## Package structure

```
CIRA.jl/
  Project.toml
  LICENSE
  README.md
  CITATION.cff
  src/
    CIRA.jl
    model.jl
    integrator.jl
    sensitivity.jl
    local_identifiability.jl
    practical_identifiability.jl
    analysis.jl
    models.jl
  scripts/
    reproduce_table1.jl
    reproduce_table2.jl
  test/
    runtests.jl
```

## Scope and limitations

This reference release implements the local structural identifiability test and the practical identifiability scoring described in the accompanying manuscript. The identifiable function count and the reparameterization dimension are reported from the numerical rank of the sensitivity matrix. The symbolic global differential elimination and the explicit construction of a minimal reparameterization are planned for a future release. Numerical results depend on the integration grid, the observation design and the evaluation point.

## Citation

See CITATION.cff for citation metadata.

## License

MIT. See LICENSE.
