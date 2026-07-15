# CIRA.jl experiments

This folder reproduces the numerical experiments reported in the manuscript. All
computations use only the CIRA package and the Julia standard library. Figures
are rendered by a small Python script from the data files produced by the Julia
run, so the plotting dependency is kept out of the package itself.

## Experiments

1. Runtime scaling. Analysis runtime as a function of the number of parameters,
   measured on a family of linear chain models. Produces `data/scaling.csv`.
2. Practical identifiability. Fisher based relative standard errors at five
   percent measurement noise. Produces `data/practical.csv`.
3. Reparameterization validation. Fisher conditioning of the full parameter set
   compared with the identifiable subspace recovered by CIRA. Produces
   `data/reparameterization.csv`.
4. Verdict agreement. CIRA structural verdicts and identifiable function counts
   compared with published results for each benchmark. Produces
   `data/agreement.csv`.
5. Robustness. Rank verdict stability and identifiable spectral gap across
   integration grids, and parameter uncertainty across noise levels. Produces
   `data/robustness_rank.csv`, `data/robustness_gap.csv`, and
   `data/robustness_noise.csv`.

## Figures

The plotting script renders three multi-panel figures, named with the journal
numbering convention and written as both high-resolution `png` (400 dpi) and
vector `pdf`:

- `figures/Fig1` (scaling). Experiment 1.
- `figures/Fig2` (identifiability). Experiment 2, experiment 3, and the noise
  sweep from experiment 5, as panels a, b, c.
- `figures/Fig3` (grid robustness). The rank stability and spectral gap results
  from experiment 5, as panels a and b.

Figure descriptions are given in `figure_captions.md`; following the target
journal style, they are not printed inside the figures.

## Reproducing

```
julia --project=. experiments/run_experiments.jl
python3 experiments/plot_figures.py
```

The first command writes the data files under `experiments/data`. The second
command reads those files and writes `png` and `pdf` figures under
`experiments/figures`. The Python script requires numpy and matplotlib.

## Note on runtime figures

Absolute runtimes depend on the machine and the language runtime. The reported
scaling exponent, obtained from a power law fit, is the reproducible quantity
and is not sensitive to these factors.
