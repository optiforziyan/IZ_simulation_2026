# From Multidimensional Biodiversity Patterns to Conservation Features in Systematic Conservation Planning

R scripts and synthetic data supporting the illustrative numerical example in the correspondence submitted to *Integrative Zoology* (manuscript INZ-2026-08-650).

The experiment distinguishes local biodiversity scores from explicit conservation features. It compares species representation targets, aggregated-index targets, and integrated planning formulations that retain hard species targets while using local biodiversity scores as spatial preferences.

## Simulation overview

- A 24 × 24 landscape of 576 planning units, each assigned a conventional cost of 1, with two existing reserves locked into every solution.
- Eight synthetic species with clustered ranges: three widespread, three intermediate, and two narrowly distributed species.
- Range-based species targets of 20%, 30%, or 50%, depending on the proportion of the landscape occupied.
- Local taxonomic diversity (TD), phylogenetic diversity (PD), functional diversity (FD), and a synthetic genetic-diversity proxy (GD).
- A local preference score combining min–max-normalized PD, FD, and GD, with equal weights in the reference analysis. Its sum over selected units is Q; TD is excluded from this score because species targets already specify taxonomic representation requirements.
- Sensitivity analyses for dimensional weights, penalty strengths, cost tolerances, and species-target policies, plus ten landscape seeds (20260925–20260934). The main example uses seed 20260925.

All species distributions, tree branches, functional entities, and genetic-diversity proxies are synthetic. Files named `*_SDM.tif` contain simulated binary occurrences, not fitted empirical species distribution models. The 19 within-species spatial strata are illustrative subdivisions, not empirically inferred populations or genetic lineages.

## Repository contents

| Path | Purpose |
| --- | --- |
| `code/run_simulation.R` | Generate landscapes, solve planning scenarios, run sensitivity analyses and repeated simulations, and export results and figures. |
| `code/audit_results.R` | Check cached solutions for binary decisions, representation targets, reserve lock-in, area limits, objective values, and solver gaps. |
| `code/make_concept_figure.R` | Draw a separate conceptual illustration of local summaries and explicit entity targets. |
| `data/` | Synthetic rasters, feature matrices, parameters, toy phylogeny, and serialized landscape. |
| `data/SDM/` | Simulated suitability and binary occurrence rasters for the eight species. |
| `figures/` | Simulation figures supplied with the repository. |
| `results/` | Created by running the simulation; contains numerical outputs, solutions, checkpoints, and software information. |

## Software requirements

Use R 4.1 or later (the script uses the native pipe operator) and the R packages `terra`, `prioritizr`, and `highs`. HiGHS is the optimization solver. Install the packages in R:

```r
install.packages(c("terra", "prioritizr", "highs"))
```

The simulation records the R session, package versions, HiGHS version, and random-number generator settings in `results/sessionInfo.txt`. The repository does not currently include a package-version lockfile; retain this session record when reproducing or extending the analysis.

## Running the analysis

Download or clone the repository. From its root directory, run the following commands in a terminal with `Rscript` available:

```sh
# Main example and sensitivity analyses
Rscript code/run_simulation.R main

# Ten landscape seeds
Rscript code/run_simulation.R replicates

# Alternatively, run both in one invocation
Rscript code/run_simulation.R all

# Audit saved solutions after a simulation run
Rscript code/audit_results.R

# Optional conceptual illustration
Rscript code/make_concept_figure.R
```

The simulation locates the repository root from the script path and creates the output directories automatically. Running the main analysis regenerates synthetic data and figures. Saved checkpoints are reused when the checked scenario inputs match. For a fresh run after changing code, move existing checkpoints out of `results/checkpoints/` before rerunning. Run the audit only after generating checkpoints.

## Outputs and interpretation

Main outputs include `main_comparison.csv`, `species_representation.csv`, `weight_epsilon_sensitivity.csv`, and `target_sensitivity.csv`. Repeated simulations produce `replicate_results.csv`. Solver diagnostics are saved to `solver_log.csv`; the audit writes `solver_audit.csv`. Output files are written under `results/`.

The integrated effective-cost and additive soft-penalty formulations are algebraically equivalent under the same constraints when the dimension-specific penalty coefficient is λ_d = ρw_d. Here, ρ controls overall penalty strength and w_d is the normalized dimensional preference weight. The hierarchical formulation first minimizes conventional cost and then maximizes Q within a proportional cost tolerance ε.

Higher Q indicates a larger sum of local preference scores. It does not by itself guarantee greater unique phylogenetic or functional coverage, or better representation of within-species strata. Explicit hard targets guarantee only the stated requirements for the features included in the optimization. See the manuscript and Supplementary Materials for the equations, scenario definitions, and interpretation.

## Associated manuscript

*From Multidimensional Biodiversity Patterns to Conservation Features in Systematic Conservation Planning.* Correspondence submitted to *Integrative Zoology*, manuscript INZ-2026-08-650.

Repository: https://github.com/optiforziyan/IZ_simulation_2026
