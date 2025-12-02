mutenm package
================
Julian Echave

<!-- README.md is generated from README.Rmd. Please edit that file -->

The `mutenm` package builds Elastic Network Models (ENMs) of proteins
and simulates mutation effects using linear response theory.

# Installation

    # install.packages("devtools")
    devtools::install_github("jechave/mutenm")

# Overview

## What is an ENM?

An Elastic Network Model represents a protein as a network of nodes
(typically Cα atoms) connected by harmonic springs. The model captures
the protein’s intrinsic flexibility through normal mode analysis,
producing a covariance matrix that describes thermal fluctuations.

## How are mutations simulated?

Mutations are simulated by perturbing the equilibrium lengths of springs
connected to the mutated site. The structural response is calculated
using linear response theory:

    δr = C × f

where `C` is the covariance matrix and `f` is the force from edge
perturbations.

Two mutation models are available:

- **lfenm**: Fast. Updates structure only; dynamics unchanged.
- **sclfenm**: Slower. Recalculates the full ENM, capturing changes in
  dynamics and entropy.

## What can you calculate?

- **Structure**: contact number, weighted contact number, distance to
  active site
- **Dynamics**: mean-square fluctuations per site or per mode
- **Energy**: minimum (stress) energy, entropic free energy, activation
  energies
- **Mutation effects**: changes in structure, dynamics, and energy
  between wild-type and mutant

# Key Functions

## Core

| Function                       | Description                                 |
|--------------------------------|---------------------------------------------|
| `enm(pdb, node, model, d_max)` | Build an ENM from a PDB structure           |
| `mutenm(wt, site_mut, ...)`    | Simulate a mutation at a specific site      |
| `mrs(wt, nmut, ...)`           | Mutation Response Scanning — scan all sites |

## Single-protein properties

| Function                                         | Description                      |
|--------------------------------------------------|----------------------------------|
| `cn(prot)`                                       | Contact number per site          |
| `wcn(prot)`                                      | Weighted contact number per site |
| `dactive(prot, pdb_site_active)`                 | Distance to active site          |
| `msfi(prot)`                                     | Mean-square fluctuation per site |
| `msfn(prot)`                                     | Mean-square fluctuation per mode |
| `v_min(prot)`                                    | Minimum (stress) energy          |
| `g_ent(prot, beta)`                              | Entropic free energy             |
| `dv_act(prot, ideal, pdb_site_active)`           | Activation energy (internal)     |
| `dg_ent_act(prot, ideal, pdb_site_active, beta)` | Activation energy (entropic)     |

## Mutation effects (comparing wild-type and mutant)

| Function                    | Description                               |
|-----------------------------|-------------------------------------------|
| `Dr2i(wt, mut)`             | Squared displacement per site             |
| `Dr2n(wt, mut)`             | Squared displacement per mode             |
| `Dmsfi(wt, mut)`            | Change in MSF per site (requires sclfenm) |
| `Dmsfn(wt, mut)`            | Change in MSF per mode (requires sclfenm) |
| `Dv_min(wt, mut)`           | Change in minimum energy                  |
| `Dg_ent(wt, mut, beta)`     | Change in entropic free energy            |
| `Ddv_act(wt, mut, ...)`     | Change in activation energy (internal)    |
| `Ddg_ent_act(wt, mut, ...)` | Change in activation energy (entropic)    |

# Quick Start

``` r
library(mutenm)
```

## 1. Build an ENM

Read a PDB file and create the elastic network model:

``` r
pdb <- bio3d::read.pdb("data-raw/2acy.pdb")
wt <- enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5)
```

The `wt` object contains the protein structure, network topology, and
normal modes.

Calculate structural and dynamic properties:

``` r
# Contact number per site
head(cn(wt))
#> [1] 15 14 14 13 19 24

# Mean-square fluctuation per site
head(msfi(wt))
#> [1] 0.27409368 0.22327783 0.24816286 0.20781144 0.15790654 0.08517998

# Distance to active site (residues 23 and 41 for 2acy)
head(dactive(wt, pdb_site_active = c(23, 41)))
#> [1] 26.57028 25.64200 26.40024 26.09031 26.31486 23.36926
```

## 2. Simulate a mutation

Create a mutant by perturbing contacts at a specific site:

``` r
mut <- mutenm(wt, site_mut = 50, mutation = 1, mut_model = "lfenm", seed = 123)
```

## 3. Analyze mutation effects

Compare wild-type and mutant properties:

``` r
# Energy change (stress energy)
dv <- Dv_min(wt, mut)
dv
#> [1] 5.58627

# Structural displacement at each site
dr2 <- Dr2i(wt, mut)
head(dr2)
#> [1] 0.0042887367 0.0018992079 0.0024924585 0.0015339962 0.0012536041
#> [6] 0.0001890627
```

## 4. Mutation Response Scanning

Scan all sites systematically with `mrs()`:

``` r
# Scan all sites with 3 mutations each
responses <- mrs(wt, nmut = 3, mut_model = "lfenm", seed = 42)

# Response tibbles
names(responses)
#> [1] "Dr2i"        "Dmsfi"       "Dr2n"        "Dmsfn"       "Dv_min"     
#> [6] "Dg_ent"      "Ddv_act"     "Ddg_ent_act" "params"

# Energy change per mutation site
head(responses$Dv_min)
#> # A tibble: 6 × 2
#>       j Dv_min
#>   <int>  <dbl>
#> 1     1   1.51
#> 2     2   1.54
#> 3     3   1.86
#> 4     4   1.52
#> 5     5   2.92
#> 6     6   4.32

# Structural response: displacement at site i due to mutation at site j
head(responses$Dr2i)
#> # A tibble: 6 × 3
#>       i     j     Dr2i
#>   <int> <int>    <dbl>
#> 1     1     1 0.0304  
#> 2     2     1 0.0119  
#> 3     3     1 0.0117  
#> 4     4     1 0.00305 
#> 5     5     1 0.00422 
#> 6     6     1 0.000403
```

See `?mutenm` for the full package documentation.
