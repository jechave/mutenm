mutenm package
================
Julian Echave

<!-- README.md is generated from README.Rmd. Please edit that file -->

The `mutenm` package contains functions to build Elastic Network Models
(ENM) of proteins and to perturb them; `mutenm` stands for Mutating
Elastic Network Models.

# Installation

Install packages `mutenm` (this package) and `jefuns` (miscellaneous
functions, some of which `mutenm` uses).

    # install.packages("devtools")

    devtools::install_github("jechave/jefuns")
    devtools::install_github("jechave/mutenm")

# Quick Start

``` r
library(mutenm)
```

## Set up the ENM for a protein

First, read a pdb file using `bio3d::read.pdb` to generate a `pdb`
object for a protein. Then, create the `prot` object, that contains the
full ENM analysis.

``` r
pdb <- bio3d::read.pdb("data-raw/2XWRa.pdb") # read a pdb file
wt <- set_enm(pdb, node = "calpha", model = "anm", d_max = 10.5)
```

`wt` created here by `set_enm()` is an object of class *prot*. In this
example, network nodes are placed at $C_\alpha$ coordinates, the model
used is Bahar’s Anisotropic Network Model (`model = "anm"`) with a
cut-off distance to define contacts of `d_max = 10.5`.

## Create a mutant

Use `get_mutant_site()` to create a mutant at a specific site:

``` r
mut <- get_mutant_site(wt, site = 50, mutation = 1, mut_model = "lfenm", mut_dl_sigma = 0.3, seed = 123)
```

This creates a mutant protein by perturbing contacts at site 50 using
the LFENM model.

## Compare wild-type and mutant

Calculate the stability change (ΔΔG) between wild-type and mutant:

``` r
ddg <- ddg_dv(wt, mut)
ddg
#> [1] 0.8618287
```

See the package documentation for more detailed examples and mutation
scanning functions.
