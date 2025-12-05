---
title: 'mutenm: An R package for simulating mutation effects on protein structure using elastic network models'
tags:
  - R
  - proteins
  - elastic network models
  - mutations
  - protein evolution
  - structural biology
authors:
  - name: Julian Echave
    orcid: 0000-0003-0678-5593
    affiliation: 1
affiliations:
  - name: Instituto de Ciencias Físicas (ICIFI-CONICET), Universidad Nacional de San Martín, Martín de Irigoyen 3100, 1650 San Martín, Buenos Aires, Argentina
    index: 1
author:
  - Julian Echave^1^
date: December 2025
bibliography: paper.bib
---

^1^ Instituto de Ciencias Físicas (ICIFI-CONICET), Universidad Nacional de San Martín, Martín de Irigoyen 3100, 1650 San Martín, Buenos Aires, Argentina

# Summary

`mutenm` is an R package for building elastic network models (ENMs) of proteins and simulating the effects of mutations using linear response theory. An ENM represents a protein as a network of nodes (typically one per residue) connected by harmonic springs. Despite their simplicity, ENMs capture the functionally relevant low-frequency collective motions of proteins through normal mode analysis.

Given a protein structure, `mutenm` constructs an ENM and calculates how mutations at each site propagate structural effects throughout the protein. The package provides Mutation Response Scanning (MRS), which systematically simulates mutations at every site and quantifies the resulting structural displacements, producing response matrices and profiles that reveal structural communication pathways within the protein.

# Statement of Need

Understanding how mutations affect protein structure and function is central to protein evolution, disease genetics, and protein engineering. Experimental approaches like deep mutational scanning generate large datasets, but computational methods are needed to interpret these results and predict effects of unstudied mutations.

Elastic network models provide a computationally efficient way to capture protein flexibility and have been widely used to study protein dynamics. However, most ENM software focuses on analyzing wild-type proteins rather than predicting mutation effects. The Perturbation Response Scanning (PRS) method in ProDy [@bakan2011; @zhang2021] is a related approach, but it applies external forces to residues rather than perturbing the ENM Hamiltonian itself.

The `mutenm` package implements the Linearly Forced ENM (LFENM) approach [@echave2008; @echave2010], which models mutations as perturbations to spring equilibrium lengths—a more physically realistic representation of how amino acid substitutions alter local interactions. These methods have been applied to study how mutation and selection shape protein structural divergence [@marcos2020; @echave2025].

# Functionality

The package implements a complete workflow for mutation analysis through four functions:

**Building the model.** The `enm()` function constructs an elastic network model from a PDB structure. It supports multiple node representations (Cα, Cβ, side-chain centroid) and several ENM variants (ANM, Ming-Wall, parameter-free ANM, and others). The resulting model contains the network topology, spring constants, and normal mode analysis results including the covariance matrix used for linear response calculations.

**Simulating single mutations.** The `mutenm()` function simulates a mutation at a specified site by perturbing spring equilibrium lengths around the mutated residue. The structural response is calculated using linear response theory: δr = C × f, where C is the covariance matrix and f is the force arising from the perturbed springs. This provides the displacement at every site caused by the mutation.

**Mutation response scanning.** The `mrs()` function systematically applies `mutenm()` to every site, averaging over multiple random perturbations per site to obtain robust estimates. The output includes a response matrix dr2ij[i,j] representing the mean squared displacement at site i caused by mutations at site j (Figure 1, left panel). From this matrix, two informative profiles are derived: the influence profile (column means—how much mutations at each site affect the protein globally) and the sensitivity profile (row means—how much each site responds to mutations anywhere).

**Visualization.** The `plot_mrs()` function creates publication-ready figures showing the response matrix as a heatmap alongside the influence and sensitivity profiles (Figure 1). Sites with high influence are structurally important—mutations there cause large displacements throughout the protein. Sites with high sensitivity are structurally responsive—they move substantially regardless of where the mutation occurs. These profiles correlate with evolutionary conservation, as sites with high influence or sensitivity tend to evolve more slowly [@echave2025].

![Mutation response analysis of acylphosphatase (PDB: 2ACY) using `mrs()` and `plot_mrs()`. Left: Response matrix showing mean squared displacement at each site (y-axis) caused by mutations at each site (x-axis), with log10 color scale. Right: Influence profile (top) showing the effect of mutating each site, and sensitivity profile (bottom) showing how each site responds to mutations elsewhere.](paper_figure.png){width=100%}

# Acknowledgements

This work was supported by CONICET (Argentina).

# References
