---
title: 'mutenm: An R package for simulating mutation effects on protein structure and dynamics using elastic network models'
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

Given a protein structure, `mutenm` constructs an ENM and provides tools for analyzing both the wild-type protein and the effects of mutations. For wild-type analysis, it calculates structural descriptors such as local packing density, dynamic properties including thermal fluctuations per residue and per normal mode, and thermodynamic quantities such as conformational entropy and activation energies. For mutation analysis, the package simulates how amino acid substitutions perturb the elastic network and quantifies the resulting structural displacements, flexibility changes, and energetic costs.

# Statement of Need

Understanding how mutations affect protein structure and function is central to protein evolution, disease genetics, and protein engineering. Experimental approaches like deep mutational scanning generate large datasets, but computational methods are needed to interpret these results and predict effects of unstudied mutations.

Elastic network models provide a computationally efficient way to capture protein flexibility and have been widely used to study protein dynamics. However, most ENM software focuses on analyzing wild-type proteins rather than predicting mutation effects. The Perturbation Response Scanning (PRS) method in ProDy [@bakan2011; @zhang2021] is a related approach, but it applies external forces to residues rather than perturbing the ENM Hamiltonian itself. Importantly, PRS only predicts structural displacements—it cannot capture how mutations affect protein dynamics.

The `mutenm` package implements the Linearly Forced ENM (LFENM) approach [@echave2008; @echave2010], which models mutations as perturbations to spring equilibrium lengths—a more physically realistic representation of how amino acid substitutions alter local interactions. The package also implements the self-consistent LFENM (scLFENM) [@echave2012], which recalculates normal modes after structural relaxation, enabling prediction of mutation effects on protein dynamics and entropy—quantities inaccessible to PRS. These methods have been applied to study how mutation and selection shape protein structural divergence [@marcos2020; @echave2025].

# Functionality

The package implements a complete workflow for mutation analysis through three core functions:

**Building the model.** The `enm()` function constructs an elastic network model from a PDB structure. The resulting model can be analyzed to obtain structural properties (contact number, weighted contact number, distance to active site), dynamic properties (mean-square fluctuations per residue and per mode), and energetic quantities (stress energy, conformational entropy, activation energies). Figure 1A shows the wild-type fluctuation profile for acylphosphatase.

**Simulating single mutations.** The `mutenm()` function simulates a mutation at a specified site by perturbing spring equilibrium lengths. Two models are available: LFENM (fast, structure only) and scLFENM (slower, recalculates dynamics). Functions quantify changes between wild-type and mutant in structure (`Dr2i`, `Dr2n`), dynamics (`Dmsfi`, `Dmsfn`), and energy (`Dv_min`, `Dg_ent`, `Ddv_act`, `Ddg_ent_act`). Figure 1B shows the change in fluctuations caused by a mutation at site 16.

**Mutation response scanning.** The `mrs()` function systematically applies `mutenm()` to every site, averaging over multiple random perturbations per site to obtain robust estimates. The output includes response matrices showing how mutations at each site affect structure and dynamics at every other site (Figure 1C). From these matrices, two informative profiles can be derived: the sensitivity profile (how much each site's dynamics change when mutations occur elsewhere) and the influence profile (how much mutations at each site affect dynamics globally). Sites with high sensitivity or influence tend to evolve more slowly, as mutations there have larger functional consequences (Figure 1D).

![Mutation effects on protein dynamics for acylphosphatase (PDB: 2ACY). A) Wild-type fluctuation profile from `enm()`. B) Change in fluctuations from a mutation at site 16 (dashed line) using `mutenm()`. C) Dynamics response matrix from `mrs()` showing flexibility change at each site (y-axis) caused by mutations at each site (x-axis). D) Sensitivity and influence profiles derived from the response matrix.](paper_figure.png){width=100%}

# Acknowledgements

This work was supported by CONICET (Argentina).

# References
