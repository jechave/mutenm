---
title: '`mutenm`: An R package to mutate elastic network models of proteins'
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
bibliography: paper.bib
---

\begin{center}
\textbf{Julian Echave}

Instituto de Ciencias Físicas (ICIFI-CONICET), Universidad Nacional de San Martín, Martín de Irigoyen 3100, 1650 San Martín, Buenos Aires, Argentina

jechave@unsam.edu.ar
\end{center}

# Summary

`mutenm` is an R package for modeling mutations within the elastic network model (ENM) framework. An ENM represents a protein as a network of nodes connected by harmonic springs; despite their simplicity, ENMs capture functionally relevant collective motions through normal mode analysis.

The core function `mutenm()` takes an ENM of a wild-type protein and generates an ENM of the mutant. This enables studying how mutations affect protein structure and dynamics: comparing wild-type and mutant properties, performing mutation response scans, or simulating evolutionary trajectories through iterated mutation.

# Statement of Need

Elastic network models are widely used to study protein dynamics. Packages such as Bio3D [@grant2006; @skjaerven2014] in R and ProDy [@bakan2011; @zhang2021] in Python provide comprehensive ENM functionality. However, these tools focus on analyzing wild-type proteins.

The `mutenm` package extends ENM analysis to mutations. It implements the Linearly Forced ENM (LFENM) [@echave2008; @echave2010], which models a mutation as a perturbation to spring equilibrium lengths around the mutation site. Given a wild-type ENM, the function `mutenm()` generates a mutant ENM with altered structure. A self-consistent variant (scLFENM) [@echave2012] additionally recalculates the mutant's normal modes, capturing changes in dynamics as well as structure.

This approach has been used to study protein evolution—the divergence of structure across sites [@marcos2015; @marcos2020; @echave2024; @echave2025], across modes [@echave2008; @echave2010], and the evolution of dynamics [@echave2012]. The package is intended for computational biologists interested in protein structure, dynamics, and evolution.

# Functionality

**Building the ENM.** The `enm()` function constructs an elastic network model from a PDB structure. Multiple node representations are supported (Cα, Cβ, side-chain centroid); the optimal choice depends on the application [@marcos2015; @echave2024]. Several force fields are available [@atilgan2001; @yang2009; @ming2005; @hinsen1998; @hinsen2000]; results are generally robust across models.

**Mutating the ENM.** The `mutenm()` function takes a wild-type ENM and generates a mutant ENM by perturbing spring equilibrium lengths around the mutation site. Two mutation models are available: LFENM perturbs the structure while keeping the Hessian fixed; scLFENM recalculates the Hessian after relaxation, capturing changes in dynamics. The mutant ENM can be analyzed like any ENM—comparing normal modes, covariance matrices, or flexibility between wild-type and mutant. By iterating `mutenm()`, users can simulate evolutionary trajectories.

**Mutation response scanning.** The `mrs()` function applies `mutenm()` systematically across all sites [@echave2021], producing a response matrix and profiles analogous to Perturbation Response Scanning (PRS) [@bakan2011] but based on a mutation model rather than arbitrary forces. The response matrix (Figure 1) reveals how mutations at each site affect structure throughout the protein. The derived influence and sensitivity profiles identify structurally important and structurally responsive sites, which correlate with evolutionary conservation [@echave2025].

**Visualization.** The `plot_mrs()` function creates publication-ready figures showing the response matrix alongside influence and sensitivity profiles (Figure 1).

# Availability

`mutenm` is available at https://github.com/jechave/mutenm with documentation and vignettes.

# Acknowledgements

This work was supported by CONICET (grant PIP-11220210100462).

# References

::: {#refs}
:::

\newpage

![Mutation response analysis of acylphosphatase (PDB: 2ACY) using `mrs()` and `plot_mrs()`. Left: Response matrix showing mean squared displacement at each site (y-axis) caused by mutations at each site (x-axis), with log10 color scale. Right: Influence profile (top) showing the effect of mutating each site, and sensitivity profile (bottom) showing how each site responds to mutations elsewhere.](paper_figure.png){width=100%}
