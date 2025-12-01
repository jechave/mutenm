# mutenm Package Overview

The `mutenm` package builds Elastic Network Models (ENMs) from protein structures and simulates the effects of mutations using linear response theory. It provides tools to analyze how mutations at specific sites affect protein structure, dynamics, and energetics.

---

## Workflow

A typical analysis has three steps:

1. **Build the ENM** from a PDB structure using `enm()`
2. **Simulate mutations** using either:
   - `mutenm()` for a single mutation at a specific site
   - `mrs()` for scanning all sites at once (mutation response scanning)
3. **Analyze the effects** using property functions (structure, dynamics, energy)

```
PDB file
    ↓
enm(pdb, node, model, d_max) → prot object (wild-type ENM)
    ↓
    ├── mutenm(wt, site_mut, ...) → prot object (mutant)
    │       ↓
    │       Compare wt and mut using D-functions
    │
    └── mrs(wt, nmut, ...) → response tibbles
```

The `prot` object contains the protein's graph (nodes and edges), the Hessian matrix, and the results of normal mode analysis (eigenvalues, eigenvectors, covariance matrix).

---

## Functions

### Core (2 functions)

| Function | Description |
|----------|-------------|
| `enm(pdb, node, model, d_max)` | Build an ENM from a PDB structure. Returns a `prot` object. |
| `mutenm(wt, site_mut, ...)` | Create a mutant by perturbing edges at a specific site. Returns a `prot` object. |

These are the core functions of the package: `enm()` builds the model, `mutenm()` simulates a mutation.

#### Mutation Models

The `mutenm()` function supports two mutation models:

| Model | Description |
|-------|-------------|
| `lfenm` | Linear Force ENM — uses linear response to estimate structural change. Only coordinates are updated. |
| `sclfenm` | Self-Consistent LFENM — recalculates the full ENM after perturbation (graph, Hessian, normal modes, covariance matrix). |

Use `lfenm` for structural effects. Use `sclfenm` when you need dynamics or entropy changes.

### Scanning (1 function)

| Function | Description |
|----------|-------------|
| `mrs(wt, nmut, ...)` | Mutation Response Scanning — loops over all sites, averages over nmut mutations per site, returns response tibbles. |

### Single Protein (9 functions)

Functions that take one `prot` object and return a property.

#### Structure

| Function | Description |
|----------|-------------|
| `cn(prot)` | Contact number (neighbors within cutoff) |
| `wcn(prot)` | Weighted contact number (Σ 1/d²) |
| `dactive(prot, pdb_site_active)` | Distance to active site |

#### Dynamics

| Function | Description |
|----------|-------------|
| `msfi(prot)` | Mean-square fluctuation at each site |
| `msfn(prot)` | MSF contribution from each mode |

#### Energy

| Function | Description |
|----------|-------------|
| `v_min(prot)` | Minimum (stress) energy |
| `g_ent(prot, beta)` | Entropic free energy |
| `dv_act(prot, ideal, pdb_site_active)` | Activation energy (internal) |
| `dg_ent_act(prot, ideal, pdb_site_active, beta)` | Activation entropy |

### Pair Comparison (8 functions)

Functions that compare wild-type and mutant: `Df(wt, mut)`.

#### Structure

| Function | Description |
|----------|-------------|
| `Dr2i(wt, mut)` | Squared displacement at each site |
| `Dr2n(wt, mut)` | Squared displacement projected onto each mode |

#### Dynamics

| Function | Description |
|----------|-------------|
| `Dmsfi(wt, mut)` | Change in MSF at each site |
| `Dmsfn(wt, mut)` | Change in MSF per mode |

Note: `Dmsfi` and `Dmsfn` require the `sclfenm` mutation model (full ENM recalculation) to be meaningful.

#### Energy

| Function | Description |
|----------|-------------|
| `Dv_min(wt, mut)` | Change in minimum energy |
| `Dg_ent(wt, mut, beta)` | Change in entropic free energy |
| `Ddv_act(wt, mut, ...)` | Change in activation energy |
| `Ddg_ent_act(wt, mut, ...)` | Change in activation entropy |

---

## Key Equations

### Symbols

| Symbol | Description |
|--------|-------------|
| `C` | Covariance matrix (3N × 3N) — encodes thermal fluctuations |
| `k_ij` | Spring constant for edge between sites i and j |
| `d_ij` | Current distance between sites i and j |
| `l_ij` | Equilibrium (rest) length of edge i-j |
| `ê_ij` | Unit vector along edge i-j |
| `δl_ij` | Perturbation to edge length (mutation effect) |
| `v0_ij` | Baseline energy of edge i-j |

### Linear Response

The structural response to a mutation is computed using linear response theory:

```
δr = C × f
```

where `C` is the covariance matrix and `f` is the force vector from edge perturbations.

### Force from Perturbation

When an edge length is perturbed by `δl_ij`, it generates a force:

```
f_ij = -k_ij × δl_ij × ê_ij
```

### Minimum Energy

The total energy of the ENM at the current conformation:

```
V_min = Σ v0_ij + ½ Σ k_ij × (d_ij - l_ij)²
```

### Stress Energy

The energy cost to deform the active site to an ideal conformation:

```
V_s = ½ Σ k_ij × (d_ij^ideal - l_ij)²
```

---

## See Also

- Run `?enm`, `?mutenm`, `?mrs` in R for detailed function documentation
- `function_grouping.md` — rationale for function organization
