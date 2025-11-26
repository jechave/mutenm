# mutenm Package Dependencies

This document describes both external package dependencies and internal code dependencies (what calls what).

---

## 1. External Dependencies

### 1.1 From DESCRIPTION (Imports)

| Package | Purpose | Used In |
|---------|---------|---------|
| `bio3d` | PDB file reading, atom selection | `enm_utils_nodes.R` |
| `Matrix` | Sparse matrix operations for efficiency | `mutscan_smrs.R` |
| `pracma` | Cross product (`pracma::cross`) | `enm_utils_nodes.R` (qb.levitt) |
| `dplyr` | Data manipulation (filter, mutate, group_by, etc.) | Throughout |
| `tibble` | Tibble data structures (tibble, as_tibble, lst) | Throughout |
| `tidyr` | Data reshaping (unnest, expand_grid) | Mutation scanning |
| `purrr` | Functional programming (map, map2, pmap) | Mutation scanning |
| `magrittr` | Pipe operator (`%>%`) | Throughout |
| `stats` | Random number generation (rnorm) | `penm.R` |
| `matrixStats` | (Listed but not actively used in code) | - |
| `jefuns` | `matrix_to_tibble`, `plot_matrix` | `mutscan_smrs.R` |

### 1.2 From bio3d (imported functions)

| Function | Used For |
|----------|----------|
| `atom.select()` | Select atoms by type (calpha, backbone, CB, etc.) |
| `combine.select()` | Combine atom selections |
| `com()` | Calculate center of mass |
| `aa321()` | Convert 3-letter to 1-letter amino acid codes |

---

## 2. Internal Code Structure

### 2.1 File Organization by Layer

```
┌─────────────────────────────────────────────────────────────────────┐
│                         USER-FACING API                              │
├─────────────────────────────────────────────────────────────────────┤
│  set_enm()          smrs()         mrs_all()      smrs_ddg()         │
│  (enm.R)        (mutscan_smrs.R) (mrs_structure.R) (mutscan_smrs_ddg.R)│
├─────────────────────────────────────────────────────────────────────┤
│                    INTERMEDIATE FUNCTIONS                            │
├─────────────────────────────────────────────────────────────────────┤
│  get_mutant_site()     generate_mutants()    mrs_all()              │
│  (penm.R)           (mrs_generate_mutants.R) (mrs_structure.R)      │
├─────────────────────────────────────────────────────────────────────┤
│              COMPARISON FUNCTIONS (wt vs mut)                        │
├─────────────────────────────────────────────────────────────────────┤
│  delta_structure_*()    delta_motion_*()      delta_energy_*()      │
│  (delta_structure_*.R)  (delta_motion_*.R)    (delta_energy.R)      │
├─────────────────────────────────────────────────────────────────────┤
│                    RESPONSE MATRIX BUILDERS                          │
├─────────────────────────────────────────────────────────────────────┤
│  mrs_structure_*()      mrs_motion_*()                              │
│  (mrs_structure_*.R)    (mrs_motion_*.R)                            │
├─────────────────────────────────────────────────────────────────────┤
│                       GETTERS & UTILITIES                            │
├─────────────────────────────────────────────────────────────────────┤
│  get_*() (enm_getters.R)    get_*() (enm_analysis.R)                │
│  utils.R                     mutscan_utils.R                         │
├─────────────────────────────────────────────────────────────────────┤
│                    ENM BUILDING BLOCKS                               │
├─────────────────────────────────────────────────────────────────────┤
│  calculate_enm_*()      kij_*()              prot_*()               │
│  (enm.R)           (enm_utils_kij_functions.R) (enm_utils_nodes.R)  │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 3. Detailed Dependency Graph

### 3.1 ENM Construction Chain

```
set_enm(pdb, node, model, d_max)
    │
    ├── create_enm()                           → creates empty prot structure
    │
    ├── set_enm_param()                        → stores node, model, d_max
    │
    ├── set_enm_nodes(pdb)
    │       └── calculate_enm_nodes(pdb, node)
    │               ├── prot_ca(pdb)           → uses atom.select(pdb, "calpha")
    │               ├── prot_sc(pdb)           → uses residue.coordinates(), residue.bfactors()
    │               │       ├── qb.micheletti()
    │               │       └── qb.levitt()    → uses pracma::cross()
    │               └── prot_cb(pdb)           → similar to prot_sc
    │
    ├── set_enm_graph()
    │       └── calculate_enm_graph(xyz, pdb_site, model, d_max)
    │               ├── dij_edge()             → calculates edge distances
    │               ├── sdij_edge()            → calculates sequence distances
    │               └── kij_*(dij, sdij)       → spring constant (model-specific)
    │                       ├── kij_anm()
    │                       ├── kij_ming_wall()
    │                       ├── kij_hnm(), kij_hnm0()
    │                       ├── kij_pfanm()
    │                       └── kij_reach()
    │
    ├── set_enm_eij()
    │       └── calculate_enm_eij(xyz, i, j)   → unit vectors along edges
    │
    ├── set_enm_kmat()
    │       └── calculate_enm_kmat(graph, eij, nsites) → Hessian matrix
    │
    └── set_enm_nma()
            └── calculate_enm_nma(kmat)        → eigen decomposition
                    → returns: mode, evalue, umat, cmat
```

### 3.2 Mutation Chain

```
get_mutant_site(wt, site_mut, mutation, mut_model, mut_dl_sigma, mut_sd_min, seed)
    │
    ├── [if mut_model == "lfenm"]
    │       └── get_mutant_site_lfenm()
    │               ├── generate_delta_lij(wt, site_mut, mut_sd_min, mut_dl_sigma)
    │               │       └── uses get_graph(wt), rnorm()
    │               ├── calculate_force(wt, delta_lij)
    │               │       └── uses get_graph(wt), get_eij(wt)
    │               └── calculate_dxyz(wt, f)
    │                       └── uses get_cmat(wt)  → δr = C × f
    │
    └── [if mut_model == "sclfenm"]
            └── get_mutant_site_sclfenm()
                    ├── (same as lfenm for initial displacement)
                    └── mutate_enm(mut)            → recalculates ENM
                            ├── mutate_graph()
                            ├── set_enm_eij()
                            ├── set_enm_kmat()
                            └── set_enm_nma()
```

### 3.3 Simulation Scanning (smrs)

```
smrs(wt, nmut, mut_dl_sigma, mut_sd_min, option, response, seed)
    │
    ├── generate_perturbations(wt, nmut, mut_dl_sigma, mut_sd_min, seed)
    │       ├── generate_delta_lij()  (for each site, each mutation)
    │       └── calculate_force()
    │       → returns: dlmat[edge, site, mutation], fmat[3N, site, mutation]
    │
    ├── [option == "site"]
    │       ├── calculate_dr2ij_smrs() → calculate_s2ij(fmat, cmat)
    │       ├── calculate_de2ij_smrs() → calculate_s2ij(fmat, cmat_sqrt)
    │       ├── calculate_df2ij_smrs() → calculate_s2ij(fmat, identity)
    │       └── calculate_dvsij_smrs()
    │
    └── [option == "mode"]
            └── calculate_*nj_smrs() → calculate_fnmat() + averaging

calculate_s2ij(fmat, amat)
    │
    └── Uses Matrix() for sparse multiplication
        Returns averaged response matrix
```

### 3.4 Slow Simulation Method (mrs_all)

```
mrs_all(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, seed)
    │
    ├── generate_mutants(wt, nmut, mut_model, ...)
    │       └── calls get_mutant_site() for each (site, mutation) pair
    │       → returns tibble with columns: wt, j, mutation, mut
    │
    └── mrs_structure_*ij(mutants)
            │
            └── For each mutant:
                    ├── delta_structure_dr2i(wt, mut)
                    ├── delta_structure_de2i(wt, mut, kmat_sqrt)
                    ├── delta_structure_df2i(wt, mut)
                    └── delta_structure_dvsi_same_topology(wt, mut)
```

### 3.5 DDG Calculations

```
smrs_ddg(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, seed)
    │
    └── Uses smrs to calculate de2ij and dvsij, then ddg = Σ(dvsij - de2ij) / 2
```

---

## 4. Getter Functions (enm_getters.R)

These are used throughout the package to access prot object components:

| Getter | Returns | Used By |
|--------|---------|---------|
| `get_enm_param(prot)` | `prot$param` | mrs_all, smrs_all |
| `get_enm_node(prot)` | `prot$param$node` | set_enm_nodes |
| `get_enm_model(prot)` | `prot$param$model` | set_enm_graph |
| `get_d_max(prot)` | `prot$param$d_max` | set_enm_graph, get_cn |
| `get_nsites(prot)` | `prot$nodes$nsites` | Throughout |
| `get_site(prot)` | `prot$nodes$site` | generate_mutants, smrs_ddg |
| `get_pdb_site(prot)` | `prot$nodes$pdb_site` | smrs_ddg, active_site_indexes |
| `get_bfactor(prot)` | `prot$nodes$bfactor` | - |
| `get_xyz(prot)` | `prot$nodes$xyz` | Multiple |
| `get_graph(prot)` | `prot$graph` | penm.R, mutscan_*.R |
| `get_eij(prot)` | `prot$eij` | penm.R, mutscan_smrs.R |
| `get_kmat(prot)` | `prot$kmat` | delta_structure_df2i |
| `get_mode(prot)` | `prot$nma$mode` | mutscan_smrs.R |
| `get_evalue(prot)` | `prot$nma$evalue` | Multiple |
| `get_umat(prot)` | `prot$nma$umat` | mutscan_smrs.R, get_cmat_sqrt |
| `get_cmat(prot)` | `prot$nma$cmat` | penm.R, mutscan_smrs.R |
| `get_nmodes(prot)` | `max(prot$nma$mode)` | mutscan_smrs.R |

---

## 5. Analysis Functions (enm_analysis.R)

| Function | Depends On | Returns |
|----------|------------|---------|
| `get_cn(prot)` | `cn_xyz()`, `get_xyz()`, `get_d_max()` | Contact number profile |
| `get_wcn(prot)` | `wcn_xyz()`, `get_xyz()` | Weighted contact number |
| `get_dactive(prot, sites)` | `dactive.xyz()`, `active_site_indexes()` | Distance to active site |
| `get_msf_site(prot)` | `get_reduced_cmat()` | Mean-square fluctuation per site |
| `get_msf_mode(prot)` | `get_evalue()` | MSF per mode (1/λ) |
| `get_mlms(prot, sdij_cut)` | `get_graph()` | Mean local mutational stress |
| `get_stress(prot)` | `get_graph()` | Stress energy per site |
| `get_rho_matrix(prot)` | `get_reduced_cmat()` | Correlation matrix |
| `get_reduced_cmat(prot)` | `get_cmat()`, `reduce_matrix()` | N×N covariance |
| `get_reduced_kmat(prot)` | `get_kmat()`, `reduce_matrix()` | N×N Kirchhoff |
| `get_msf_site_mode(prot)` | `get_umat2()`, `get_evalue()` | MSF[site,mode] matrix |
| `get_umat2(prot)` | `get_umat()` | |u_i,n|² matrix |

---

## 6. Utility Functions (utils.R)

| Function | Used By |
|----------|---------|
| `my_as_xyz(r)` | Throughout (converts vector to 3×N matrix) |
| `wcn_xyz(xyz)` | `get_wcn()` |
| `cn_xyz(xyz, d_max)` | `get_cn()` |
| `cn_graph(graph)` | - |
| `reduce_matrix(m)` | `get_reduced_cmat()`, `get_reduced_kmat()` |
| `dactive.xyz(xyz, site_active)` | `get_dactive()` |
| `xyz_indices_site(site)` | `active_site_indexes()` |
| `my_quad_form(x, m, y)` | `dgact_dv()` |
| `tr(m)` | `dbhat()`, `rwsip()`, `dh()` |
| `logdet(m)` | `dbhat()` |
| `dbhat(ca, cb)` | `delta_motion_dbhati()` |
| `rwsip(ca, cb)` | `delta_motion_rwsipi()` |
| `dh(ca, cb)` | `delta_motion_dhi()` |

---

## 7. Energy Functions

### enm_energy.R

| Function | Depends On |
|----------|------------|
| `enm_v_min(prot)` | `get_graph()`, `v_dij()` |
| `enm_g_entropy(prot, beta)` | `get_evalue()`, `enm_g_entropy_mode()` |
| `v_dij(dij, v0ij, kij, lij)` | (internal) |
| `enm_g_entropy_mode(energy, beta)` | (internal) |

### enm_energy_activation.R

| Function | Depends On |
|----------|------------|
| `dgact_dv(prot, ideal, pdb_site_active)` | `kmat_asite()`, `dxyz_asite()`, `my_quad_form()` |
| `dgact_tds(prot, ideal, pdb_site_active, beta)` | `kmat_asite()`, `enm_g_entropy_mode()` |
| `active_site_indexes(prot, pdb_site_active)` | `get_site()`, `get_pdb_site()`, `xyz_indices_site()` |
| `kmat_asite(prot, pdb_site_active)` | `active_site_indexes()`, `get_cmat()`, `solve()` |
| `dxyz_asite(prot, ideal, pdb_site_active)` | `active_site_indexes()`, `get_xyz()` |

### delta_energy.R

| Function | Depends On |
|----------|------------|
| `ddg_dv(wt, mut)` | `enm_v_min()` |
| `ddg_tds(wt, mut, beta)` | `enm_g_entropy()` |
| `delta_energy_dvs(wt, mut, ideal)` | `calculate_vs()` |
| `ddgact_dv(wt, mut, ideal, pdb_site_active)` | `dgact_dv()` |
| `ddgact_tds(wt, mut, ideal, pdb_site_active, beta)` | `dgact_tds()` |

---

## 8. File Summary

| File | Purpose | Exports |
|------|---------|---------|
| `enm.R` | ENM construction | `set_enm` |
| `enm_getters.R` | Access prot components | `get_*` (multiple) |
| `enm_analysis.R` | Derived properties | `get_cn`, `get_wcn`, `get_msf_*`, etc. |
| `enm_energy.R` | ENM energies | `enm_v_min`, `enm_g_entropy` |
| `enm_energy_activation.R` | Activation energies | `dgact_dv`, `dgact_tds` |
| `enm_utils_nodes.R` | Node coordinate calculation | (internal) |
| `enm_utils_kij_functions.R` | Spring constant models | (internal) |
| `penm.R` | Mutation perturbation | `get_mutant_site` |
| `mutscan_smrs.R` | Simulation scanning | `smrs`, `smrs_all` |
| `mutscan_smrs_ddg.R` | Simulation ΔΔG | `smrs_ddg` |
| `mutscan_smrs_ddgact.R` | Simulation ΔΔG‡ | `smrs_ddgact` |
| `mutscan_utils.R` | Matrix sqrt helpers | (internal) |
| `mrs_generate_mutants.R` | Generate mutant tibble | `generate_mutants` |
| `mrs_structure.R` | Wrapper for mrs_all | `mrs_all` |
| `mrs_structure_by_site.R` | Structure response matrices | `mrs_structure_*ij` |
| `mrs_structure_by_mode.R` | Mode-space structure | `mrs_structure_*nj` |
| `mrs_motion_by_site.R` | Dynamics response matrices | `mrs_motion_*ij` |
| `mrs_motion_by_mode.R` | Mode-space dynamics | `mrs_motion_*nj` |
| `delta_structure_by_site.R` | wt-mut structure comparison | `delta_structure_*i` |
| `delta_structure_by_mode.R` | Mode-space comparison | `delta_structure_*n` |
| `delta_motion_by_site.R` | wt-mut dynamics comparison | `delta_motion_*i` |
| `delta_motion_by_mode.R` | Mode-space dynamics | `delta_motion_*n` |
| `delta_energy.R` | Energy differences | `ddg_*`, `delta_energy_*` |
| `utils.R` | General utilities | (internal) |
| `penm-imports.R` | Package imports | - |
| `mutenm-package.R` | Package documentation | - |
| `defunct_functions.R` | Deprecated functions | - |

---

*Generated for mutenm package development reference*
