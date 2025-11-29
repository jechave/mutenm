# mutenm Package Dependencies

This document describes both external package dependencies and internal code dependencies (what calls what).

---

## 1. External Dependencies

### 1.1 From DESCRIPTION (Imports)

| Package | Purpose | Used In |
|---------|---------|---------|
| `bio3d` | PDB file reading, atom selection | `enm_utils_nodes.R` |
| `pracma` | Cross product (`pracma::cross`) | `enm_utils_nodes.R` (qb.levitt) |
| `dplyr` | Data manipulation (filter, mutate, group_by, etc.) | Throughout |
| `tibble` | Tibble data structures (tibble, as_tibble, lst) | Throughout |
| `tidyr` | Data reshaping (pivot_longer, expand_grid, unnest) | `mrs.R`, `enm.R` |
| `magrittr` | Pipe operator (`%>%`) | Throughout |
| `stats` | Random number generation (rnorm) | `mutenm.R` |
| `jefuns` | `beta_boltzmann()` for thermodynamic calculations | `delta_energy.R`, `enm_energy_activation.R` |

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
│  enm()                  mrs()                                        │
│  (enm.R)               (mrs.R)                                       │
├─────────────────────────────────────────────────────────────────────┤
│                    INTERMEDIATE FUNCTIONS                            │
├─────────────────────────────────────────────────────────────────────┤
│  mutenm()                                                            │
│  (mutenm.R)                                                          │
├─────────────────────────────────────────────────────────────────────┤
│              COMPARISON FUNCTIONS (wt vs mut)                        │
├─────────────────────────────────────────────────────────────────────┤
│  delta_structure_dr2i()    delta_motion_dmsfi()    Dv_min()         │
│  delta_structure_dr2n()    delta_motion_dmsfn()    Dg_ent()          │
│                                                    delta_energy_*()  │
├─────────────────────────────────────────────────────────────────────┤
│                       GETTERS & ANALYSIS                             │
├─────────────────────────────────────────────────────────────────────┤
│  get_*() (enm_getters.R)    get_*() (enm_analysis.R)                │
│  utils.R                                                             │
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
enm(pdb, node, model, d_max)
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
mutenm(wt, site_mut, mutation, mut_model, mut_dl_sigma, mut_sd_min, seed)
    │
    ├── [if mut_model == "lfenm"]
    │       └── mutenm_lfenm()
    │               ├── generate_delta_lij(wt, site_mut, mut_sd_min, mut_dl_sigma)
    │               │       └── uses get_graph(wt), rnorm()
    │               ├── calculate_force(wt, delta_lij)
    │               │       └── uses get_graph(wt), get_eij(wt)
    │               └── calculate_dxyz(wt, f)
    │                       └── uses get_cmat(wt)  → δr = C × f
    │
    └── [if mut_model == "sclfenm"]
            └── mutenm_sclfenm()
                    ├── (same as lfenm for initial displacement)
                    └── mutate_enm(mut)            → recalculates ENM
                            ├── mutate_graph()
                            ├── set_enm_eij()
                            ├── set_enm_kmat()
                            └── set_enm_nma()
```

### 3.3 Mutation Response Scanning (mrs)

```
mrs(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, responses, seed)
    │
    ├── Validate responses
    │
    ├── For each site j:
    │       └── For each mutation m:
    │               ├── mutenm(wt, j, m, ...)
    │               │
    │               ├── Calculate requested responses:
    │               │       ├── delta_structure_dr2i(wt, mut)
    │               │       ├── delta_structure_dr2n(wt, mut)
    │               │       ├── delta_motion_dmsfi(wt, mut)
    │               │       ├── delta_motion_dmsfn(wt, mut)
    │               │       ├── Dv_min(wt, mut)
    │               │       ├── Dg_ent(wt, mut)
    │               │       ├── delta_energy_dvs(wt, mut)
    │               │       ├── delta_energy_act_dv(wt, mut)
    │               │       └── delta_energy_act_tds(wt, mut)
    │               │
    │               └── Discard mutant (rm(mut))
    │
    └── Return tibbles with averaged responses
```

### 3.4 Energy Calculations

```
Dv_min(wt, mut)               → v_min(mut) - v_min(wt)
Dg_ent(wt, mut)               → g_ent(mut) - g_ent(wt)
delta_energy_act_dv(wt, mut)  → dgact_dv(mut) - dgact_dv(wt)
delta_energy_act_tds(wt, mut) → dgact_tds(mut) - dgact_tds(wt)
```

---

## 4. Getter Functions

### 4.1 Exported Getters (enm_getters.R)

| Getter | Returns | Used By |
|--------|---------|---------|
| `get_enm_param(prot)` | `prot$param` | mrs |
| `get_nsites(prot)` | `prot$nodes$nsites` | Throughout |
| `get_site(prot)` | `prot$nodes$site` | mrs |
| `get_pdb_site(prot)` | `prot$nodes$pdb_site` | active_site_indexes |
| `get_bfactor(prot)` | `prot$nodes$bfactor` | - |
| `get_xyz(prot)` | `prot$nodes$xyz` | Multiple |

### 4.2 Internal Getters (enm_getters.R)

| Getter | Returns | Used By |
|--------|---------|---------|
| `get_enm_node(prot)` | `prot$param$node` | set_enm_nodes |
| `get_enm_model(prot)` | `prot$param$model` | set_enm_graph |
| `get_d_max(prot)` | `prot$param$d_max` | set_enm_graph, get_cn |
| `get_graph(prot)` | `prot$graph` | mutenm.R, enm_analysis.R |
| `get_eij(prot)` | `prot$eij` | mutenm.R |
| `get_kmat(prot)` | `prot$kmat` | delta functions |
| `get_mode(prot)` | `prot$nma$mode` | - |
| `get_evalue(prot)` | `prot$nma$evalue` | Multiple |
| `get_umat(prot)` | `prot$nma$umat` | delta functions |
| `get_cmat(prot)` | `prot$nma$cmat` | mutenm.R |
| `get_nmodes(prot)` | `max(prot$nma$mode)` | mrs |

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

## 6. Utility Functions (utils.R, internal)

| Function | Used By |
|----------|---------|
| `my_as_xyz(r)` | Throughout (converts vector to 3×N matrix) |
| `wcn_xyz(xyz)` | `get_wcn()` |
| `cn_xyz(xyz, d_max)` | `get_cn()` |
| `reduce_matrix(m)` | `get_reduced_cmat()`, `get_reduced_kmat()` |
| `dactive.xyz(xyz, site_active)` | `get_dactive()` |
| `xyz_indices_site(site)` | `active_site_indexes()` |
| `my_quad_form(x, m, y)` | `dgact_dv()` |

---

## 7. Energy Functions

### 7.1 enm_energy.R

| Function | Depends On |
|----------|------------|
| `v_min(prot)` | `get_graph()`, `v_dij()` |
| `g_ent(prot, beta)` | `get_evalue()`, `enm_g_entropy_mode()` |
| `v_dij(dij, v0ij, kij, lij)` | (internal) |
| `enm_g_entropy_mode(energy, beta)` | (internal) |

### 7.2 enm_energy_activation.R

| Function | Depends On |
|----------|------------|
| `dgact_dv(prot, ideal, pdb_site_active)` | `kmat_asite()`, `dxyz_asite()`, `my_quad_form()` |
| `dgact_tds(prot, ideal, pdb_site_active, beta)` | `kmat_asite()`, `enm_g_entropy_mode()` |
| `active_site_indexes(prot, pdb_site_active)` | `get_site()`, `get_pdb_site()`, `xyz_indices_site()` |
| `kmat_asite(prot, pdb_site_active)` | `active_site_indexes()`, `get_cmat()`, `solve()` |
| `dxyz_asite(prot, ideal, pdb_site_active)` | `active_site_indexes()`, `get_xyz()` |

### 7.3 delta_energy.R

| Function | Depends On |
|----------|------------|
| `Dv_min(wt, mut)` | `v_min()` |
| `Dg_ent(wt, mut, beta)` | `g_ent()` |
| `delta_energy_dvs(wt, mut, ideal)` | `calculate_vs()` |
| `delta_energy_act_dv(wt, mut, ideal, pdb_site_active)` | `dgact_dv()` |
| `delta_energy_act_tds(wt, mut, ideal, pdb_site_active, beta)` | `dgact_tds()` |

---

## 8. File Summary

| File | Purpose | Exports |
|------|---------|---------|
| `enm.R` | ENM construction | `enm` |
| `enm_getters.R` | Access prot components | `get_enm_param`, `get_nsites`, `get_site`, `get_pdb_site`, `get_bfactor`, `get_xyz` |
| `enm_analysis.R` | Derived properties | `get_cn`, `get_wcn`, `get_dactive`, `get_msf_site`, `get_msf_mode`, `get_mlms`, `get_stress`, `get_rho_matrix`, `get_reduced_cmat`, `get_reduced_kmat`, `get_msf_site_mode`, `get_umat2` |
| `enm_energy.R` | ENM energies | `v_min`, `g_ent` |
| `enm_energy_activation.R` | Activation energies | `dgact_dv`, `dgact_tds` |
| `enm_utils_nodes.R` | Node coordinate calculation | (internal) |
| `enm_utils_kij_functions.R` | Spring constant models | (internal) |
| `mutenm.R` | Mutation perturbation | `mutenm` |
| `mrs.R` | Mutation response scanning | `mrs` |
| `delta_structure_by_site.R` | Structure response (per site) | `delta_structure_dr2i` |
| `delta_structure_by_mode.R` | Structure response (per mode) | `delta_structure_dr2n` |
| `delta_motion_by_site.R` | Dynamics response (per site) | `delta_motion_dmsfi` |
| `delta_motion_by_mode.R` | Dynamics response (per mode) | `delta_motion_dmsfn` |
| `delta_energy.R` | Energy differences | `Dv_min`, `Dg_ent`, `delta_energy_dvs`, `delta_energy_act_dv`, `delta_energy_act_tds` |
| `utils.R` | General utilities | (internal) |
| `mutenm-imports.R` | Package imports | - |
| `mutenm-package.R` | Package documentation | - |

---

*Generated for mutenm package development reference*
