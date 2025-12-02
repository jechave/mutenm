# mutenm Package Summary

## Overview

The `mutenm` package builds Elastic Network Models (ENMs) from protein structures and simulates mutation effects using linear response theory. It's a development version being trimmed from a larger `penm` package.

## Package Structure

### Files in R/

| File | Purpose | Exports |
|------|---------|---------|
| `enm.R` | ENM construction pipeline | `enm()` |
| `mutenm.R` | Mutation simulation | `mutenm()` |
| `mrs.R` | Mutation response scanning | `mrs()` |
| `single_protein.R` | Properties of one protein | `cn`, `wcn`, `dactive`, `msfi`, `msfn`, `v_min`, `g_ent`, `dv_act`, `dg_ent_act` |
| `pair_comparison.R` | Compare wt vs mut | `Dr2i`, `Dr2n`, `Dmsfi`, `Dmsfn`, `Dv_min`, `Dg_ent`, `Ddv_act`, `Ddg_ent_act` |
| `internals_getters.R` | Access prot components | (internal) |
| `internals_utils.R` | General utilities | (internal) |
| `internals_nodes.R` | Node coordinate calculation | (internal) |
| `internals_kij.R` | Spring constant models | (internal) |
| `mutenm-imports.R` | Package imports | - |
| `mutenm-package.R` | Package documentation | - |

### Data Flow

```
PDB file → bio3d::read.pdb() → pdb object
    ↓
enm(pdb, node, model, d_max) → prot object
    ↓
    ├── mutenm(wt, site_mut, ...) → mutant prot object
    │       ↓
    │       Compare using D-functions (Dr2i, Dv_min, etc.)
    │
    └── mrs(wt, nmut, ...) → response tibbles
```

### The `prot` Object

A list with class "prot" containing:
- `param`: ENM parameters (node, model, d_max)
- `nodes`: Site information (nsites, site, pdb_site, bfactor, xyz)
- `graph`: Tibble of edges (edge, i, j, v0ij, sdij, lij, kij, dij)
- `eij`: Unit vectors along edges (nedges × 3 matrix)
- `kmat`: Hessian matrix (3N × 3N)
- `nma`: Normal mode analysis (mode, evalue, umat, cmat)

## Exported Functions (20 total)

### Core (3)
- `enm(pdb, node, model, d_max)` - Build ENM
- `mutenm(wt, site_mut, mutation, mut_model, ...)` - Create mutant
- `mrs(wt, nmut, ...)` - Scan all sites

### Single Protein (9)
**Structure:** `cn`, `wcn`, `dactive`
**Dynamics:** `msfi`, `msfn`
**Energy:** `v_min`, `g_ent`, `dv_act`, `dg_ent_act`

### Pair Comparison (8)
**Structure:** `Dr2i`, `Dr2n`
**Dynamics:** `Dmsfi`, `Dmsfn`
**Energy:** `Dv_min`, `Dg_ent`, `Ddv_act`, `Ddg_ent_act`

## Dependencies

**External:**
- `bio3d` - PDB file reading, atom selection
- `pracma` - Cross product for Levitt's CB approximation
- `dplyr`, `tibble`, `tidyr`, `magrittr` - Data manipulation
- `stats` - Random number generation

## ENM Models Supported

| Model | Function | Description |
|-------|----------|-------------|
| anm | `kij_anm` | Standard ANM (k=1 for contacts) |
| ming_wall | `kij_ming_wall` | Ming & Wall 2005 (backbone penalty) |
| hnm | `kij_hnm` | Hinsen model |
| hnm0 | `kij_hnm0` | Exponential Hinsen |
| pfanm | `kij_pfanm` | Parameter-free ANM (1/d²) |
| reach | `kij_reach` | Reach et al. model |

## Node Types

| Type | Aliases | Description |
|------|---------|-------------|
| ca | calpha | Alpha carbon positions |
| sc | side_chain | Side chain center of mass |
| cb | beta | Beta carbon (or Levitt approx for Gly) |

## Mutation Models

| Model | Description |
|-------|-------------|
| lfenm | Linear Force ENM - only updates coordinates |
| sclfenm | Self-Consistent LFENM - recalculates full ENM |

## Testing

Tests in `tests/testthat/` with reference data in `tests/data/`:
- `test_enm.R` - ENM construction
- `test_mutenm.R` - Mutation (lfenm)
- `test_mutenm_sc.R` - Mutation (sclfenm)
- `test_mrs.R` - Mutation response scanning
- `test_enm_analysis.R` - Single protein properties
- `test_enm_energy.R` - Energy functions
- `test_delta_structure_motion.R` - Structure/dynamics comparison
- `test_delta_energy.R` - Energy comparison
