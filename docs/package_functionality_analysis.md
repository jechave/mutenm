# mutenm Package Functionality Analysis

This document provides a comprehensive analysis of what the mutenm package can do, based on examination of the source code.

---

## 1. ENM Construction

### `enm(pdb, node, model, d_max)`

Creates a `prot` object containing the Elastic Network Model.

**Parameters:**

| Parameter | Type | Description | Values |
|-----------|------|-------------|--------|
| `pdb` | bio3d object | PDB structure from `bio3d::read.pdb()` | - |
| `node` | character | Network node placement | `"ca"`, `"calpha"` (α-carbons), `"sc"`, `"side_chain"` (side chain centers), `"cb"`, `"beta"` (β-carbons) |
| `model` | character | ENM model type | `"anm"`, `"ming_wall"`, `"pfanm"`, `"hnm"`, `"hnm0"`, `"reach"` |
| `d_max` | numeric | Distance cutoff for contacts (Å) | Typical: 10.5 for Cα, 12.5 for side chains |

**Output (`prot` object components):**

| Component | Description |
|-----------|-------------|
| `param` | List: `node`, `model`, `d_max` |
| `nodes` | List: `nsites`, `site` (sequential 1:N), `pdb_site` (PDB residue numbers), `bfactor`, `xyz` (coordinates) |
| `graph` | Tibble: `edge`, `i`, `j`, `v0ij`, `sdij` (sequence distance), `lij` (equilibrium length), `kij` (spring constant), `dij` (actual distance) |
| `eij` | Matrix (n_edges × 3): unit vectors along contacts |
| `kmat` | Matrix (3N × 3N): Kirchhoff/Hessian matrix |
| `nma` | List: `mode`, `evalue` (eigenvalues), `umat` (eigenvectors), `cmat` (covariance matrix) |

---

## 2. Mutation Perturbation

### `mutenm(wt, site_mut, mutation, mut_model, mut_dl_sigma, mut_sd_min, seed)`

Core function for introducing mutations and calculating structural responses. (Called by `mrs()`, also exported for direct use.)

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wt` | prot | - | Wild-type protein object |
| `site_mut` | integer | - | Site to mutate (sequential index, not PDB number) |
| `mutation` | integer | 0 | Mutation index (0 = no mutation, returns wt) |
| `mut_model` | character | `"lfenm"` | Mutation model: `"lfenm"` (fast, linear) or `"sclfenm"` (self-consistent, recalculates ENM) |
| `mut_dl_sigma` | numeric | 0.3 | Standard deviation for edge length perturbations (Å) |
| `mut_sd_min` | integer | 2 | Minimum sequence distance for edges to be perturbed |
| `seed` | integer | 241956 | Random seed for reproducibility |

**Mutation Models:**

| Model | Description | Speed | Accuracy |
|-------|-------------|-------|----------|
| `lfenm` | Linear Force ENM - uses linear response approximation | Fast | Approximate |
| `sclfenm` | Self-Consistent LFENM - recalculates full ENM after perturbation | Slow | More accurate |

**How mutations work:**
1. Edges connected to `site_mut` with sequence distance ≥ `mut_sd_min` are selected
2. Each selected edge length is perturbed by `δl ~ N(0, mut_dl_sigma)`
3. Forces are calculated: `f_ij = -k_ij × δl_ij`
4. Structural response: `δr = C × f` (using covariance matrix)

---

## 3. Mutation Response Scanning

### `mrs(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, responses, seed)`

Memory-efficient function for calculating mutation-response matrices. Uses process-and-discard approach: generates one mutant at a time, calculates responses, discards immediately.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wt` | prot | - | Wild-type protein |
| `nmut` | integer | - | Number of mutations per site to simulate |
| `mut_model` | character | `"lfenm"` | Mutation model |
| `mut_dl_sigma` | numeric | 0.3 | Perturbation magnitude (Å) |
| `mut_sd_min` | integer | 2 | Minimum sequence distance |
| `responses` | character | `"dr2ij"` | Vector of responses to calculate |
| `seed` | integer | NULL | Random seed |

**Available responses:**

| Response | Type | Description | Requires |
|----------|------|-------------|----------|
| `dr2ij` | site matrix (N×N) | Squared displacement at site i due to mutation at j | lfenm or sclfenm |
| `dr2nj` | mode matrix (M×N) | Squared displacement in mode n due to mutation at j | lfenm or sclfenm |
| `dmsfij` | site matrix (N×N) | MSF change at site i due to mutation at j | sclfenm |
| `dmsfnj` | mode matrix (M×N) | MSF change in mode n due to mutation at j | sclfenm |
| `ddg_dv` | profile (N) | Minimum energy change for mutations at j | lfenm or sclfenm |
| `ddg_tds` | profile (N) | Entropic free energy change for mutations at j | lfenm or sclfenm |
| `Ddv_act` | profile (N) | Activation energy change (energy part) | lfenm or sclfenm |
| `Ddg_ent_act` | profile (N) | Activation energy change (entropy part) | lfenm or sclfenm |

**Returns:** List with:
- Requested site matrices as tibbles: `i`, `j`, `<response>`
- Requested mode matrices as tibbles: `n`, `j`, `<response>`
- Requested scalar profiles as tibbles: `j`, `<response>`
- `$params` — ENM and mutation parameters for reproducibility

---

## 4. Pairwise Comparison Functions

These functions compare a single wild-type/mutant pair. Used internally by `mrs()`.

### 4.1 Structural comparisons

| Function | Returns | Description |
|----------|---------|-------------|
| `Dr2i(wt, mut)` | vector[i] | Squared displacement per site |
| `Dr2n(wt, mut)` | vector[n] | Squared displacement per mode |

### 4.2 Motion comparisons

| Function | Returns | Description |
|----------|---------|-------------|
| `Dmsfi(wt, mut)` | vector[i] | MSF change per site |
| `Dmsfn(wt, mut)` | vector[n] | MSF change per mode |

---

## 5. Energy Functions

### 5.1 ENM energy calculations

| Function | Description |
|----------|-------------|
| `v_min(prot)` | Minimum (stress) energy of the ENM |
| `g_ent(prot, beta)` | Entropic free energy contribution |

### 5.2 Pairwise energy differences (ΔΔG)

| Function | Description |
|----------|-------------|
| `Dv_min(wt, mut)` | Minimum energy difference |
| `Dg_ent(wt, mut, beta)` | Entropic free energy difference |

### 5.3 Activation energy functions (ΔΔG‡)

| Function | Description |
|----------|-------------|
| `dv_act(prot, ideal, pdb_site_active)` | Energy contribution to activation |
| `dg_ent_act(prot, ideal, pdb_site_active, beta)` | Entropy contribution to activation |
| `Ddv_act(wt, mut, ideal, pdb_site_active)` | Activation energy difference (energy part) |
| `Ddg_ent_act(wt, mut, ideal, pdb_site_active, beta)` | Activation energy difference (entropy part) |

---

## 6. Analysis Functions (Exported)

### 6.1 Site profiles

| Function | Returns |
|----------|---------|
| `cn(prot)` | Contact number profile |
| `wcn(prot)` | Weighted contact number profile |
| `dactive(prot, pdb_site_active)` | Distance to active site profile |
| `msfi(prot)` | Mean-square fluctuation per site |
| `get_mlms(prot, sdij_cut)` | Mean local mutational stress profile |

### 6.2 Mode profiles

| Function | Returns |
|----------|---------|
| `msfn(prot)` | MSF per mode (1/λ) |

---

## 7. Key Equations

### Linear Response
```
δr = C × f
```
where `C` is the covariance matrix, `f` is the force vector from edge perturbations.

### Force from perturbation
```
f_ij = -k_ij × δl_ij × ê_ij
```

### Minimum energy
```
V_min = Σ v0_ij + ½ Σ k_ij × (d_ij - l_ij)²
```

### Stress energy
```
V_s = ½ Σ k_ij × (d_ij^ideal - l_ij)²
```

---

*Generated for mutenm package development reference*
