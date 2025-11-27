# mutenm Package Functionality Analysis

This document provides a comprehensive analysis of what the mutenm package can do, based on examination of the source code.

---

## 1. ENM Construction

### `set_enm(pdb, node, model, d_max)`

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

### `get_mutant_site(wt, site_mut, mutation, mut_model, mut_dl_sigma, mut_sd_min, seed)`

Core function for introducing mutations and calculating structural responses.

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
| `dvs` | profile (N) | Stress energy change for mutations at j | lfenm or sclfenm |
| `ddgact_dv` | profile (N) | Activation energy change (energy part) | lfenm or sclfenm |
| `ddgact_tds` | profile (N) | Activation energy change (entropy part) | lfenm or sclfenm |

**Returns:** List with:
- Requested site matrices as tibbles: `i`, `j`, `<response>`
- Requested mode matrices as tibbles: `n`, `j`, `<response>`
- Requested scalar profiles as tibbles: `j`, `<response>`
- `$params` — ENM and mutation parameters for reproducibility

---

## 4. Response Types

### 4.1 Structural Responses (available via `mrs()`)

| Response | Matrix Element `M[i,j]` | Description |
|----------|-------------------------|-------------|
| `dr2ij` | `<(r_i^mut - r_i^wt)²>` | Mean squared displacement of site i, averaged over mutations at j |
| `dr2nj` | `<(q_n^mut - q_n^wt)²>` | Mean squared displacement in mode n, averaged over mutations at j |

### 4.2 Motion Responses (require sclfenm, available via `mrs()`)

| Response | Matrix Element | Description |
|----------|----------------|-------------|
| `dmsfij` | `<ΔMSF_i>_j` | MSF change at site i due to mutations at j |
| `dmsfnj` | `<ΔMSF_n>_j` | MSF change in mode n due to mutations at j |

### 4.3 Energy Profiles (available via `mrs()`)

| Response | Description |
|----------|-------------|
| `ddg_dv` | Minimum energy change (per mutation site) |
| `ddg_tds` | Entropic free energy change (per mutation site) |
| `dvs` | Stress energy change (per mutation site) |

---

## 5. Pairwise Comparison Functions

These functions compare a single wild-type/mutant pair. Used internally by `mrs()`.

### 5.1 Structural comparisons (`delta_structure_*`)

| Function | Returns | Description |
|----------|---------|-------------|
| `delta_structure_dr2i(wt, mut)` | vector[i] | Squared displacement per site |
| `delta_structure_dr2n(wt, mut)` | vector[n] | Squared displacement per mode |

### 5.2 Motion comparisons (`delta_motion_*`)

| Function | Returns | Description |
|----------|---------|-------------|
| `delta_motion_dmsfi(wt, mut)` | vector[i] | MSF change per site |
| `delta_motion_dmsfn(wt, mut)` | vector[n] | MSF change per mode |

---

## 6. Stability Prediction (ΔΔG)

### 6.1 Pairwise energy functions

| Function | Description |
|----------|-------------|
| `ddg_dv(wt, mut)` | Minimum energy difference |
| `ddg_tds(wt, mut, beta)` | Entropic free energy difference |
| `delta_energy_dvs(wt, mut, ideal)` | Stress energy difference |

---

## 7. Activation Energy (ΔΔG‡)

### 7.1 Pairwise activation energy functions

| Function | Description |
|----------|-------------|
| `ddgact_dv(wt, mut, ideal, pdb_site_active)` | Energy contribution to activation change |
| `ddgact_tds(wt, mut, ideal, pdb_site_active, beta)` | Entropy contribution to activation change |

---

## 8. Key Equations

### Linear Response
```
δr = C × f
```
where `C` is the covariance matrix, `f` is the force vector from edge perturbations.

### Force from perturbation
```
f_ij = -k_ij × δl_ij × ê_ij
```

### Strain energy
```
V_s = ½ Σ k_ij × δl_ij²
```

### Stability change (ΔΔG)
```
ΔΔG ≈ ΔV_s - ΔE_def = ½ Σ k_ij δl_ij² - ½ Σ k_i (δr_i)²
```

---

*Generated for mutenm package development reference*
