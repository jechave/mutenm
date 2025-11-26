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

### 3.1 Analytical Method: `amrs(wt, mut_dl_sigma, mut_sd_min, option, response)`

Calculates mutation-response matrices using closed-form analytical formulas.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wt` | prot | - | Wild-type protein |
| `mut_dl_sigma` | numeric | - | Perturbation magnitude (Å) |
| `mut_sd_min` | integer | - | Minimum sequence distance |
| `option` | character | `"site"` | `"site"` (N×N matrix) or `"mode"` (M×N matrix) |
| `response` | character | `"dr2"` | Response type: `"dr2"`, `"de2"`, `"df2"`, `"stress"` |

**Returns:** Matrix where columns are mutated sites, rows are responding sites or modes.

### 3.2 Simulation Method: `smrs(wt, nmut, mut_dl_sigma, mut_sd_min, option, response, seed)`

Calculates mutation-response matrices using Monte Carlo simulation.

**Additional Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nmut` | integer | - | Number of mutations per site to simulate |
| `seed` | integer | 1024 | Random seed |

---

## 4. Response Types

### 4.1 Structural Responses

| Response | Matrix Element `M[i,j]` | Description |
|----------|-------------------------|-------------|
| `dr2` | `<(r_i^mut - r_i^wt)²>` | Mean squared displacement of site i, averaged over mutations at j |
| `de2` | `<½ k_i (δr_i)²>` | Deformation energy at site i |
| `df2` | `<f_i²>` | Squared force magnitude at site i |

### 4.2 Energy Responses (from `amrs_all`/`smrs_all`)

| Response | Description |
|----------|-------------|
| `dvsij` | Stress energy: `½ Σ k_ij × δl_ij²` at site i due to mutation at j |
| `dvmij` | Mode energy: `dvsij - de2ij` |

### 4.3 Profiles (sums over matrices)

| Profile | Calculation | Meaning |
|---------|-------------|---------|
| `dr2j`, `de2j`, `df2j`, `dvsj`, `dvmj` | Sum over i | **Influence**: total effect of mutating site j |
| `dr2i`, `de2i`, `df2i`, `dvsi`, `dvmi` | Mean over j | **Sensitivity**: average response of site i |

---

## 5. Dynamics/Motion Responses

These require `sclfenm` model (full ENM recalculation) to compare wt vs mutant ensembles.

### 5.1 Single-pair comparisons (`delta_motion_*` functions)

| Function | Returns | Description |
|----------|---------|-------------|
| `delta_motion_dmsfi(wt, mut)` | vector[i] | Change in mean-square fluctuations: `MSF_i^mut - MSF_i^wt` |
| `delta_motion_dhi(wt, mut)` | vector[i] | Change in site entropy |
| `delta_motion_dbhati(wt, mut)` | vector[i] | Bhattacharyya distance between site distributions |
| `delta_motion_rwsipi(wt, mut)` | vector[i] | Root-weighted square inner product (ensemble similarity) |

### 5.2 Motion response matrices (`mrs_motion_*` functions)

These take a tibble of mutants from `generate_mutants()`:

| Function | Matrix Element | Description |
|----------|----------------|-------------|
| `mrs_motion_dmsfij(mutants)` | `<ΔMSF_i>_j` | MSF change at i, averaged over mutations at j |
| `mrs_motion_dhij(mutants)` | `<Δh_i>_j` | Entropy change at i |
| `mrs_motion_dbhatij(mutants)` | `<d_bhat(i)>_j` | Bhattacharyya distance at i |
| `mrs_motion_rwsipij(mutants)` | `<rwsip_i>_j` | RWSIP similarity at i |

---

## 6. Stability Prediction (ΔΔG)

### 6.1 `amrs_ddg(wt, mut_model, mut_dl_sigma, mut_sd_min)`

Calculates stability change profile analytically.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `wt` | prot | Wild-type protein |
| `mut_model` | character | Must be `"lfenm"` (currently) |
| `mut_dl_sigma` | numeric | Perturbation magnitude |
| `mut_sd_min` | integer | Minimum sequence distance |

**Returns:** Tibble with columns: `site`, `pdb_site`, `ddg`

**Calculation:** `ΔΔG = ΔV_s - ΔE_def = ½Σ(dvsij - de2ij)`

### 6.2 Direct energy functions

| Function | Description |
|----------|-------------|
| `ddg_dv(wt, mut)` | Minimum energy difference |
| `ddg_tds(wt, mut, beta)` | Entropic free energy difference |
| `delta_energy_dvs(wt, mut, ideal)` | Stress energy difference |

---

## 7. Activation Energy (ΔΔG‡)

### 7.1 `amrs_ddgact(wt, pdb_site_active, mut_model, mut_dl_sigma, mut_sd_min)`

Calculates activation energy change profile.

**Additional Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_site_active` | integer vector | PDB residue numbers defining the active site |

**Returns:** Tibble with columns: `site`, `pdb_site`, `is_active`, `ddgact`

**Calculation:** Energy needed to deform mutant's active site back to wild-type conformation.

### 7.2 Direct activation energy functions

| Function | Description |
|----------|-------------|
| `ddgact_dv(wt, mut, ideal, pdb_site_active)` | Energy contribution to activation change |
| `ddgact_tds(wt, mut, ideal, pdb_site_active, beta)` | Entropy contribution to activation change |

---

## 8. Summary: What's in README vs What's Not

### Covered in README:
- Basic `set_enm()` with `node="calpha"`, `model="anm"`
- `amrs()` and `smrs()` for single-site scanning
- Structural responses: `dr2`, `de2`, `df2`
- `option="site"` and `option="mode"`

### NOT Covered in README:
- **Other node types:** `"sc"` (side chains), `"cb"` (beta carbons)
- **Other ENM models:** `"ming_wall"`, `"pfanm"`, `"hnm"`, `"hnm0"`, `"reach"`
- **sclfenm mutation model**
- **Energy responses:** `dvs`, `dvm`
- **Dynamics/motion responses:** `dmsf`, `dh`, `dbhat`, `rwsip`
- **Motion response matrices:** `mrs_motion_*` functions
- **Stability functions:** `amrs_ddg()`, `ddg_dv()`, `ddg_tds()`
- **Activation energy functions:** `amrs_ddgact()`, `ddgact_dv()`, `ddgact_tds()`
- **Comprehensive scanning:** `amrs_all()`, `smrs_all()`
- **`mut_sd_min` parameter explanation:** sequence distance filtering for edges

---

## 9. Key Equations

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
