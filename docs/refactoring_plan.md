# Function Renaming Plan

## Naming Convention

- `D` prefix = comparison between two proteins (wt vs mut)
- `d` prefix = inherent difference within a single protein (e.g., activation energy)
- Underscore separates components (e.g., `v_min`, `dv_act`, `g_ent`)

## Renames

| Current | New | Status |
|---------|-----|--------|
| `enm_v_min` | `v_min` | ✓ done |
| `delta_energy_dv` | `Dv_min` | ✓ done |
| `enm_g_entropy` | `g_ent` | pending |
| `delta_energy_tds` | `Dg_ent` | pending |
| `dgact_dv` | `dv_act` | pending |
| `dgact_tds` | `dg_ent_act` | pending |
| `delta_energy_act_dv` | `Ddv_act` | pending |
| `delta_energy_act_tds` | `Ddg_ent_act` | pending |
| `get_msf_site` | `msfi` | pending |
| `get_msf_mode` | `msfn` | pending |
| `delta_motion_dmsfi` | `Dmsfi` | pending |
| `delta_motion_dmsfn` | `Dmsfn` | pending |
| `delta_structure_dr2i` | `Dr2i` | pending |
| `delta_structure_dr2n` | `Dr2n` | pending |

**Deferred:** stress-related functions (`calculate_vs`, `delta_energy_dvs`, `get_stress`)

## Process for each rename

0. Check there're appropriate tests for the functions to be renamed
1. Rename function definition in R/ file
2. Update roxygen `@details` if it mentions the function name
3. Update all callers (grep in R/ and tests/)
4. Run `devtools::document()`
5. Run `devtools::test()`
6. Commit
7. Update docs/ files if needed (package_dependencies.md, package_functionality_analysis.md)

## Suggested groupings (families)

- **g_ent family**: `enm_g_entropy` → `g_ent`, `delta_energy_tds` → `Dg_ent`
- **dv_act family**: `dgact_dv` → `dv_act`, `delta_energy_act_dv` → `Ddv_act`
- **dg_ent_act family**: `dgact_tds` → `dg_ent_act`, `delta_energy_act_tds` → `Ddg_ent_act`
- **msf family**: `get_msf_site` → `msfi`, `get_msf_mode` → `msfn`, `delta_motion_dmsfi` → `Dmsfi`, `delta_motion_dmsfn` → `Dmsfn`
- **structure family**: `delta_structure_dr2i` → `Dr2i`, `delta_structure_dr2n` → `Dr2n`

## Notes

- 2acy active site residues (PDB numbering): Arg23, Asn41 → `pdb_site_active = c(23, 41)`
- lfenm mutants have same eigenvalues as wt (entropy diff = 0)
- sclfenm mutants recalculate eigenvalues (entropy diff ≠ 0)
