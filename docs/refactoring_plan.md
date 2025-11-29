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
| `enm_g_entropy` | `g_ent` | ✓ done |
| `enm_g_entropy_mode` | `g_ent_mode` | ✓ done |
| `delta_energy_tds` | `Dg_ent` | ✓ done |
| `dgact_dv` | `dv_act` | ✓ done |
| `dgact_tds` | `dg_ent_act` | ✓ done |
| `delta_energy_act_dv` | `Ddv_act` | ✓ done |
| `delta_energy_act_tds` | `Ddg_ent_act` | ✓ done |
| `get_msf_site` | `msfi` | ✓ done |
| `get_msf_mode` | `msfn` | ✓ done |
| `delta_motion_dmsfi` | `Dmsfi` | ✓ done |
| `delta_motion_dmsfn` | `Dmsfn` | ✓ done |
| `delta_structure_dr2i` | `Dr2i` | ✓ done |
| `delta_structure_dr2n` | `Dr2n` | ✓ done |

**Deferred:** stress-related functions (`calculate_vs`, `delta_energy_dvs`, `get_stress`)

## Process for each rename

0. Check there're appropriate tests for the functions to be renamed
1. Rename function definition in R/ file
2. Update roxygen `@details` if it mentions the function name
3. Update all callers (grep in R/ and tests/)
4. Update docs/ files (package_dependencies.md, package_functionality_analysis.md)
5. Run `devtools::document()`
6. Run `devtools::test()`
7. Grep for any remaining occurrences of old name
8. Single commit with all changes (code + docs together)

## Suggested groupings (families)

- ~~**g_ent family**: `enm_g_entropy` → `g_ent`, `enm_g_entropy_mode` → `g_ent_mode`, `delta_energy_tds` → `Dg_ent`~~ **DONE**
- ~~**dv_act family**: `dgact_dv` → `dv_act`, `delta_energy_act_dv` → `Ddv_act`~~ **DONE**
- ~~**dg_ent_act family**: `dgact_tds` → `dg_ent_act`, `delta_energy_act_tds` → `Ddg_ent_act`~~ **DONE**
- ~~**msf family**: `get_msf_site` → `msfi`, `get_msf_mode` → `msfn`, `delta_motion_dmsfi` → `Dmsfi`, `delta_motion_dmsfn` → `Dmsfn`~~ **DONE**
- ~~**structure family**: `delta_structure_dr2i` → `Dr2i`, `delta_structure_dr2n` → `Dr2n`~~ **DONE**

## Notes

- 2acy active site residues (PDB numbering): Arg23, Asn41 → `pdb_site_active = c(23, 41)`
- lfenm mutants have same eigenvalues as wt (entropy diff = 0)
- sclfenm mutants recalculate eigenvalues (entropy diff ≠ 0)
