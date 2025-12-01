# mutenm Function Grouping (Alternative)

This document organizes the 20 exported functions by: Core, Single Protein, Pair Comparison.

---

## Quick Reference

### Core (3 functions)

| Function | Description |
|----------|-------------|
| `enm` | Build ENM from PDB |
| `mutenm` | Create single mutant |
| `mrs` | Mutation response scanning |

### Single Protein (9 functions)

Functions that take one `prot` object and return a property.

#### Structure

| Function | Description |
|----------|-------------|
| `cn` | Contact number |
| `wcn` | Weighted contact number |
| `dactive` | Distance to active site |

#### Dynamics

| Function | Description |
|----------|-------------|
| `msfi` | Mean-square fluctuation per site |
| `msfn` | Mean-square fluctuation per mode |

#### Energy

| Function | Description |
|----------|-------------|
| `v_min` | Minimum (stress) energy |
| `g_ent` | Entropic free energy |
| `dv_act` | Activation energy (internal) |
| `dg_ent_act` | Activation entropy |

### Pair Comparison (8 functions)

Functions that compare wild-type and mutant: `Df(wt, mut)`.

#### Structure

| Function | Description |
|----------|-------------|
| `Dr2i` | Squared displacement per site |
| `Dr2n` | Squared displacement per mode |

#### Dynamics

| Function | Description |
|----------|-------------|
| `Dmsfi` | Change in MSF per site |
| `Dmsfn` | Change in MSF per mode |

#### Energy

| Function | Description |
|----------|-------------|
| `Dv_min` | Change in minimum energy |
| `Dg_ent` | Change in entropic free energy |
| `Ddv_act` | Change in activation energy |
| `Ddg_ent_act` | Change in activation entropy |

---

## Rationale

1. **Core** — functions that build/mutate models
2. **Single Protein** — properties of one structure
3. **Pair Comparison** — properties that require wt and mut

Within Single Protein and Pair Comparison, functions are grouped by physical concept (structure, dynamics, energy).
