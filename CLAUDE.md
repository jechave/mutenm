# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

The `mutenm` package (Mutate Elastic Network Models) is an R package for building Elastic Network Models (ENMs) of proteins and calculating mutation-response matrices using the LFENM perturbation model.

**Note:** This is a development version being trimmed from the full `penm` package to create a minimal, focused release.

## Development Commands

### Building and Installing
```bash
R CMD INSTALL .
R CMD build .
R CMD check mutenm_*.tar.gz

# Using devtools (recommended)
Rscript -e "devtools::install()"
Rscript -e "devtools::check()"
```

### Testing
```bash
Rscript -e "devtools::test()"
Rscript -e "testthat::test_file('tests/testthat/test_enm.R')"
```

### Documentation
```bash
Rscript -e "devtools::document()"
```

## Architecture

### Core Components (to keep in minimal release)

1. **ENM Creation (`R/enm.R`)**
   - `set_enm()`: Creates a `prot` object containing ENM structure
   - Target models for v1.0: "anm", "ming_wall", "pfanm", "hnm"
   - Node types: "ca", "cb", "sc"

2. **Mutation Perturbation (`R/penm.R`)**
   - `get_mutant_site()`: Core mutation function
   - Target for v1.0: `lfenm` model only

3. **Mutation Response Scanning**
   - `amrs()`: Analytical single-site mutation response
   - `smrs()`: Simulation-based single-site mutation response

4. **Response Types**: `dr2`, `de2`, `df2`

### To be removed during trimming

- `sclfenm` mutation model
- Models: "hnm0", "reach"
- `amrs_ddg`, `smrs_ddg`, `amrs_ddgact`, `smrs_ddgact`

### Data Flow

1. PDB file → `bio3d::read.pdb()` → pdb object
2. pdb object → `set_enm()` → prot object
3. prot object → scanning functions → response matrices

### Key Dependencies

- **bio3d**: PDB file reading
- **jefuns**: Utility functions
- **Matrix**: Sparse matrix operations
- **pracma**: Numerical analysis

### Testing

Tests in `tests/testthat/` with reference data in `tests/data/`.
Run `Rscript -e "devtools::test()"` after changes.

## Important: Keeping CLAUDE.md in Sync

When making changes to the package, always update this file in the same commit:
- Adding or removing exported functions: Update the Architecture section
- Changing function signatures or API: Update relevant documentation
- Removing functionality during trimming: Remove references from all sections

This ensures CLAUDE.md stays accurate and useful for future development sessions.
