# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Context Files

Before starting work, read:
- `docs/package_functionality_analysis.md` — what the package can do
- `docs/package_dependencies.md` — code structure and dependencies

## Package Overview

The `mutenm` package (Mutate Elastic Network Models) is an R package for building Elastic Network Models (ENMs) of proteins and calculating mutation-response matrices using the LFENM perturbation model.

**Note:** This is a development version being iteratively trimmed from the full `penm` package. The final scope will emerge through the trimming process itself.

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

### Current Components

Components are being evaluated during the trimming process. See the R/ directory for current code.

### Data Flow

1. PDB file → `bio3d::read.pdb()` → pdb object
2. pdb object → `set_enm()` → prot object
3. prot object → scanning functions → response matrices

### Key Dependencies

- **bio3d**: PDB file reading
- **jefuns**: Utility functions
- **pracma**: Numerical analysis (cross product)

### Testing

Tests in `tests/testthat/` with reference data in `tests/data/`.
Run `Rscript -e "devtools::test()"` after changes.

## Trimming Process

This package is being trimmed iteratively. Decisions about what to keep or remove are made step by step as the process unfolds.
