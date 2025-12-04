# JOSS Submission Requirements

## Paper Requirements
- **Format**: `paper.md` in the repository
- **Content must include**:
  - Title, author names, affiliations
  - Summary (high-level functionality for non-specialists)
  - Statement of need (problems solved, target audience, comparison to similar tools)
  - Key references
- **Word count**: 250-1000 words (short papers)
- **Should NOT include**: API documentation or detailed function lists

## Software Requirements

| Requirement | Status for mutenm |
|-------------|-------------------|
| OSI-approved license (actual LICENSE file) | Need to verify |
| >1000 lines of code (or justify if less) | Need to check |
| Hosted on clonable repository | ✓ GitHub |
| Issue tracker readable without registration | ✓ GitHub |
| pip/CRAN installable or documented installation | ✓ devtools::install_github() |

## Documentation Requirements
- Installation instructions | ✓ README
- Example usage | ✓ README + vignettes
- API documentation | ✓ roxygen docs
- Community guidelines (how to contribute/report issues) | May need to add

## Tests
- Automated test suite with CI preferred | ✓ testthat tests (CI not set up)

## Substantial Scholarly Effort
- Minimum ~3 months individual work
- Flagged if <1000 LOC, rejected if <300 LOC

---

## Review Criteria Checklist

### Paper
- [ ] Author names and affiliations
- [ ] High-level functionality summary for non-specialists
- [ ] Clear statement of need
- [ ] Comparison to similar packages
- [ ] References to research using the software
- [ ] Key references with software archive link

### Software License
- [ ] OSI-approved license in repository
- [ ] Actual LICENSE or COPYING file present (not just "MIT license" in README)

### Documentation
- [ ] Statement of need (problems solved, target audience)
- [ ] Installation instructions with dependencies
- [ ] Example usage demonstrating real-world problems
- [ ] API documentation (functions documented with examples)
- [ ] Community guidelines (contribution, issue reporting, support)

### Functionality
- [ ] Reviewers can install and verify core functionality works

### Tests
- [ ] Automated test suite with CI (good)
- [ ] Or documented manual testing steps (OK)

---

**Note**: The review criteria state the paper "should not include software documentation such as API functionality." Function lists in the paper may be seen as too detailed.
