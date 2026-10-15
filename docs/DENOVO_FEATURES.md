# De novo feature review

Additive features based on upstream history. Each dependent branch adds one reviewable module group.

## F01: TRGT parental-allele candidate hypotheses

Branch: `feature/01-trgt-candidates`. Depends on `main`.

Modules: `trio_alleles.py`, `trio_prepare.py`.

## F02: Read-local C3 evidence and complete molecule extraction

Branch: `feature/02-read-locus`. Depends on `feature/01-trgt-candidates`.

Modules: `trio_locus.py`, `trio_io.py`.

## F03: Complete local allele comparison and strict evidence review

Branch: `feature/03-local-alleles`. Depends on `feature/02-read-locus`.

Modules: `trio_local.py`, `trio_finalize.py`, `trio_review.py`.
