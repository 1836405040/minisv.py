# Structural junction review

```bash
minisv denovo structural --manifest family.json --output results/structural1 --workers 16 --rate 0.08
```

The manifest needs `reference`, `father`, `mother`, and `children`. Each member
has `sample` and an indexed `bam`; each child has `candidate_vcf` and/or
`candidate_vcfs`. All paths are explicit. These are caller inputs, not truth
coordinates. Reuse original caller VCFs to recover BND records removed by older
adapters. The full `denovo run` command also invokes this reviewer automatically.
No TRGT, parental assemblies, or C3 cache is required by the standalone command.

Set `structural_exclude_bed` to a BED mask on the same reference (for example
`data/chm13v2.cen-mask.bed`). Nonautosomal events are excluded. For intrachromosomal
events, overlap of the breakpoint span with this BED excludes the event; for
interchromosomal events, each endpoint is checked. Without an explicit mask,
autosomal centromeric intervals are not excluded automatically.

## Models and evidence

* INV: both reference-relative inversion junctions must be reproduced.
* DUP: test the tandem-duplication junction. A symbolic DUP does not identify
  the destination of a dispersed copy; copy-number and destination evidence
  remain unresolved.
* BND: parse all four bracket orientations using zero-based interbase positions.
  Canonical endpoint pairs deduplicate exact reciprocal records. Mate links must
  resolve within the same source and geometry. Caller IDs survive preparation.
  Inserted sequence or complex BND anchors are explicitly unresolved.

Each junction joins 200 reference bases from each retained side. The reviewer
extracts complete available primary/supplementary read sequences near endpoints,
requires MAPQ >=20, excludes secondary/duplicate/QC-failed records, and counts
query names once per junction. It compares both sequence orientations using
semi-global edit distance. Alternative-specific evidence requires edit rate
<=0.08 (configurable) and at least 10 edits of advantage over both reference
contexts. Child evidence needs at least 3 molecules and both sequence orientations.
Both parental reference contexts need at least 10 specific molecules per parent;
any parental alternative support is a conflict. Low-discrimination junctions,
imprecise breakpoints, and exact junctions shared by siblings retain explicit flags.

The current model is experimental: junction support alone does not establish
parental haplotype callability or an independently validated de novo mutation.
No INV/DUP/BND is automatically promoted to high confidence. A fully supported
junction can enter REVIEW_REQUIRED; conflicts or missing evidence enter
INSUFFICIENT_OR_CONFLICTING_EVIDENCE. This policy does not change INS/DEL tiers.

## Outputs and reproducibility

`structural_input_records.jsonl` (standalone) or `structural_pending.jsonl`
(full pipeline) preserves input records. `structural_reviewed.jsonl` contains
junction models and molecule names for each family member. `structural_unresolved.jsonl`
retains unsupported, filtered, or out-of-scope records. Exact canonical junction
geometry is deduplicated; nearby breakpoint representations are not fuzzy-merged.

Balanced events use `delta: null`; BND also has `length: null`. Use `junctions`
for both chromosomes, positions and retained sides rather than interpreting
`start/end` as a single genomic interval. `variant_class` identifies structural
events; it does not represent independent de novo confirmation.

`summary.json`, `status.json` and `provenance.json` record counts and run inputs.
Output directories must be new. Synthetic tests cover orientations, paired BND
identity, inherited events, missing parent coverage, incomplete inversion support,
region masks, parser boundaries, gzip, and integrated INS/DEL nonregression.
