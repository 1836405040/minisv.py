# Family de novo SV and indel workflow

```bash
minisv denovo init --father dad --mother mom --child kid -o family.json
# Replace paths with prepared family inputs.
minisv denovo check -m family.json
minisv denovo run --manifest family.json --output results/run1 --workers 16 --error-rate 0.08
minisv denovo export -i results/run1/all_reviewed.jsonl -o results/run1/events.tsv
```

The manifest supplies the reference, indexed family HiFi BAMs, four parental
haplotype FASTAs, caller VCFs with RNAMES, and TRGT genotype/spanning files.
Existing local evidence and C3 files may be supplied explicitly. Caller and
assembly execution are upstream inputs, not launched by this command.

Set `minimum_variant_length` to 30 to include 30–49 bp de novo indels;
the default remains 50 for SVs. Output `variant_class` distinguishes these
classes. `basis` distinguishes reference-relative and parent-relative events.
Do not combine their lengths without sequence-aware reconciliation.

Only `confidence_policy: strict` is supported. The optional
`parent_readability_policy` is `marker_assisted` by default, or `strict` for
the original audit. Concrete INS review requires `parent_haplotype_alignments`
with F0/F1/M0/M1 entries, each with `fasta` and reference `paf` including cs tags.
Marker-assisted review additionally requires parental BAM HP/PS tags.
Absent phase evidence leaves the original unresolved classification intact.

The workflow requires child-specific molecules, parental coverage, target
absence in parents, and branch-specific allele readability. C3 alone cannot
promote a candidate. The marker fallback preserves every other conflict and
requires directly linked markers and matching bidirectional consensuses.

Outputs include `strict_reviewed.jsonl`, `marker_assisted_reviewed.jsonl`,
`all_reviewed.jsonl`, tier-specific JSONL files, `summary.json`, `status.json`,
and `provenance.json`. A failed run retains its logs; use a new output directory
for a rerun. Input evidence is not modified. Reference coordinates are 0-based
half-open in event JSON; insertion boundaries have start=end.

Unit and integration tests run separately in GitHub Actions. Integration tests
generate synthetic BAM/FASTA/VCF data and execute actual minimap2 alignment.
No patient data or truth coordinates are required for discovery.

INV/DUP/BND use the [structural junction review](STRUCTURAL.md), with separate
outputs and explicit unresolved reasons. Use original caller VCFs when older
preparation steps removed these types.
