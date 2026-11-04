# Fig. 2 audit: how to fill `audit_rubric.csv`

One row per **figure panel x rater**. Two raters score every panel
independently; `01_fig2_audit.R` computes agreement (Cohen's kappa) and then
summarises the panels where the raters agree (disagreements are listed for a
consensus pass; enter the consensus with `rater = consensus`).

**Sampling (fix before scoring).** PubMed, one calendar year, journals chosen in
advance, articles whose figures show genome-browser tracks of ChIP-seq,
CUT&RUN/CUT&Tag or ATAC-seq at named example genes. Record the query and the
random seed used to draw the sample in the Methods.

| column | values | meaning |
|---|---|---|
| panel_id | text | unique id, e.g. `P017_Fig3C` |
| assay | ChIP / CUT&RUN / CUT&Tag / ATAC / other | |
| n_loci_shown | integer | example loci in the panel |
| selection_rule_stated | yes / partial / no | text says how the loci were chosen |
| quantification_shown | yes / no | genome-wide or group-level quantification next to the examples |
| normalization_stated | yes / no | track normalization stated (CPM, spike-in, ...) |
| global_change_addressed | yes / no / NA | if signal differs between conditions: is a genome-wide shift shown or ruled out? NA when one condition |
| peak_gene_rule_stated | yes / no | how peaks were linked to genes |
| replicates_shown | yes / no | individual replicates or replicate counts shown |
| qc_mentioned | yes / no | mappability / blacklist / signal-gap checks mentioned |
