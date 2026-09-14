# Gene Expression Analysis in Alzheimer's Disease

**A complete RNA-seq differential expression workflow identifying genes and a candidate hub 
gene associated with Alzheimer's disease, built from a public dataset (GEO: GSE104704).**

This repository is self-contained: everything needed to understand the project — question, 
methods, results, code, and limitations — is documented below. The full written report 
(PDF) is also included for a more detailed narrative treatment of the same analysis.

---

## Project Summary

| | |
|---|---|
| **Question** | Which genes are differentially expressed between Alzheimer's disease (AD) and age-matched healthy brain tissue, and what might these changes suggest about biological processes associated with neurodegeneration and candidate genes for further investigation? |
| **Dataset** | [GSE104704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704) — human lateral temporal lobe RNA-seq (10 aged healthy vs. 12 AD samples) |
| **Method** | Differential expression (limma) → functional enrichment (g:Profiler) → protein interaction network (STRING) → candidate gene validation (UniProt, PDB, literature) |
| **Key result** | 179 significant DEGs (73 up, 106 down); **BDNF identified as the top candidate hub gene for further investigation** — not a confirmed therapeutic target — supported by network connectivity and independent published literature |
| **Full report** | See `AD_Gene_Expression_Report.pdf` in this repository for the complete written analysis with all figures |

---

## Repository Structure
- `scripts/` — R code:
  - `load_data.R` — loads the expression matrix and assigns sample groups
  - `phase4_limma_analysis.R` — differential expression analysis (limma), volcano plot, heatmap
  - `phase3_qc.R` — quality control: PCA, sample correlation heatmap, outlier sensitivity check
- `figures/` — Volcano plot, DEG heatmap, PCA plot, correlation heatmap, STRING network image, candidate protein structure
- `results/` — DEG list, enrichment results (g:Profiler), and STRING interaction/annotation tables
- `archive_first_attempt/` — An initial DESeq2-based analysis, preserved for transparency (see Methodological Note below)

---

## Methods

**Dataset.** RNA-seq gene expression data was obtained from NCBI GEO accession GSE104704 
(lateral temporal lobe, human brain). Sample-to-group labels were verified against GEO 
sample metadata by matching individual sample IDs directly, rather than assuming column 
order (an initial mismatch was caught and corrected during this step).

**Data type and normalization.** The original study (Nativio et al., 2018) generated raw 
counts via FeatureCounts and ran DESeq2 internally. The publicly deposited supplementary 
file, however, represents processed (non-raw) expression data — confirmed by its 
non-integer values and by GEO's own note that raw data are hosted separately via SRA. This 
is why the present reanalysis used **limma**, appropriate for normalized continuous data, 
rather than DESeq2, which requires raw counts.

**Quality control.** Expression data was filtered to remove low-expression genes (mean 
log2 expression ≤ 1), retaining 21,563 of 27,130 genes. PCA and a sample correlation 
heatmap were used to check sample structure. One sample (X11.10T.RNA, Old group) was 
flagged as an outlier by both methods. A sensitivity check confirmed the main findings 
were not dependent on this sample (88% of original DEGs remained significant after its 
removal) — see full report for details.

**Differential expression.** Analyzed using limma's moderated t-statistic with empirical 
Bayes variance shrinkage (eBayes). Significance threshold: adjusted p < 0.05, |log2FC| > 1.

**Functional enrichment & network analysis.** Significant genes were submitted to 
g:Profiler (GO/KEGG) and STRING (protein interaction network). Hub gene status (network 
connectivity) was treated as supporting evidence for further investigation, not as proof 
of disease relevance on its own.

**Candidate validation.** The top hub gene (BDNF) was checked against UniProt (function), 
RCSB PDB (structure), and published literature for independent support.

---

## Methodological Note: DESeq2 → limma Correction

An initial differential expression analysis was performed using DESeq2. This was later 
identified as methodologically incorrect: DESeq2 requires raw integer counts, but the 
GSE104704 supplementary file contains pre-normalized, continuous values. Using normalized 
data with DESeq2 violates its statistical model and produces unreliable significance 
estimates.

This was corrected by switching to limma. The original DESeq2-based analysis is preserved 
in `archive_first_attempt/` for transparency, alongside the corrected limma-based analysis 
in `scripts/` and `results/`.

---

## Key Findings

- **179 significant DEGs**: 73 upregulated, 106 downregulated in AD
- Most significant genes: SGO1, KLF15 (up); RPH3A, VGF, CRH, NEUROD6 (down)
- Enrichment highlighted immune/receptor-related and neuronal signaling processes
- **BDNF** identified as the top hub gene in the STRING network (14 connections, more than 
  double any other gene) — **a candidate hub gene warranting further investigation**, 
  supported independently by published literature on reduced BDNF in AD brain tissue 
  (Phillips et al., 1991; Connor et al., 1997). **BDNF should not be interpreted as a 
  confirmed therapeutic target based on this analysis alone** — network connectivity and 
  differential expression indicate association, not a demonstrated causal or therapeutic role.

---

## Limitations

- **Small sample size** (10 Old vs. 12 AD) — findings are preliminary, not definitive
- **Secondary analysis** of existing public data, not new wet-lab samples
- **Unadjusted confounders**: RNA integrity, post-mortem interval, sex, and sequencing 
  batch were not incorporated. Metadata retrieval via GEOquery was attempted but blocked 
  by a system-level security restriction; group balance on these variables could not be 
  confirmed. This is the most significant open limitation of this analysis.
- **Enrichment background**: g:Profiler was run against the default whole-genome 
  background rather than an expression-restricted background, due to a technical issue 
  reproducing the custom-background configuration — may inflate apparent term significance
- **No experimental validation** — BDNF was identified computationally; no functional, 
  cellular, or animal experiments were performed
- **Association, not causation** — see above
- **Bulk tissue analysis** — expression reflects an average across cell types; some signal 
  may reflect cell-composition shifts rather than within-cell-type change
- **Outlier sample** (X11.10T.RNA) retained after a sensitivity check showed it did not 
  drive the main findings

*(Full detail on each point is in the written report, `AD_Gene_Expression_Report.pdf`.)*

---

## Software and Package Versions
- R version 4.6.1 (2026-06-24), Windows 11 x64
- limma 3.68.4
- ggplot2 4.0.3
- pheatmap 1.0.13

## References
Connor, B., Young, D., Yan, Q., Faull, R.L., Synek, B., & Dragunow, M. (1997). 
Brain-derived neurotrophic factor is reduced in Alzheimer's disease. 
*Molecular Brain Research*, 49(1-2), 71-81.

Phillips, H.S., Hains, J.M., Armanini, M., Laramee, G.R., Johnson, S.A., & 
Winslow, J.W. (1991). BDNF mRNA is decreased in the hippocampus of individuals 
with Alzheimer's disease. *Neuron*, 7(5), 695-702.

## Reproducing this analysis
1. Download `GSE104704_RNA-Seq_Table.txt.gz` from 
   [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704)
2. Run `scripts/load_data.R` to load and label the data
3. Run `scripts/phase3_qc.R` to reproduce quality control checks
4. Run `scripts/phase4_limma_analysis.R` to reproduce the differential expression 
   analysis, volcano plot, and heatmap
5. See `results/limma_significant_DEGs.csv` for the full significant gene list

## Contact
Maha Saleemi — [linkedin.com/in/maha-s-7bb503419](https://linkedin.com/in/maha-s-7bb503419) — mahaamjad929@gmail.com
