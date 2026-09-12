# Gene Expression Analysis in Alzheimer's Disease

## Research Question
Which genes are differentially expressed between aged Alzheimer's disease brains 
and age-matched healthy brains, and what do those genes suggest about the 
biological mechanisms driving neurodegeneration and potential therapeutic targets?

## Dataset
- **Source:** GEO accession [GSE104704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704) 
  (human lateral temporal lobe, RNA-seq)
- **Samples used:** 10 aged healthy controls ("Old") vs. 12 Alzheimer's disease ("AD") 
  (8 "Young" samples were excluded from this comparison to isolate age-matched, 
  disease-specific differences)
- **Reference paper:** PubMed ID [29507413](https://pubmed.ncbi.nlm.nih.gov/29507413/)
- **AD vs. healthy definition:** AD was defined as clinically diagnosed Alzheimer's 
  disease in donor brain tissue from the lateral temporal lobe, compared against 
  younger and elderly cognitively normal controls matched by brain region.

## Repository Structure
- `scripts/` — R code:
  - `load_data.R` — loads the expression matrix and assigns sample groups
  - `phase4_limma_analysis.R` — differential expression analysis (limma), volcano plot, heatmap
  - `phase3_qc.R` — quality control: PCA, sample correlation heatmap, outlier sensitivity check
- `figures/` — Volcano plot, DEG heatmap, PCA plot, correlation heatmap, STRING network image, candidate protein structure
- `results/` — DEG list, enrichment results (g:Profiler), and STRING interaction/annotation tables
- `archive_first_attempt/` — An initial DESeq2-based analysis, preserved for transparency 
  (see Methodological Note below)

## Methods

**Dataset.** RNA-seq gene expression data was obtained from NCBI GEO accession GSE104704. 
Sample-to-group labels were verified against GEO sample metadata by matching individual 
sample IDs directly, rather than assuming column order (an initial mismatch was caught 
and corrected during this step).

**Data type confirmation.** The GEO series record for GSE104704 explicitly states that 
raw sequencing data are available separately via SRA, while the deposited supplementary 
file represents processed expression data. Consistent with this, expression values in 
the supplementary file were continuous, non-integer numbers (e.g., FPKM-style values), 
confirming this dataset required a normalized-data-appropriate method rather than a 
raw-count method.

**Quality control.** Expression data was filtered to remove low-expression genes 
(mean log2 expression ≤ 1) prior to analysis. Principal Component Analysis (PCA) and 
a sample correlation heatmap were used to assess overall sample structure and detect 
outliers. One sample (X11.10T.RNA, Old group) was identified as an outlier by both 
methods independently. A sensitivity analysis was performed by re-running the 
differential expression analysis without this sample: 248 genes reached significance 
without it (vs. 179 with it included), and 158 of the original 179 genes (88%) 
remained significant, indicating the main findings are not driven by this single 
sample. The outlier was retained in the primary analysis to preserve sample size.

**Differential expression analysis.** Differential expression was analyzed using 
limma, appropriate for the normalized, continuous expression data available for this 
dataset. Expression values were log2-transformed prior to analysis. Genes were 
considered significantly differentially expressed at an adjusted p-value 
(Benjamini-Hochberg) below 0.05 and an absolute log2 fold change greater than 1.

**Functional enrichment.** The significant gene list was submitted to g:Profiler 
(biit.cs.ut.ee/gprofiler) to identify enriched GO terms and KEGG pathways (Homo sapiens).

**Protein interaction network.** The same gene list was submitted to STRING 
(string-db.org) to construct a protein-protein interaction network. Hub genes were 
identified by connectivity (node degree) within this network; high connectivity was 
treated as one line of supporting evidence for further investigation, not as proof of 
disease relevance on its own.

**Candidate validation.** The top hub gene was investigated using UniProt (function) 
and RCSB PDB (structure).

## Methodological Note: DESeq2 to limma Correction
An initial differential expression analysis was performed using DESeq2. This was 
later identified as methodologically incorrect: DESeq2 requires raw integer sequencing 
counts, but the GSE104704 supplementary file contains pre-normalized, continuous 
expression values (confirmed via non-integer decimal values and the GEO series 
description, which notes raw counts are hosted separately via SRA). Using normalized 
data with DESeq2 violates its statistical model and produces unreliable significance 
estimates.

This was corrected by switching to limma, which is designed for normalized, 
continuous expression data. The original DESeq2-based analysis is preserved in 
`archive_first_attempt/` for transparency, alongside the corrected limma-based 
analysis in `scripts/` and `results/`.

## Key Findings
- **179 significant DEGs** (73 upregulated, 106 downregulated in AD)
- Most significant genes included SGO1, KLF15 (up) and RPH3A, VGF, CRH, NEUROD6 (down)
- Enrichment highlighted immune/receptor-related processes (e.g., interleukin-8 
  receptor activity, lipoxygenase activity) alongside neuronal signaling genes
- **BDNF** identified as the top hub gene in the STRING network (14 connections, 
  more than double any other gene), emerging as a **candidate gene warranting further 
  investigation** — supported independently by existing literature on reduced BDNF 
  in AD brain tissue (see References)

## Limitations
- **Small sample size:** 10 Old vs. 12 AD samples is modest for RNA-seq differential 
  expression analysis; findings should be considered preliminary.
- **Secondary analysis of public data:** no new wet-lab samples were generated.
- **Pre-normalized data:** required the DESeq2-to-limma correction described above.
- **No experimental validation:** BDNF was identified computationally based on network 
  connectivity and known literature. No functional, cellular, or animal-model 
  experiments were conducted as part of this project, and no drug interaction was tested.
- **Association, not causation:** network connectivity and differential expression 
  indicate association with AD status, not a demonstrated causal or therapeutic role.
- **Bulk tissue analysis:** expression differences reflect an average across multiple 
  cell types (neurons, glia, immune cells); some signal may reflect shifts in cell-type 
  composition rather than transcriptional change within a single cell type.
- **Tool result variability:** g:Profiler enrichment results can vary modestly between 
  runs due to periodic background database updates; results reported here reflect the 
  analysis run on 2026-09-09.
- **Outlier sample:** one sample (X11.10T.RNA) was flagged as an outlier via PCA and 
  correlation analysis but retained in the main analysis (see Quality Control above).

## Software and Package Versions
- R version 4.6.1 (2026-06-24)
- Platform: Windows 11 x64
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
