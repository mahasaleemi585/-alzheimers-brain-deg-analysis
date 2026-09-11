# Gene Expression in Alzheimer's Disease

## Research Question
Which genes are differentially expressed between aged Alzheimer's disease brains 
and age-matched healthy brains, and what do those genes suggest about the 
biological mechanisms driving neurodegeneration and potential therapeutic targets?

## Dataset
- **Source:** GEO accession [GSE104704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704) 
  (human lateral temporal lobe, RNA-seq)
- **Samples used:** 10 aged healthy controls ("Old") vs. 12 Alzheimer's disease ("AD") 
  (8 "Young" samples were excluded from this comparison)
- **Reference paper:** PubMed ID [29507413](https://pubmed.ncbi.nlm.nih.gov/29507413/)
- **AD vs. healthy definition:** AD was defined as clinically diagnosed Alzheimer's 
  disease in donor brain tissue from the lateral temporal lobe, compared against 
  younger and elderly cognitively normal controls matched by brain region.

## Repository Structure
- `scripts/` — R code (`load_data.R` loads and labels the data; `phase4_limma_analysis.R` 
  runs the differential expression analysis and generates figures)
- `figures/` — Volcano plot, heatmap, STRING network image, and candidate protein structure
- `results/` — DEG list, enrichment results, and STRING interaction/annotation tables
- `archive_first_attempt/` — An initial DESeq2-based analysis, kept for transparency 
  (see Methods below for why this was revised)

## Methods
1. Loaded processed/normalized expression data (FPKM-style values) from GEO
2. Verified sample-to-group labels using GEO sample metadata, matched by sample ID 
   (not column order, which was misaligned initially — see below)
3. **Differential expression analysis using limma.** An initial analysis used DESeq2, 
   but was corrected to limma after determining the available data was pre-normalized 
   rather than raw sequencing counts — DESeq2 assumes raw counts, and using normalized 
   data with it violates its statistical assumptions. The DESeq2-based results are 
   preserved in `archive_first_attempt/` for transparency.
4. Thresholds: adjusted p-value < 0.05, |log2 fold change| > 1
5. Enrichment analysis via g:Profiler (GO terms + KEGG pathways)
6. Protein-protein interaction network via STRING
7. Candidate target validated via UniProt (function) and RCSB PDB (structure)

## Key Findings
- **179 significant DEGs** (73 upregulated, 106 downregulated in AD)
- Enrichment highlighted immune/receptor-related processes (e.g., interleukin-8 
  receptor activity, lipoxygenase activity) alongside neuroactive signaling
- **BDNF** identified as the top hub gene in the STRING network — 14 connections, 
  more than double any other gene in the network
- Final candidate therapeutic target: **BDNF** (Brain-Derived Neurotrophic Factor) — 
  supported by its central, well-established role in synaptic plasticity, long-term 
  potentiation (LTP) and depression (LTD), and neuronal survival, plus a solved 
  structure (PDB: [1BND](https://www.rcsb.org/structure/1BND))

## Limitations
- Small sample size (10 vs. 12 samples)
- Secondary analysis of existing public data, not new wet-lab data
- No functional/experimental validation of BDNF as a therapeutic target
- The available RNA-seq data was pre-normalized rather than raw counts, requiring 
  the DESeq2-to-limma correction described above
- g:Profiler enrichment results can vary slightly between runs due to periodic 
  updates to its background annotation database; results reported here reflect 
  the analysis run on 2026-09-09

## Reproducing this analysis
1. Download `GSE104704_RNA-Seq_Table.txt.gz` from 
   [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704)
2. Run `scripts/load_data.R` in R/RStudio to load and label the data
3. Run `scripts/phase4_limma_analysis.R` to reproduce the differential expression 
   analysis, volcano plot, and heatmap
4. See `results/limma_significant_DEGs.csv` for the full significant gene list
