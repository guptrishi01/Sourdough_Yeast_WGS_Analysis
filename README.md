# Experimental Evolution of Sourdough Yeast: Genomic Variant Analysis

## Acknowledgements

I would like to thank Professor Caiti Heil and the Heil Lab at North Carolina State University
for allowing me to be part of this amazing project, as well as providing consent to post my work here. 
This was my first introduction into bioinformatics and I am extremely grateful for all of the tools and 
skills that I have gained during my time working here.

I am listed as a co-author in **Demystifying Domestication of Sourdough Yeast: Evolution, Diversity and Metabolism**,
which was submitted to **The International Conference on Yeast Genetics and Molecular Biology**.


## Overview

This project investigates genetic variants in *Saccharomyces cerevisiae* that confer 
fitness advantages in the bread dough environment through an **experimental evolution 
framework**. Using whole-genome sequencing data from sourdough and nature-isolated yeast 
populations, this project identifies and characterizes three major classes of genomic 
variation:

- **Single Nucleotide Variants (SNVs)**
- **Copy Number Variants (CNVs)**
- **Loss of Heterozygosity (LOH) Events**

Variants were identified relative to ancestral strains in both experimental populations, 
allowing us to distinguish evolved, potentially adaptive mutations from standing genetic 
variation.

---

## Biological Context

Yeast strains used in sourdough baking experience strong, consistent selective pressures — 
high sugar concentrations, acidity, and temperature fluctuation. By evolving yeast populations in this environment 
and comparing their genomes to ancestral strains and nature-isolated counterparts, we can identify mutations that 
arose specifically in response to the sourdough niche.

| Population | Samples | Ancestor |
| --- | --- | --- |
| Sourdough (SD) | 46 evolved clones | 892B |
| Nature (NAT) | 26 evolved clones | ERR1308867 |

---

## Repository Structure
```
├── scripts/
│   |── sample_fastqc.sh
│   ├── sample_fastq_arrayjob.sh
│   ├── map_filter_pereads.sh
│   ├── call_variants.sh
│   ├── calling_variants_arrayjob.sh
│   ├── genomicsDBImport.sh
│   ├── genotypeGVCF.sh
│   ├── gatherVCFs.sh
│   ├── subset_filter_variants.sh
│   ├── variantFiltration.sh
│   ├── variantAnnotation.sh
│   └── variantSnpSiftFiltering.sh
│
├── LOHet/
│   ├── all_chromosomes_combined.vcf
│   ├── Calculate_LOHet.ipynb
│   ├── Plotting_ROHet.ipynb
│   └── Output/
│      ├── heterozygosity_892B_vs_sourdough_samples/
│      ├── heterozygosity_ERR1308867_vs_nature_samples/
│      ├── ROHet_final_heatmap.png
│      ├── ROHet_NAT.png
│      └── ROHet_SD.png
│
|── Phylogeny/
│   ├── zll_chromosomes_combined.vcf
│   ├── all_chromosomes_combined.gds
│   ├── new_filtered_variants.vcf
│   ├── nat_samples.txt
│   ├── sd_samples.txt
│   ├── sample_ancestors.csv
│   ├── PCA.py
│   ├── SNPrelate.R
│   └── Plots + Output/
│      ├── Phylogenetic_Neighbor_Tree.png
│      ├── Circular_Phylogenetic_Tree.png
│      ├── Annotated_Phylogenetic_Tree.png
│      ├── annotated_tree_legend.png
│      └── PCA_plot.png
│
|── CNV/
|   ├── newref.fasta
|   ├── newref.fai
|   ├── sd_samples.txt
|   ├── nat_samples.txt
|   ├── all_chromosomes_combined.vcf
|   ├── CNV_Chromosomal_Evaluation.ipynb
|   ├── CNV_Heatmap_Plotting.ipynb
|   ├── Ancestor_Descendant_CNV_Comparison.ipynb
|   ├── CNV_Gain_Loss/
|      ├── Coverage_nature_sourdough_aligned.png
|      ├── whole_genome_coverage_892B_vs_clones.pdf
|      ├── whole_genome_coverage_ERR1308867_vs_clones.pdf
|      ├── Coverage Files/
│         └── [sample]_coverage.txt / _coverage_values.txt
|      |── Output/
|         ├── [sample].BED files
|         ├── combined_nat_gain.xlsx
|         ├── combined_nat_loss.xlsx
|         ├── combined_sd_gain.xlsx
|         └── combined_sd_loss.xlsx
```

---

## Pipeline Overview

Raw paired-end FASTQ files were processed through a multi-step bioinformatics pipeline 
to produce a final, annotated variant callset. The following flowchart describes the full bioinformatics pipeline used to process 
raw sequencing reads into a filtered, annotated variant callset.

```mermaid
%%{init: {'flowchart': {'nodeSpacing': 80, 'rankSpacing': 100, 'wrappingWidth': 900}, 'themeVariables': {'fontSize': '18px', 'fontFamily': 'arial'}}}%%
flowchart TD

    A["📂 Raw Input
    ───────────────────────────────
    Paired-end FASTQ files  _R1.fastq.gz + _R2.fastq.gz
    46 Sourdough + 26 Nature + Commercial Samples"]

    B["🔍 Step 1 — Quality Control
    ───────────────────────────────
    Tool: FastQC
    Scripts: sample_fastqc.sh + sample_fastq_arrayjob.sh
    Output: Per-sample QC reports — HTML + zip — for all runs"]

    C["🗺️ Step 2 — Read Mapping
    ───────────────────────────────
    Tool: BWA-MEM
    Script: map_filter_pereads.sh
    Map paired-end reads to S. cerevisiae reference genome — newref.fasta
    Output: SAM file per sample"]

    D["⚙️ Step 3 — BAM Processing
    ───────────────────────────────
    Tools: SAMtools + Picard     Script: map_filter_pereads.sh
    1. SAMtools view — SAM to BAM
    2. SAMtools sort + index
    3. SAMtools flagstat — alignment QC statistics
    4. Picard MarkDuplicates — flag PCR duplicates
    5. Picard AddOrReplaceReadGroups
    6. SAMtools view — filter BAM — mapQ 20 — remove unmapped, secondary, QC-fail
    7. SAMtools sort + index
    Output: sample_filtered.sorted.bam per sample"]

    E["🧬 Step 4 — Per-Sample Variant Calling
    ───────────────────────────────
    Tool: GATK HaplotypeCaller
    Scripts: call_variants.sh + calling_variants_arrayjob.sh
    Run in GVCF mode — ERC GVCF — on each filtered BAM
    Output: sample.g.vcf.gz — one GVCF per sample"]

    F["🗄️ Step 5 — Consolidate GVCFs
    ───────────────────────────────
    Tool: GATK GenomicsDBImport     Script: genomicsDBImport.sh
    Create TSV sample map — consolidate all GVCFs per chromosome
    into a shared GenomicsDB workspace
    Output: CHROM_gdb workspace — one per chromosome"]

    G["🔗 Step 6 — Joint Genotyping
    ───────────────────────────────
    Tool: GATK GenotypeGVCFs     Script: genotypeGVCF.sh
    Jointly genotype all samples from GenomicsDB workspace
    Run per chromosome — leverages cohort-wide information
    Output: CHROM_genotyped_variants.vcf — one per chromosome"]

    H["📦 Step 7 — Merge Chromosomal VCFs
    ───────────────────────────────
    Tool: Picard GatherVcfs     Script: gatherVCFs.sh
    Combine all 16 per-chromosome genotyped VCFs into one cohort-wide callset
    ⭐ Output: all_chromosomes_combined.vcf"]

    I["🚦 Step 8 — Hard Filter Variants
    ───────────────────────────────
    Tool: GATK VariantFiltration     Script: variantFiltration.sh
    QD less than 2.0 — FS greater than 60.0 — SOR greater than 3.0
    MQ less than 50.0 — MQRankSum less than -5.0 — ReadPosRankSum less than -4.0
    grep PASS to retain only passing variants
    Output: new_filtered_variants.vcf"]

    J["🧹 Step 9 — Filter Ancestral Variants
    ───────────────────────────────
    Tools: bcftools + htslib     Script: subset_filter_variants.sh
    bgzip + index VCF
    bcftools view — split into SD ancestor, NAT ancestor, SD descendants, NAT descendants
    bcftools isec — remove variants shared between ancestor and descendants
    Output: Variants private to SD evolved clones + NAT evolved clones"]

    K["🏷️ Step 10 — Variant Annotation
    ───────────────────────────────
    Tool: SnpEff     Script: variantAnnotation.sh
    Annotate functional effects — missense, synonymous, stop gained, frameshift, etc.
    Output: Functionally annotated VCF — SD private variants + NAT private variants"]

    L["🔎 Step 11 — Annotation Filtering
    ───────────────────────────────
    Tool: SnpSift     Script: variantSnpSiftFiltering.sh
    Filter annotated variants by functional impact — HIGH and MODERATE priority
    Output: Final filtered and annotated variant sets of interest"]

    M1["📊 LOH Analysis
    ───────────────────────────────
    Input: all_chromosomes_combined.vcf
    Notebooks: Calculate_LOHet.ipynb + Plotting_ROHet.ipynb
    Output: LOH heatmaps + ancestor vs descendant heterozygosity plots"]

    M2["🌳 Phylogeny and PCA
    ───────────────────────────────
    Input: all_chromosomes_combined.vcf + new_filtered_variants.vcf
    Scripts: SNPRelate.R + PCA.py
    Output: Neighbor-joining tree, circular tree, annotated tree, PCA plot"]

    M3["📈 CNV Analysis
    ───────────────────────────────
    Input: Per-sample coverage files + all_chromosomes_combined.vcf
    Notebooks: CNV_Chromosomal_Evaluation + CNV_Heatmap_Plotting + Ancestor_Descendant_CNV
    Output: CNV heatmaps, BED files, gain and loss tables, coverage plots"]

    A --> B
    B --> C
    C --> D
    D --> E
    E --> F
    F --> G
    G --> H
    H --> I
    I --> J
    J --> K
    K --> L
    K --> M1
    K --> M2
    K --> M3

    classDef default fill:#2E86AB,color:#fff,stroke:#1a5c7a,padding:20px
    classDef highlight fill:#27AE60,color:#fff,stroke:#1a7a42,padding:20px
    classDef downstream fill:#8E44AD,color:#fff,stroke:#5d2d73,padding:20px
    classDef annotation fill:#16A085,color:#fff,stroke:#0e6b59,padding:20px

    style A fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style B fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style C fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style D fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style E fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style F fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style G fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style H fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style I fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style J fill:#2E86AB,color:#fff,stroke:#1a5c7a
    style K fill:#16A085,color:#fff,stroke:#0e6b59
    style L fill:#E67E22,color:#fff,stroke:#a85515
    style M1 fill:#8E44AD,color:#fff,stroke:#5d2d73
    style M2 fill:#8E44AD,color:#fff,stroke:#5d2d73
    style M3 fill:#8E44AD,color:#fff,stroke:#5d2d73
```

### Tools Used

| Step | Tool |
|---|---|
| Quality Control | FastQC |
| Read Mapping | BWA-MEM |
| SAM/BAM Processing | SAMtools, Picard |
| Variant Calling | GATK HaplotypeCaller |
| Joint Genotyping | GATK GenomicsDBImport, GenotypeGVCFs |
| VCF Merging | Picard GatherVcfs |
| Variant Filtering | GATK VariantFiltration, bcftools |
| Variant Annotation | SnpEff |
| Variant Subsetting | SnpSift, bcftools isec |
| CNV Analysis | Custom Python (Jupyter Notebooks) |
| LOH Analysis | Custom Python (Jupyter Notebooks) |
| Phylogenetics | SNPRelate (R), Custom Python PCA |

---

---

## Analysis Modules

### 🔬 Variant Calling Pipeline (`scripts/`)

Shell scripts implementing the full GATK best-practices pipeline from raw FASTQ files 
to a filtered, annotated multi-sample VCF (`all_chromosomes_combined.vcf`).

**Key steps:**
1. FastQC quality assessment of raw reads
2. BWA-MEM mapping to *S. cerevisiae* reference genome
3. SAMtools/Picard BAM processing and duplicate marking
4. GATK HaplotypeCaller in GVCF mode (per-sample)
5. GenomicsDBImport consolidation (per-chromosome)
6. Joint genotyping with GenotypeGVCFs
7. Picard GatherVcfs merging into `all_chromosomes_combined.vcf`
8. GATK VariantFiltration (hard filtering, GATK best practices)
9. bcftools filtering to remove ancestral variants
10. SnpEff annotation and SnpSift filtering

---

---

### Loss of Heterozygosity (`LOHet/`)

**`Calculate_LOHet.ipynb`**
Detects LOH events genome-wide using variant call data from `all_chromosomes_combined.vcf`. 
The pipeline:
- Identifies heterozygous sites in the ancestral strain
- Tracks those sites in evolved clones
- Flags regions with 3+ consecutive homozygous sites (windowed in 500 bp bins)
- Characterizes LOH events by size, location, and gene content

**`Plotting_ROHet.ipynb`**
Generates chromosome-wide heterozygosity profiles for all samples and comparative 
visualizations between SD and NAT populations.

**Key Outputs:**
- Per-sample ancestor vs. descendant heterozygosity plots (SD and NAT)
- Full genome-wide LOH heatmap across all 16 chromosomes (10,000 bp windows, chrM excluded)
- Separate SD and NAT heatmaps

---

### Phylogenetics & PCA (`Phylogeny/`)

**`SNPrelate.R`**
- Converts `all_chromosomes_combined.vcf` to GDS format
- Computes an Identity-By-State (IBS) distance matrix
- Builds a neighbor-joining phylogenetic tree

**`PCA.py`**
- Reads `new_filtered_variants.vcf`
- Converts genotypes to alternate allele counts
- Imputes missing values and assigns sample groups (SD, NAT, commercial)
- Generates PCA plot separating populations

**Key Outputs:**
- Neighbor-joining phylogenetic tree (full layout and circular)
- Annotated phylogenetic tree with color-coded sample groups
- PCA plot of SD, NAT, and commercial samples

---

### Copy Number Variation (`CNV/`)

**`CNV_Chromosomal_Evaluation.ipynb`**
Analyzes chromosome-level CNV to detect whole-chromosome duplications or deletions 
using normalized read depth.

**`CNV_Heatmap_Plotting.ipynb`**
Visualizes genome-wide CNV patterns across all samples as heatmaps.

**`Ancestor_Descendant_CNV_Comparison.ipynb`**
Directly compares CNV profiles between ancestral and evolved strains to identify 
patterns of genomic gain and loss over the course of experimental evolution.

**Key Outputs:**
- CNV gain/loss BED files per sample
- Compiled gain/loss tables for SD and NAT populations (`.xlsx`)
- Genome-wide coverage comparison plots
- CNV heatmaps

---
