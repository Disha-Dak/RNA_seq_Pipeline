# RNA_seq_Pipeline



├── featureCounts/          ← “Step 1”: raw FASTQ → featureCounts
│   ├── rna_seq.sh          ← Your Bash pipeline
│   ├── data/               ← (optional) example FASTQ or links
│ 
│
└── deseq2/                 ← “Step 2”: featureCounts → DESeq2 analysis
    ├── _targets.R          ← targets pipeline
    ├── results/            ← outputs (PCA, boxplots, etc.)

    This repository contains two modular pipelines:

1. **featureCounts** (Bash):  
   - Downloads SRA FASTQ  
   - QC (FastQC) → Trim (Trimmomatic) → Align (HISAT2) → Quantify (featureCounts)  
2. **DESeq2** (R + targets):  
   - Imports raw counts  
   - Applies two normalization methods (median-of-ratios & TMM)  
   - Filters zero-variance genes  
   - Generates PCA & boxplots via a reusable, scalable `targets` script

## Quick Start

```bash
# Clone repo
git clone https://github.com/Disha-Dak/RNA_seq_Pipeline.git


# Step 1: featureCounts
cd featureCounts
bash rna_seq.sh
Fetch SRA runinfo → download FASTQ via xargs wget
Raw QC (FastQC)
Trim adapters (Trimmomatic)
Post-trim QC (FastQC)
Align to GRCh38 (HISAT2)
Count features (featureCounts)

# Step 2: DESeq2 analysis
cd ../deseq2
Rscript -e "targets::tar_make()"
Reads data/GSE60424_raw_counts_…tsv
Cleans NA genes, filters by Sample Name & celltype
Builds DESeqDataSet, runs DESeq2
Two normalizations:DESeq2 median-of-ratios,edgeR TMM
Filters zero-variance genes for PCA
Outputs PCA plots (results/PCA_*.png) and combined boxplot
   
