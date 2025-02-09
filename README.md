### Introduction

This workflow tests filters and tests variants shared in families.

### Usage

The typical command looks like the following. `--cohorts` is a required input. 
Different versions of the workflow can be called using `-r` and output directed
to `--output_dir`.

```bash
nextflow run houlstonlab/identify-family-sharing \
    -r main \
    --output_dir results/ \
    --cohorts input/cohorts_input.csv
```

### Inputs & Parameters

- `genome`    : genome version (default is'hg38') 
- `style`     : chromosome names style 'UCSC' or 'NCBI'
- `categories`: selection categories. One or more of `'Pathogenic,Damaging,Splicing,High,PTV,Stop'`
- `GQ`        : minimum genotyping quality. Default `> 10`
- `DP`        : minimum read depth. Default `> 5`
- `VAF`       : variant allele frequency. Default `> 0.2`
- `AF`        : maximum allele frequence. Default `> 0.5` 
- `AF_COL`    : allele frequency column
- `MAX_AF`    : maximum group allele frequence. Default `> 0.5`
- `DS`        : minimum delta score for spliceAI
- `CADD`      : minimum CADD score. Defalt `> 5`
- `type`      : one or more of `'complete'` and `'partial'` sharing

### Output

- `pedigree` : pedigrees of individual families
- `subsets`  : subsetted VCF files
- `filtered` : filtered variants in VCF
- `extracted`: variant annotations
- `plinked`  : variants in plink format
- `markers`  : pedigrees with attached genotypes
- `sharing`  : a table of variant sharing among family members
- `reports`  : reports of relationships and mendelian errors
- `summary`  : combined sharing files
- `plots`    : Plots of shared variants
