### Introduction

This workflow tests filters and tests variants shared in families.

### Usage

The typical command looks like the following. `--variants`, `--pedigrees` and 
`--phenotypes` are required inputs. Different versions of the workflow can be 
called using `-r` and output directed to `--output_dir`

```bash
nextflow run houlstonlab/identify-family-sharing \
    -r main \
    --output_dir results/ \
    --variants *.variants.vcf.gz{,.tbi} \
    --pedigrees *.families.ped \
    --phenotypes *.phenotypes.tsv
```

### Inputs & Parameters

- `variants`  : a VCF file with variants
- `pedigrees` : a ped file with family pedigrees
- `phenotypes`: a tsv file with the following columns: famid, id, aff, carrier 
and pheno
- `genome`: genome version (default is'hg38') 
- `style`: chromosome names style 'UCSC' or 'NCBI'
- `categories`: one or more filtering categories 'Damaging, HIGH, Rare, etc'
    
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
