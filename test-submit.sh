#!/bin/bash

#SBATCH -o test/test.out
#SBATCH -e test/test.err
#SBATCH -J test
#SBATCH -p master-worker
#SBATCH -t 120:00:00

# Setup test directory
mkdir -p test/ test/input
cd test/

# # Download test data
# URL="https://figshare.com/ndownloader/files"

# wget -c $URL/50487621 -O input/pheno.variants.vcf.gz
# wget -c $URL/50487624 -O input/pheno.variants.vcf.gz.tbi
# wget -c $URL/50427141 -O input/pheno.families.ped

# Run nextflow
module load Nextflow

# nextflow run houlstonlab/identify-family-sharing -r main \
nextflow run ../main.nf \
    --output_dir ./results/ \
    -profile local,gha \
    -resume

# usage: nextflow run [ local_dir/main.nf | git_url ]  
# These are the required arguments:
#     -r            {main,dev} to run specific branch
#     -profile      {local,cluster} to run using differens resources
#     -params-file  params.json to pass parameters to the pipeline
#     -resume       To resume the pipeline from the last checkpoint

mv .nextflow.log nextflow.log