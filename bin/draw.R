#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

family_id   <- args[1]
pheno       <- args[2]
genotype    <- args[3]
phenotype   <- args[4]
pedigree    <- args[5]
which       <- args[6]
markers     <- args[7]

# family_id   <- 'FAM_01'
# pheno       <- 'pheno'
# genotype    <- 'identify-family-sharing/test/results/family/FAM_01.pheno.genotype.tsv'
# phenotype   <- 'identify-family-sharing/test/results/family/FAM_01.pheno.phenotype.tsv'
# pedigree    <- 'identify-family-sharing/test/results/family/FAM_01.pheno.ped'
# which  <- 'gene'
# markers <- 'identify-family-sharing/test/results/sharing/FAM_01.pheno.gene.sharing.tsv'

# Load pedigree
p <- pedtools::readPed(pedigree)

phenotype <- readr::read_tsv(phenotype)
affected <- dplyr::filter(phenotype, aff == 2)
affected <- dplyr::pull(affected, id)

# Load markers
markers <- readr::read_tsv(markers)

if (which == 'variant') {
  mrks <- markers$variant
  names(mrks) <- mrks
} else if (which == 'gene') {
  mrks <- purrr::map(strsplit(markers$variant, ','), unique)
  names(mrks) <- markers$gene
}

# Plot pedigree
purrr::imap(
  mrks,
  ~{
    png(filename = paste(family_id, pheno, which, .y, 'png', sep = '.'),
        width = nrow(phenotype), height = nrow(phenotype), unit = 'in', res = 300)
    plot(p, hatched = affected, marker = .x)
    dev.off()
  }
)
