#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

family_id   <- args[1]
pheno       <- args[2]
genotype    <- args[3]
phenotype   <- args[4]
pedigree    <- args[5]
which       <- args[6]

# family_id   <- 'FAM_01'
# pheno       <- 'pheno'
# genotype    <- 'identify-family-sharing/test/results/family/FAM_01.pheno.genotype.tsv'
# phenotype   <- 'identify-family-sharing/test/results/family/FAM_01.pheno.phenotype.tsv'
# pedigree    <- 'identify-family-sharing/test/results/family/FAM_01.pheno.ped'
# which  <- 'gene'

# Load phenotype
phenotype <- readr::read_tsv(phenotype)
affected <- dplyr::filter(phenotype, aff == 2)
affected <- dplyr::pull(affected, id)

# Load pedigree
p <- pedtools::readPed(pedigree)
p2 <- pedtools::as_kinship2_pedigree(p)
p2$affected <- phenotype$aff[match(phenotype$id, p2$id)] - 1

# Sharing Probability for One Family, One Variant
prob <- RVS::RVsharing(p2, useAffected = TRUE)

genotype <- readr::read_tsv(genotype)
genotype <- dplyr::select(genotype, gene, variant = name)

g <- pedtools::getGenotypes(p)

m <- dplyr::as_tibble(g)
m <- dplyr::mutate(m, id = rownames(g))
m <- tidyr::gather(m, variant, genotype, -id)
m <- dplyr::mutate(m, pheno = ifelse(id %in% affected, 'affected', 'not_affected'))
m <- dplyr::right_join(genotype, m)

if ( which == 'variant' ) {
  m <- dplyr::group_by(m, gene, variant, pheno, genotype)
  m <- dplyr::summarise(m, n = length(id), id = paste(id, collapse = ','))
} else if ( which == 'gene' ) {
  # TODO: generate similar table for genes
  m <- dplyr::group_by(m, gene, pheno, genotype)
  m <- dplyr::summarise(m, variant = paste(variant, collapse = ','), n = length(id), id = paste(id, collapse = ','))
}

m <- dplyr::ungroup(m)
m <- tidyr::pivot_wider(m, names_from = 'pheno', values_from = c('n', 'id'))
m <- dplyr::mutate_at(m, dplyr::vars(dplyr::starts_with('n_')), ~ifelse(is.na(.x), 0, .x))

# TODO: add probability of sharing
# m <- dplyr::mutate(m, famid = family_id, prob = prob)
m <- dplyr::mutate(m, famid = family_id)
m <- dplyr::select(m, famid, dplyr::everything())

sharing <- dplyr::filter(m, genotype != 'a/a')

# TODO: filter by n_affected and return all
# sharing <- dplyr::filter(m, genotype != 'a/a', n_affected > 1)
sharing <- head(sharing, 10)

readr::write_tsv(sharing, paste(family_id, pheno, which, 'sharing', 'tsv', sep = "."))
