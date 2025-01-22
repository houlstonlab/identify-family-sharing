#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
famid           <- args[1]
phenotype       <- args[2]
category        <- args[3]
type            <- args[4]
ped_file        <- args[5]
phenotypes_file <- args[6]
variants_file   <- args[7]
       
# Load data
# genotypes
pedigree <- pedtools::readPed(ped_file)
pedigree$FAMID <- famid
genotypes <- tibble::as_tibble(pedigree)
genotypes <- dplyr::mutate(genotypes, famid = famid)
genotypes <- tidyr::gather(genotypes, variant, genotype, tidyr::starts_with('chr'))

# phenotypes
phenotypes <- readr::read_tsv(phenotypes_file)
phenotypes <- dplyr::mutate(phenotypes, carrier = ifelse(is.na(carrier), 0, carrier))

# annotations
variants <- readr::read_tsv(variants_file)
annotations <- dplyr::select(variants, ensembl = Gene, gene = SYMBOL, variant = SNP, consequence = Consequence)

# merge tables
d <- tibble::tibble(famid, phenotype, category)
d <- dplyr::left_join(d, genotypes)
d <- dplyr::left_join(d, phenotypes)
d <- dplyr::left_join(d, annotations)

# modify individual id and phenotype, and status 
d <- dplyr::mutate(d, individual = ifelse(!is.na(pheno), paste0(id, ' (', pheno, ')'), id))
d <- dplyr::mutate(d, 
  individual_type = dplyr::case_when(
  aff == 2 ~ 'affected',
  aff != 2 & carrier == 1 ~ 'carrier',
  aff != 2 & carrier == 0 ~ 'not_affected'
  ))

# modify genotypes
d <- dplyr::mutate(d,
  genotype2 = dplyr::case_when(
    genotype == '-/-' ~ 'missing',
    genotype == 'a/a' ~ 'alt_hom',
    genotype == 'a/b' ~ 'alt_het',
    genotype == 'b/b' ~ 'ref_hom',
  ))

# Complete sharing
if ( type == 'complete' ) {
  family_summary <- dplyr::group_by(d, famid, phenotype, category, ensembl, gene, consequence, variant, genotype2)
  family_summary <- dplyr::summarise(family_summary, n = dplyr::n())
  family_summary <- tidyr::pivot_wider(family_summary, names_from = 'genotype2', values_from = 'n')
  family_summary <- dplyr::mutate(family_summary, n_variant = 1)
  
  individual_summary <- dplyr::filter(d, genotype2 == 'alt_het')
  individual_summary <- dplyr::group_by(individual_summary, famid, phenotype, category, ensembl, gene, consequence, variant, individual_type)
  individual_summary <- dplyr::summarise(individual_summary, n = dplyr::n(), individual = paste(unique(individual), collapse = ', '))
  individual_summary <- tidyr::pivot_wider(individual_summary, names_from = 'individual_type', values_from = c('n', 'individual'))
} else if ( type == 'partial' ) {
  family_summary <- dplyr::group_by(d, famid, phenotype, category, ensembl, gene, genotype2)
  family_summary <- dplyr::summarise(
    n_variant = length(unique(variant)),
    variant = paste(unique(variant), collapse = ','),
    consequence = paste(unique(consequence), collapse = ','),
    family_summary, n = dplyr::n()
  )
  family_summary <- tidyr::pivot_wider(family_summary, names_from = 'genotype2', values_from = 'n')
  
  individual_summary <- dplyr::filter(d, genotype2 == 'alt_het')
  individual_summary <- dplyr::group_by(individual_summary, famid, phenotype, category, ensembl, gene, individual_type)
  individual_summary <- dplyr::summarise(individual_summary, n = dplyr::n(), individual = paste(unique(individual), collapse = ', '))
  individual_summary <- tidyr::pivot_wider(individual_summary, names_from = 'individual_type', values_from = c('n', 'individual'))
  }

# Merge
sharing <- dplyr::left_join(family_summary, individual_summary)

# Format output
col_names <- c(
  "famid", "phenotype", "category", "ensembl", "gene", "n_variant", "variant", "consequence",   
  "alt_het", "alt_hom", "ref_hom","missing",
  "n_affected", "n_carrier", "n_not_affected", 
  "individual_affected", "individual_carrier", "individual_not_affected"
)

output.df <- data.frame(
  matrix(
    character(),
    ncol = length(col_names),
    dimnames = list(c(), col_names)),
  stringsAsFactors = FALSE
)

sharing <- dplyr::mutate_all(sharing, as.character)
sharing <- dplyr::full_join(output.df, sharing)
sharing <- dplyr::select(sharing, col_names)

# Write output to file
readr::write_tsv(
  sharing,
  paste(famid, phenotype, category, type, 'tsv', sep = '.')
)
