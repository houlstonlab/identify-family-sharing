#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

family_id   <- args[1]
pedigree    <- args[2]
pheno       <- args[3]
category    <- unlist(strsplit(args[4], ','))
genotypes   <- unlist(strsplit(args[5], ','))

# family_id   <- 'FAM_01'
# pedigree    <- 'identify-family-sharing/test/results/family/FAM_01.ped'
# pheno       <- 'Splicing'
# genotypes   <- 'identify-family-sharing/test/results/variants/FAM_01.pheno.Splicing.genotypes.tsv'

# Load pedigree
family_ped <- readr::read_tsv(pedigree)
family_ped <- dplyr::filter(family_ped, famid == family_id)

# Load phenotype
phenotype <- dplyr::select(family_ped, id, aff)
readr::write_tsv(phenotype, paste(family_id, pheno, 'phenotype', 'tsv', sep = '.'))

# Load genotypes
# Rename variants
d <- purrr::map(genotypes, readr::read_tsv, col_names = FALSE)
d <- purrr::set_names(d, category)
d <- dplyr::bind_rows(d, .id = 'category')
d <- dplyr::mutate(d, name = paste(category, dplyr::row_number(), sep = '.'))
d <- tidyr::gather(d, X, genotype, -X1, -X2, -category, -name)
d <- dplyr::mutate(d, gene = X1, variant = paste(category, X2, sep = '.'))
d <- dplyr::select(d, gene = X1, variant = X2, category, name, genotype)
d <- tidyr::separate(d, genotype, c('id', 'genotype'), sep = '=')
d <- tidyr::separate(d, variant, c('chrom', 'pos', 'ref', 'alt'), sep = ':', remove = FALSE) 
readr::write_tsv(d, paste(family_id, pheno, 'genotype', 'tsv', sep = '.'))

# Create pedigree
p <- with(family_ped, pedtools::ped(id, fid, mid, sex))

# Select and format genotypes
mrks <- dplyr::group_split(d, name)
mrks <- purrr::map(mrks,
    ~{
    g <- d$genotype[match(pedtools::untypedMembers(p), .x$id)]
    g <- stringr::str_replace_all(g, '0', 'a')
    g <- stringr::str_replace_all(g, '1', 'b')
    
    pedtools::marker(p,
        geno = g,
        chrom = unique(.x$chrom),
        posMb = as.numeric(unique(.x$pos))/1000000,
        # alleles = c(unique(.x$ref), unique(.x$alt)),
        NAstrings = NA,
        name = unique(.x$name)
    )
  })
if (length(mrks) > 0) {
  # Add genotypes to pedigree
  p2 <- pedtools::setMarkers(p, mrks)
} else {
  p2 <- p
}

pedtools::writePed(p2, paste(family_id, pheno, sep = '.'))
