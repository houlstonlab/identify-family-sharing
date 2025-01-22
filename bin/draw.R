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

# phenotypes
phenotypes <- readr::read_tsv(phenotypes_file)

# variants
variants <- readr::read_tsv(variants_file)
variants <- dplyr::filter(variants, !is.na(ensembl))

# extract info
if ( type == 'complete' ) {
  variants <- dplyr::filter(variants, alt_het > 1 & is.na(alt_hom))
} else if ( type == 'partial' ) {
  variants <- dplyr::filter(variants, n_variant > 1, n_affected > 1 & is.na(alt_hom))
}

purrr::map(
  dplyr::group_split(variants, gene, variant),
  ~{
    d <- .x
    marker <- unique(unlist(strsplit(d$variant, ',')))
    
    title <- paste(c(d$famid, d$gene, marker), collapse = '\n')
    
    labs <- phenotypes$id
    names(labs) <- ifelse(
      !is.na(phenotypes$pheno),
      paste0(phenotypes$id, ' (', phenotypes$pheno, ')'),
      phenotypes$id
    )
    
    aff <- phenotypes$id[phenotypes$aff == 2]
    carr <- phenotypes$id[phenotypes$carrier == 1]
    
    file_name <- with(
      d,
      ifelse(
        length(marker) > 1,
        paste(famid, phenotype, category, type, gene, 'png', sep = '.'),
        paste(famid, phenotype, category, type, gene, variant, 'png', sep = '.')
      )
    )    
    size <- pedtools::pedsize(pedigree)
    
    png(file_name,
        width = size/1.5, height = size + length(marker)/2,
        units = 'in', res = 300)
    
    plot(pedigree,
         aff = aff,
         carrier = carr,
         labs = labs,
         title = title,
         marker = marker,
         margins = c(0.6, 1, length(marker) + 4, 1)
    )
    
    dev.off()
  }
)
