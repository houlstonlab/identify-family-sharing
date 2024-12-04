#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

bim       <- args[1]
bed       <- args[2]
fam       <- args[3]
ped_file  <- args[4]
out_file  <- args[5]

# bim       <- 'identify-family-sharing/test/results/plinked/FAM_01.pheno.Damaging.bim'
# bed       <- 'identify-family-sharing/test/results/plinked/FAM_01.pheno.Damaging.bed'
# fam       <- 'identify-family-sharing/test/results/plinked/FAM_01.pheno.Damaging.fam'
# ped_file  <- 'identify-family-sharing/test/results/pedigree/FAM_01.ped'
# out_file  <- 'identify-family-sharing/test/results/plots/FAM_01.pheno.Damaging'

# Load files
pdg <- pedtools::readPed(ped_file, colSkip = 5)
plnk <- snpStats::read.plink(bed, bim, fam)

convert_genotypes <- function(s) {
  dplyr::case_when(
    s == '00' ~ '-/-',
    s == '01' ~ 'a/a',
    s == '02' ~ 'a/b',
    s == '03' ~ 'b/b'
  )
}

create_markers <- function(ped, genotypes) {
  markers <- purrr::map(
    1:ncol(genotypes),
    ~{
      # Get genotypes
      g <- unclass(genotypes[, .x])
      s <- as.character(g)
      s <- convert_genotypes(s)
      names(s) <- rownames(g)
      
      # Get variant information
      name <- colnames(g)
      var <- unlist(strsplit(name, ':'))
      
      # Create markers
      pedtools::marker(
        ped,
        geno = s,
        chrom = as.character(var[1]),
        posMb = as.integer(var[2])/1000000,
        name = name
      )
    }
  )
  names(markers) <- colnames(genotypes)
  
  return(markers)
}

# Convert genotypes to markers, and attach to pedigree
mms <- create_markers(pdg, plnk$genotypes)
pdg2 <- pedtools::setMarkers(pdg, mms)

# Write pedigree
pedtools::writePed(pdg2, out_file)
