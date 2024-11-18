#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

family_id   <- args[1]
pheno       <- args[2]
pedigree    <- args[3]

# Load pedigree
p <- pedtools::readPed(pedigree)

# Check for mendelian errors
errors <- pedtools::mendelianCheck(p)
readr::write_lines(errors, paste(family_id, pheno, "mendelian_errors.txt", sep = "."))
