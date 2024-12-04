#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

ped_file  <- args[1]
test_file <- args[2]
out_file  <- args[3]

# ped_file  <- 'identify-family-sharing/test/results/markers/FAM_01.pheno.Damaging.marked.ped'
# test_file <- 'identify-family-sharing/test/results/tests/FAM_01.pheno.Damaging.tdt'
# out_file  <- 'identify-family-sharing/test/results/plots/FAM_01.pheno.Damaging'

# Load files
tst <- read.table(test_file, header = TRUE)
pdg <- pedtools::readPed(ped_file)

# TODO: Get affected individuals
# Extract markers
markers <- dplyr::filter(tst, P < 1)
markers <- dplyr::pull(markers, SNP)
markers <- head(markers, 3)

# Plot pedigree
purrr::map(
  markers,
  ~{
    png(filename = paste(out_file, .x, 'png', sep = '.'),
        width = pedtools::pedsize(pdg), height = pedtools::pedsize(pdg),
        unit = 'in', res = 300)
    
    plot(
      pdg,
      marker = .x,
      title = .x
    )
    
    dev.off()
  }
)
