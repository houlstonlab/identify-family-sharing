#' Read pedigree from file
#'
#' @param ped_file A character. Path to pedigree file
#' @param famid A logical (TRUE) to add family id
#'
#' @return A pedigree object
#'
#' @examples
#' # read a separate pedigree file
#' fl <- 'workflow/lib/extdata/FACT0246.pedigree.tsv'
#' readPed(fl)
#' 
#' # read a combined pedigree file
#' fl <- 'workflow/lib/extdata/combined.pedigrees.tsv'
#' readPed(fl)
#' 
#' @importFrom reader read_tsv
#' @importFrom kinship2 pedigree
#' 
#' @export
readPed <- function(ped_file, famid = TRUE) {
  # identify and rename columns
  col.names <- c('famid', 'id', 'sex', 'fid', 'mid', 'status', 'affected')
  
  # read family data
  fam <- readr::read_tsv(ped_file, col_names = col.names)
  # fam$affected <- ifelse(fam$affected == 0, NA, fam$affected - 1)
  
  # create pedigree object
  fams <- split(fam, fam$famid)
  fam.ped <- purrr::map(
    fams, 
    ~{
      f <- with(.x, kinship2::pedigree(id, fid, mid, sex, affected))
      # add family id when true
      f$famid <- unique(with(.x, famid))
      f
    }
  )
  if (length(fam.ped) == 1) {
    fam.ped <- fam.ped[[1]]
  }
  
  # return pedigree object
  return(fam.ped)
}

#' Read pedigree and genotype (snp matrix) similar to PLINK
#'
#' @param ped_file A character. Path to file in pedigree format
#' @param genotype_file A character. Path to genotype in vcf format
#'
#' @return A list
#' 
#' @examples
#' # read a combined pedigree and genotype 
#' ped_file <- 'workflow/lib/extdata/combined.pedigrees.tsv'
#' genotype_file <- 'workflow/lib/extdata/combined.variants.vcf.gz'
#' 
#' readPedGT(ped_file, genotype_file)
#' 
#' @importFrom reader read_tsv
#' @importFrom VariantAnnotation readVcf
#' @importFrom VariantAnnotation genotypeToSnpMatrix
#' 
#' @export
readPedGT <- function(ped_file, genotype_file) {
  # create pedigree object
  ped <- readPed(ped_file)
  
  # create famInfo object
  cols <- c('pedigree', 'member', 'father', 'mother', 'sex', 'affected')
  fam <- readr::read_tsv(ped_file, col_names = cols)
  
  # create a snp matrix
  vcf <- VariantAnnotation::readVcf(genotype_file)
  gt <- VariantAnnotation::genotypeToSnpMatrix(vcf[, colnames(vcf) %in% fam$member])
  
  # add ped and fam to the snp matrix
  gt$fam <- fam[fam$member %in% colnames(vcf),]
  gt$ped <- ped
  
  return(gt)
}

#' Create combinations of carriers
#'
#' @param ped A pedigree object
#'
#' @return A list of sets of carriers. Items are carrier.sets, carrier.sets.len,
#' and carrier.sets.prob.
#'
#' @examples
#' # read pedigree file
#' fl <- 'workflow/lib/extdata/FACT0246.pedigree.tsv'
#' ped <- readPed(fl)
#' 
#' # create a carrier set
#' carrierSet(ped)
#' 
#' @export
carrierSet <- function(ped) {
  # extract carriers/affected
  carriers <- with(ped, id[affected == 1 & !is.na(affected)])
  
  # get combinations
  carrier.sets <- list()
  for (i in length(carriers):1) {
    carrier.sets <- c(carrier.sets, combn(carriers, i, simplify=FALSE))
  }
  
  # get length of each set
  carrier.sets.len <- list()
  for (i in 1:length(carrier.sets)) {
    carrier.sets.len <- unlist(c(carrier.sets.len, length(carrier.sets[[i]])))
  }
  
  # calculate probability of sharing in the carrier sets
  carrier.sets.prob <- list()
  for (i in 1:length(carrier.sets)) {
    carrier.sets.prob <- unlist(c(carrier.sets.prob, RVS::RVsharing(ped, carriers = carrier.sets[[i]])))
  }
  
  res <- list(
    carrier.sets = carrier.sets,
    carrier.sets.len = carrier.sets.len,
    carrier.sets.prob = carrier.sets.prob
  )
  
  return(res)
}

#' Get phenotype of pedigree members
#'
#' @param ped A pedigree object
#' @param affected A logical (TRUE) get ids of affected members
#'
#' @return A character vector
#' 
#' @examples
#' # example code
#' 
#' # read pedigree file
#' fl <- 'workflow/lib/extdata/FACT0246.pedigree.tsv'
#' ped <- readPed(fl)
#' 
#' # affected members
#' getPhenotype(ped)
#' 
#' # non-affected members
#' getPhenotype(ped, FALSE)
#' 
#' @export
getPhenotype <- function(ped, affected = TRUE) {
  # extract ids and affected status
  id <- ped$id
  aff <- ped$affected
  
  if ( affected ) {
    # return affects
    res <- id[aff == 1 & !is.na(aff)]
  } else {
    # return non affected
    res <- id[aff == 0 & !is.na(aff)]
  }
  
  return(res)
}

#' Subset genotype matrix to pedigree members
#'
#' @param ped A pedigree object
#' @param genotypes A genotype matrix
#'
#' @return A genotype matrix
#'
#' @examples
#' # read a combined pedigree and genotype 
#' ped_file <- 'workflow/lib/extdata/combined.pedigrees.tsv'
#' genotype_file <- 'workflow/lib/extdata/combined.variants.vcf.gz'
#' 
#' sm <- readPedGT(ped_file, genotype_file)
#' subsetGT(sm$ped[[1]], sm$genotypes)
#' 
#' @export
subsetGT <- function(ped, genotypes) {
  # subset to pedigree members
  ind <- intersect(ped$id, rownames(genotypes))
  gt <- genotypes[ind,]
  
  # remove if no variance in affected subjects
  aff <- getPhenotype(ped)
  ind <- apply(gt[aff,], 2, function(x) any(x == '02'))
  gt <- gt[, ind]
  
  return(gt)
}

#' Get sharing pattern from a pedigree and a genotype matrix
#'
#' @param ped A pedigree object
#' @param genotypes A genotype matrix
#' @param affected A logical (TRUE) get ids of affected members
#' @param complete A logical (TRUE) return complete sharing
#'
#' @return A logical vector for complete and a character vector of sample names otherwise
#'
#' @examples
#' # read a combined pedigree and genotype 
#' ped_file <- 'workflow/lib/extdata/combined.pedigrees.tsv'
#' genotype_file <- 'workflow/lib/extdata/combined.variants.vcf.gz'
#' 
#' sm <- readPedGT(ped_file, genotype_file)
#' getSharing(sm$ped[[1]], sm$genotypes, TRUE, FALSE)
#' 
#' @export
getSharing <- function(ped, genotypes, affected = TRUE, complete = TRUE) {
  # get a subset of subjects
  aff <- getPhenotype(ped, affected)
  
  if ( complete ) {
    # return logical when complete
    res <- apply(genotypes[aff,], 2, function(x) all(x == '02'))
  } else {
    # return names otherwise
    res <- apply(genotypes[aff,], 2, function(x) names(x)[x == '02'])
  }
  
  return(res)
}

#' Evaluate sharing probability form a pedigree and a genotype matrix
#'
#' @param ped A pedigree object
#' @param faminfo A data.frame of pedigree subjects
#' @param genotypes A genotype matrix
#' @param genes A character vector of gene names
#' @param type A character vector sharing type
#'
#' @return A data.frame
#' 
#' @examples
#' # read a combined pedigree and genotype 
#' ped_file <- 'workflow/lib/extdata/combined.pedigrees.tsv'
#' genotype_file <- 'workflow/lib/extdata/combined.variants.vcf.gz'
#' 
#' sm <- readPedGT(ped_file, genotype_file)
#' 
#' # get genes annotation
#' genes <- getGenes(genotype_file)
#' 
#' evalSharing(sm$ped, sm$fam, sm$genotypes[, 1:100], genes[1:100], type = 'complete')
#' evalSharing(sm$ped, sm$fam, sm$genotypes[, 1:100], genes[1:100], type = 'partial')
#' evalSharing(sm$ped, sm$fam, sm$genotypes, genes, type = 'gene')
#' 
#' @importFrom purrr map
#' @importFrom purrr map_df
#' @importFrom purrr map2_df
#' @importFrom dplyr tibble
#' @importFrom dplyr reduce
#' @importFrom dplyr bind_rows
#' @importFrom RVS RVsharing
#' @importFrom RVS multipleVariantPValue
#' 
#' @export
evalSharing <- function(ped, faminfo, genotypes, genes, type = c('complete', 'partial', 'gene')) {
  
  if ( type == 'complete') {
    # genotypes <- genotypes[, 1:100] # TODO: to be removed
    # get sharing
    shared <- purrr::map(ped, ~getSharing(.x, genotypes))
    shared_tidy <- purrr::map_df(
      shared, ~dplyr::tibble(variant = names(.x)[.x], gene = genes[.x]),
      .id = 'family'
    )

    # calculate probability per family
    share_prob <- RVS::RVsharing(ped, useAffected = TRUE)
    share_prob_tidy <- dplyr::tibble(prob = share_prob, family = names(share_prob))
    
    # calculate probability across families
    share_prob_all <- RVS::multipleVariantPValue(genotypes, faminfo, share_prob)
    share_prob_all_tidy <- dplyr::tibble(
      variant = names(share_prob_all$pvalues),
      potential_pvalues = share_prob_all$potential_pvalues,
      across_pvalue = share_prob_all$pvalues
    )
    
    # tidy output
    res <- purrr::reduce(
      list(shared_tidy, share_prob_tidy, share_prob_all_tidy),
      dplyr::left_join
    )
  } else if (type == 'partial') {
    # genotypes <- genotypes[, 1:100] # TODO: to be removed
    
    # per family
    res <- purrr::map_df(
      ped, 
      function(x) {
        # create a carrier set
        carrier_sets <- carrierSet(x)
        probs <- carrier_sets$carrier.sets.prob
        n <- carrier_sets$carrier.sets.len
        
        # get the sharing pattern
        aff <- getPhenotype(x, affected = TRUE)
        shared <- getSharing(x, genotypes, complete = FALSE)
        
        # calculate probability of observed event
        df <- purrr::imap(
          shared,
          ~ {
            if (length(.x) > 0 & length(.x) < length(aff)) {
              # get gene name
              g <- genes[which(colnames(genotypes) == .y)]

              # calculate probability
              prob <- RVS::RVsharing(x, .x, useAffected = TRUE)
              
              # calculate pvalue
              p <- sum(probs[probs <= prob & n >= length(.x)])
              
              # tidy output
              dplyr::tibble(
                variant = .y,
                gene = g,
                prob = prob,
                pvalue = p
              )
            }
          })
        dplyr::bind_rows(df)
      },
      .id = 'family')
  } else if (type == 'gene') {
    # tidy gene and variant info
    genes <- dplyr::tibble(
      gene = genes,
      variant = colnames(genotypes)
    )
    genes <- genes$gene
    
    # subset genotyp by pedigree and no variance
    sms <- purrr::map(ped, ~subsetGT(.x, genotypes))
    
    # make a list of sites, and
    # map to gene and family
    sites <- purrr::map_df(sms,
      function(x) {
      df <- dplyr::tibble(
        variant = colnames(x),
        index = 1:ncol(x),
        gene = genes[na.omit(match(colnames(x), colnames(genotypes)))]
      )
      },
      .id = 'family'
    )

    # split sites
    sites <- split(sites, sites$gene)
    # sites <- sites[1:5] # TODO: to be removed
    
    # create carrier sets, and calculate probability of sharing in each
    carrier_set <- purrr::map(ped, carrierSet)
    ppl <- purrr::map(carrier_set, ~.x$carrier.sets.prob)
    npl <- purrr::map(carrier_set, ~.x$carrier.sets.len)
    
    # calculate sharing probability per gene
    res <- purrr::map_df(
      sites, 
      function(x) {
        df <- RVgene(sms, ped,
                     pattern.prob.list = ppl, N.list = npl,
                     sites = x$index,
                     fams = x$family)
        
        dplyr::tibble(
          pvalue = df$p,
          pall = df$pall,
          potential_pvalue = df$potentialp,
          family = names(df$pall)
        )
        
      },
      .id = 'gene'
    )
  }
  return(res)
}

#' Get gene annotation from genotype file
#'
#' @param genotype_file A character. Path to genotype in vcf format
#'
#' @return A character vector
#'
#' @examples
#' # read gene annotation
#' genotype_file <- 'workflow/lib/extdata/combined.variants.vcf.gz'
#' 
#' getGenes(genotype_file)
#' 
#' @importFrom VariantAnnotation readInfo
#' 
#' @export
getGenes <- function(genotype_file) {
  as.character(
    VariantAnnotation::readInfo(genotype_file, 'SYMBOL')
  )
}

#' Read genotype and pedigree in a snp matrix
#'
#' @param ped_file A character. Path to file in pedigree format
#' @param genotype_file A character. Path to genotype in vcf format
#'
#' @return A list
#' @export
readPedGT <- function(ped_file, genotype_file) {
  # create pedigree object
  ped <- readPed(ped_file)
  
  # create famInfo object
  cols <- c('pedigree', 'member', 'father', 'mother', 'sex', 'affected')
  fam <- read_tsv(ped_file, col_names = cols)
  
  # create a snp matrix
  vcf <- readVcf(genotype_file)
  gt <- genotypeToSnpMatrix(vcf[, colnames(vcf) %in% fam$member])
  
  # add ped and fam to the snp matrix
  gt$fam <- fam[fam$member %in% colnames(vcf),]
  gt$ped <- ped
  
  return(gt)
}

#' Create combinations of carriers
#'
#' @param ped A pedigree object
#'
#' @return A list of sets of carriers. Items are carrier.sets, carrier.sets.len,
#' and carrier.sets.prob.
#'
#' @examples
#' # read pedigree file
#' fl <- 'workflow/lib/extdata/FACT0246.pedigree.tsv'
#' ped <- readPed(fl)
#' 
#' # create a carrier set
#' carrierSet(ped)
#' 
#' @export
carrierSet <- function(ped) {
  # extract carriers/affected
  carriers <- with(ped, id[affected == 1 & !is.na(affected)])
  
  # get combinations
  carrier.sets <- list()
  for (i in length(carriers):1) {
    carrier.sets <- c(carrier.sets, combn(carriers, i, simplify=FALSE))
  }
  
  # get length of each set
  carrier.sets.len <- list()
  for (i in 1:length(carrier.sets)) {
    carrier.sets.len <- unlist(c(carrier.sets.len, length(carrier.sets[[i]])))
  }
  
  # calculate probability of sharing in the carrier sets
  carrier.sets.prob <- list()
  for (i in 1:length(carrier.sets)) {
    carrier.sets.prob <- unlist(c(carrier.sets.prob, RVsharing(ped, carriers = carrier.sets[[i]])))
  }
  
  res <- list(
    carrier.sets = carrier.sets,
    carrier.sets.len = carrier.sets.len,
    carrier.sets.prob = carrier.sets.prob
  )
  
  return(res)
}

#' Get sharing pattern
#'
#' @param gt A genotype matrix
#' @param ped A pedigree object
#' @param sites A character. Names of variant/s
#' @param affected A logical (TRUE) to get sharing for affected subjects only
#' @param return A character. 'complete' for sharing by all subjects and
#' 'partial' for partial sharing in a subset
#'
#' @return A character. The names of subjects sharing a variant
#' 
#' @examples
#' # read pedigree file
#' fl <- 'workflow/lib/extdata/FACT0246.pedigree.tsv'
#' ped <- readPed(fl)
#' 
#' # load variants genotype
#' fl <- 'workflow/lib/extdata/FACT0246.rare.filtered.vcf.gz'
#' gt <- VariantAnnotation::readGeno(fl, 'GT')
#' 
#' sites <- rownames(gt)[1]
#' 
#' # returns NA when no complete sharing is present
#' getSharing(gt, ped, sites, TRUE, 'complete')
#' 
#' # returns subset of samples when no partial sharing is present
#' getSharing(gt, ped, sites, TRUE, 'partial')
#' 
#' # returns multiple variants sharing pattern when gene is used
#' sites <- rownames(gt)[c(3, 4)]
#' 
#' getSharing(gt, ped, sites, TRUE, 'complete')
#' 
#' @export
getSharing <- function(gt, ped, sites, affected = TRUE, return = c('complete', 'partial')) {
  if (affected) {
    ind <- with(ped, id[affected == 1 & !is.na(affected)])
  } else {
    ind <- with(ped, id[affected == 0 & !is.na(affected)])
  }
  # subset to sites
  m <- gt[sites, ]
  
  # m not ref
  if (is.matrix(m)) {
    n <- m[,ind] == '0/1' | m[,ind] == '0|1'
    n <- apply(n, 2, any)
  } else {
    n <- m[ind] == '0/1' | m[ind] == '0|1'
  }
  
  # return rowname
  if (return == 'complete') {
    # when all is not ref
    if (all(n)) {
      res <- names(n)
    } else {
      res <- NA
    }
  } else if (return == 'partial') {
    # when sum is not ref
    if (sum(n) > 0) {
      res <- names(n)[n]
    } else {
      res <- NA
    }
  }
  
  return(res)
}

#' Apply sharing calculations
#'
#' @param gt A genotype matrix
#' @param ped A pedigree object
#' @param gene 
#' @param return A character. 'complete', 'partial' or gene
#'
#' @return A tibble
#' 
#' @importFrom purrr map
#' @importFrom dplyr tibble
#' 
#' @examples
#' # read pedigree file
#' fl <- 'workflow/lib/extdata/FACT0246.pedigree.tsv'
#' ped <- readPed(fl)
#' 
#' # load variants genotype
#' fl <- 'workflow/lib/extdata/FACT0246.rare.filtered.vcf.gz'
#' gt <- VariantAnnotation::readGeno(fl, 'GT')
#' 
#' # extract gene info
#' gene <- as.character(readInfo(fl, 'SYMBOL'))
#' 
#' # returns complete, partial and gene sharing info
#' applySharing(gt, ped, gene, 'complete')
#' applySharing(gt, ped, gene, 'partial')
#' applySharing(gt, ped, gene, 'gene')
#' 
#' @export
applySharing <- function(gt, ped, gene, return = c('complete', 'partial', 'gene')) {
  # get affected and un affected subjects
  aff <- with(ped, id[affected == 1 & !is.na(affected)])
  notaff <- with(ped, id[affected == 0 & !is.na(affected)])
  
  # get partial sharing
  if (return == 'partial') {
    # get sites
    sites <- rownames(gt)
    names(sites) <- gene
    
    # get sharing in affected subjects and subset to positive
    affected <- map(sites, ~getSharing(gt, ped, .x, TRUE, 'partial'))
    
    # remove completely shared variants
    ind <- !is.na(affected) & lengths(affected) < length(aff)
  } else {
    # get complete sharing
    # get sites
    sites <- rownames(gt)
    names(sites) <- gene
    
    # get sharing in affected subjects and subset to positive
    affected <- map(sites, ~getSharing(gt, ped, .x, TRUE, 'complete'))
    ind <- !is.na(affected)
    
    if (return == 'gene') {
      # remove completely shared variants
      complete <- rownames(gt) %in% sites[ind]
      sites <- rownames(gt)[!complete]
      gene <- gene[!complete]
      
      # get sites
      sites <- split(sites, gene)
      sites <- sites[lengths(sites) > 1]
      
      affected <- map(sites, ~getSharing(gt, ped, .x, TRUE, 'complete'))
      ind <- !is.na(affected)
    }
  }
  
  
  # subset affected
  affected <- affected[ind]
  n.affected <- lengths(affected)
  
  # subset site
  sites2 <- sites[ind]
  
  # get sharing in affected
  notaffected <- map(sites2, ~getSharing(gt, ped, .x, FALSE, 'partial'))
  n.notaffected <- lengths(notaffected)
  
  # return results
  res <- tibble(
    gene = names(sites2),
    variant = map_chr(sites2, ~paste(.x, collapse = ',')),
    family = with(ped, unique(famid)),
    n.affected = n.affected,
    n.notaffected = n.notaffected,
    affected = map_chr(affected, ~paste(.x, collapse = ',')),
    notaffected = map_chr(notaffected, ~paste(.x, collapse = ','))
  )
  
  return(res)
}

#' Calculates the probability of completely shared variants/genes
#'
#' @param ped A list of pedigree objects
#' @param variants A list of data.frames of sharing produced by applySharing
#'
#' @return A data.frame. Modified variants file combined
#' 
#' @importFrom dplyr tibble
#' @importFrom dplyr left_join
#' @importFrom purrr map_df
#' @importFrom RVS RVsharing
#' @importFrom RVS multipleFamilyPValue
#' 
#' @examples
#' # load pedigrees
#' fls <- c(
#'  'workflow/lib/extdata/FACT0246.pedigree.tsv',
#'  'workflow/lib/extdata/FACT5731.pedigree.tsv',
#'  'workflow/lib/extdata/FACT6322.pedigree.tsv'
#' )
#' 
#' peds <- purrr::map(fls, readPed)
#' names(peds) <- c('FACT0246', 'FACT5731', 'FACT6322')
#' 
#' # load variants
#' fls2 <- c(
#' 'workflow/lib/extdata/FACT0246.gene.rare.sharing.tsv',
#' 'workflow/lib/extdata/FACT5731.gene.rare.sharing.tsv',
#' 'workflow/lib/extdata/FACT6322.gene.rare.sharing.tsv'
#' )
#' 
#' variants <- map(fls2, read_tsv)
#' 
#' # apply sharing probability 
#' res <- RVCompleteSharing(peds, variants)
#' res
#' @export
RVCompleteSharing <- function(ped, variants) {
  # remove empty input
  ind <- map_int(variants, nrow) > 0
  
  if (sum(ind) > 0) {
    ped <- ped[ind]
    variants <- variants[ind]
    variants <- bind_rows(variants)
    
    # sharing probability for each family
    sharing_prob <- RVsharing(ped)
    
    # sharing probability and pvalue for single variants single family
    sharing_prob_df <- tibble(
      family = names(sharing_prob),
      prob = sharing_prob,
      p.value = sharing_prob
    )
    
    # sharing probability across families
    # as the sum of all sharing probabilities across families at most as large as 
    # the sharing probability observed.
    sharing_pall_df <- with(variants, split(family, variant)) %>%
      map_df(~{
        v <- names(sharing_prob) %in% .x
        names(v) <- names(sharing_prob)
        p <- multipleFamilyPValue(sharing_prob, v)
        tibble(p.all = p)
      }, .id = 'variant')
    
    # write results
    res <- variants %>%
      left_join(sharing_prob_df) %>%
      left_join(sharing_pall_df)
  } else {
    res <- tibble(
      gene = as.character(), variant = as.character(), family = as.character(),
      n.affected = as.numeric(), n.notaffected = as.numeric(),
      affected = as.character(), notaffected = as.character(),
      prob = as.numeric(), p.value = as.numeric(), p.all = as.numeric()
    )
  }
  
  return(res)
}

#' Calculates the probability of partially shared variants/genes
#'
#' @param ped A pedigree object
#' @param pat A pattern of sharing produced by 
#' @param tidy 
#'
#' @return
#' 
#' @examples
#' # load pedigrees
#' fls <- c(
#'  'workflow/lib/extdata/FACT0246.pedigree.tsv',
#'  'workflow/lib/extdata/FACT5731.pedigree.tsv',
#'  'workflow/lib/extdata/FACT6322.pedigree.tsv'
#' )
#' 
#' peds <- purrr::map(fls, readPed)
#' names(peds) <- c('FACT0246', 'FACT5731', 'FACT6322')
#' 
#' # load variants
#' fls <- c(
#' 'workflow/lib/extdata/FACT0246.partial.rare.sharing.tsv',
#' 'workflow/lib/extdata/FACT5731.partial.rare.sharing.tsv',
#' 'workflow/lib/extdata/FACT6322.partial.rare.sharing.tsv'
#' )
#' 
#' variants <- map(fls, readr::read_tsv)
#' 
#' # apply partial sharing probability
#' RVPartialSharing(peds, variants)
#' 
#' @importFrom purrr map_int
#' @importFrom purrr map2_df
#' @importFrom purrr map
#' @importFrom purrr imap
#' @importFrom dplyr tibble
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom RVS RVsharing
#'
#' @export
RVPartialSharing <- function(ped, variants) {
  # remove empty input
  ind <- map_int(variants, nrow) > 0
  ped <- ped[ind]
  variants <- variants[ind]

  res <- map2_df(
    ped, variants,
    function(p, v) {
      # extract sharing pattern
      pat <- with(v, split(affected, variant)) %>%
        map(function(x) strsplit(x, ','))
      
      # create carrier sets
      carriers <- carrierSet(p)
      
      # calculate the probability and pvalue of observed sharing
      prob <- imap(
        pat, 
        function(x, .y) { 
          tibble(
            variant = .y,
            prob = RVsharing(p, carriers = x),
            p.value = sum(carriers$carrier.sets.prob[carriers$carrier.sets.prob <= prob & carriers$carrier.sets.len >= length(x)])
          )
        }
      )
      prob <- bind_rows(prob)
      left_join(v, prob)
    }
  )
  
  return(res)
}
