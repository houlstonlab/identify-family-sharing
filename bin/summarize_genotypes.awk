#!/usr/bin/awk -f

BEGIN {
    FS = "\t"; OFS = "\t"
}

{
    gene = $1
    variant = $2
    het_count = 0
    hom_count = 0
    het_samples = ""
    hom_samples = ""
    for (i = 3; i <= NF; i++) {
        split($i, arr, "=")
        sample = arr[1]
        genotype = arr[2]
        if (genotype == "0/1" || genotype == "0|1" || genotype == "1/0" || genotype == "1|0") {
            het_count++
            het_samples = het_samples ? het_samples "," sample : sample
        } else if (genotype == "1/1" || genotype == "1|1") {
            hom_count++
            hom_samples = hom_samples ? hom_samples "," sample : sample
        }
    }
    print gene, variant, het_count, het_samples, hom_count, hom_samples
}