process PEDIGREE {
    tag "${famid}"

    label 'simple'

    publishDir("${params.output_dir}/pedigree", mode: 'copy')

    input:
    tuple val(famid), val(id), val(fid), val(mid), val(sex), val(aff)

          
    output:
    tuple val(famid),
          path("${famid}.samples.txt"),
          path("${famid}.parents.txt"),
          path("${famid}.sex.txt"),
          path("${famid}.aff.txt"),
          path("${famid}.ped")

    script:
    """
    #!/bin/bash
    # Write the values as columns in a text file
    paste \
        <(echo -e "${id.join('\n')}") \
        <(echo -e "${fid.join('\n')}") \
        <(echo -e "${mid.join('\n')}") \
        <(echo -e "${sex.join('\n')}") \
        <(echo -e "${aff.join('\n')}") \
        > ${famid}.ped
    
    cat ${famid}.ped | awk '{print \$1}' > ${famid}.samples.txt
    cat ${famid}.ped | awk '{print "${famid}\t"\$1,\$2,\$3}' > ${famid}.parents.txt
    cat ${famid}.ped | awk '{print "${famid}\t"\$1,\$4}' > ${famid}.sex.txt
    cat ${famid}.ped | awk '{print "${famid}\t"\$1,\$5}' > ${famid}.aff.txt
    """
}
