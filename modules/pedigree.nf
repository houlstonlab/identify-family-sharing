process PEDIGREE {
    tag "${famid}"

    label 'simple'

    publishDir("${params.output_dir}/pedigree", mode: 'copy')

    input:
    tuple val(famid), val(id), val(fid), val(mid), val(sex), val(aff), val(family)

          
    output:
    tuple val(famid), path("${famid}.ped")

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
        <(echo -e "${family.join('\n')}") \
        > ${famid}.ped
    """
}
