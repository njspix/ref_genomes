nextflow.enable.dsl=2

process download {
    executor 'pbs'
    queue params.queue
    cpus 1
    time '72h'

    input:
        tuple val(name), val(target) 

    output:
        tuple val(name), path("*.2bit"), emit: twobit optional true
        tuple val(name), path("*.gz"), emit: gzip optional true

    script:
    """
        wget -O ${name} ${target}
    """
}

process unpack {
    executor 'pbs'
    queue params.queue
    cpus 1
    time '72h'
    conda params.conda_path+'/2bit'

    input:
        tuple val(name), path(twobit) 

    output:
        tuple val(id), path("*.f*a") 

    script:
    id = file(name).getSimpleName().toString()
    """
        twoBitToFa ${twobit} ${id}.fa
    """
}

process unzip {
    executor 'pbs'
    queue params.queue
    cpus 1
    time '72h'

    input:
        tuple val(name), path(gzip) 

    output:
        tuple val(id), path("*.fa") 

    script:
    id = file(name).getSimpleName().toString()
    """
        gunzip -f ${gzip}
    """
}

process adjoin_lambda {
    executor 'pbs'
    queue params.queue
    cpus 1
    time '72h'

    input:
        tuple val(ref_name), path(ref), path(lambda) 

    output:
        tuple val(new_ref_name), path("*.fa") 

    script:
    new_ref_name = ref_name+"_lambda"
    """
        cat ${ref} ${lambda} > ${new_ref_name}.fa
    """
}

process biscuit_index {
    executor 'pbs'
    queue params.queue
    cpus 1
    time '72h'
    conda params.conda_path+'/seq_qc'
    publishDir "${params.out_dir}/${id}/biscuit", mode: 'copy'

    input:
        tuple val(id), path(ref)

    output:
        path("*", includeInputs: true)

    script:
    """
        biscuit index ${ref}
    """
}

process bwa_index {
    executor 'pbs'
    queue params.queue
    cpus 1
    time '72h'
    conda params.conda_path+'/bwa'
    publishDir "${params.out_dir}/${id}/bwa", mode: 'copy'

    input:
        tuple val(id), path(ref)

    output:
        path("*", includeInputs: true)

    script:
    """
        bwa index ${ref}
    """
}


params.queue = 'laird'
params.conda_path = '/varidata/research/projects/laird/nathan/tools/miniconda3/envs/'
params.out_dir = "."

download_files = Channel.of(
    ['lambda.fa.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/840/245/GCA_000840245.1_ViralProj14204/GCA_000840245.1_ViralProj14204_genomic.fna.gz'],
    ['hg38-noalt-d1.fa.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz'],
    ['hg38-alt-d1.fa.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz'],
    ['mm39.2bit', 'ftp://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit'],
    ['mm10.2bit', 'ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit']
)

workflow {
    download(download_files)
    unpack(download.out.twobit)
    unzip(download.out.gzip)

    unpack.out.mix(unzip.out).branch{
        lambda: it[0] =~ /lambda/
        other: true
        }.set{refs}
    refs.other.combine(refs.lambda.map{it-> return it[1]}).set{ refs1 }

    adjoin_lambda(refs1)
    biscuit_index(adjoin_lambda.out.mix(refs.other))
    bwa_index(adjoin_lambda.out.mix(refs.other))
}

