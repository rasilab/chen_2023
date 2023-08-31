import pandas as pd
import itertools as it

sra_annotations = pd.read_table("../../annotations/sra_annotations.tsv")

rule all:
  input:
    singularity_container = "sra-tools_3.0.0.sif",
    fastq = [f'../../data/fastq/{srr}.fastq' for srr in sra_annotations['srr']]

rule download_image:
    """Download singularity image
    """
    # We are manually downloading this image because 
    # the image does not have bash as currently required by Snakemake
    # see https://github.com/snakemake/snakemake/issues/1521
    input: 
    output: 'sra-tools_3.0.0.sif'
    envmodules: "Singularity"
    shell:
        """
        singularity pull --disable-cache --force docker://docker.io/ncbi/sra-tools:3.0.0
        """        

rule get_concatenated_paired_end_fastq:
    """Download fastq from SRA"""
    # We are manually executing this image because 
    # the image does not have bash as currently required by Snakemake
    # see https://github.com/snakemake/snakemake/issues/1521
    input: 
      sra = 'sra-tools_3.0.0.sif',
      #annotations = '../../annotations/sra_annotations.tsv'
    output: '../../data/fastq/{srr}.fastq'
    params:
        directory = '../../data/fastq/'
    threads: 36
    envmodules: "Singularity"
    shell:
        """
        singularity exec -B $(dirname $(dirname $(pwd))) sra-tools_3.0.0.sif \
        fasterq-dump --concatenate-reads --threads {threads} --outdir {params.directory} {wildcards.srr}
        """