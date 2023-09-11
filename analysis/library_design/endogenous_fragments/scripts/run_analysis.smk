rule all:
  input:
    "yeast_fragments_design.nbconvert.ipynb",
    '../tables/yeastorffrags.csv'

rule extract_fragments:
  input:
    db = '../db/saccharomyces_cerevisiae_R64-1-1_20110208.gff',
    fasta = '../db/orf_coding_all_R64-1-1_20110203.fasta',
    expression = '../weinberg2016/GSE53313_Cerevisiae_RNA_RPF.txt',
    notebook = "yeast_fragments_design.ipynb",
  output:
    "yeast_fragments_design.nbconvert.ipynb",
    '../tables/yeastorffrags.csv'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=python3 {input.notebook}
    """