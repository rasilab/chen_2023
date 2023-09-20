rule all:
  input:
    '../tables/yeastorffrags.csv'

rule extract_fragments:
  input:
    fasta = '../db/orf_coding_all_R64-1-1_20110203.fasta',
    annotated_expression = '../weinberg2016/GSE53313_Cerevisiae_RNA_RPF_annotated.csv',
    notebook = 'design_yeast_fragments.ipynb'
  output:
    '../tables/yeastorffrags.csv'
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """