rule all:
  input:
    "get_sgd_orf_annotations.nbconvert.ipynb",
    "../data/sgd/cds_coding.fasta",
    "calculate_stall_scores.nbconvert.ipynb",
    "../data/stallscores/posbulky_motif_scores_sgd_orfs.tsv",
    "../data/stallscores/dipeptide_scores_sgd_orfs.tsv",
    "evaluate_motif_scores.nbconvert.ipynb",
    "../data/motifs/stall_motifs_controls.tsv"

rule write_sgd_annotations_fasta:
  input:
    db = "../db/S288C_reference_genome_R64-3-1_20210421/sgd_gene_annotations.gff3",
    notebook = "get_sgd_orf_annotations.ipynb",
  output:
    "get_sgd_orf_annotations.nbconvert.ipynb",
    "../data/sgd/cds_coding.fasta"
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """

rule calculate_stall_scores:
  input:
    notebook = "calculate_stall_scores.ipynb",
    fasta = '../data/sgd/cds_coding.fasta'
  output:
    "calculate_stall_scores.nbconvert.ipynb",
    "../data/stallscores/posbulky_motif_scores_sgd_orfs.tsv",
    "../data/stallscores/dipeptide_scores_sgd_orfs.tsv"
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=python3 {input.notebook}
    """

rule subset_genes_to_test:
  input:
    notebook = "evaluate_motif_scores.ipynb",
    db = "../db/S288C_reference_genome_R64-3-1_20210421/sgd_gene_annotations.gff3",
    table1 = "../data/stallscores/posbulky_motif_scores_sgd_orfs.tsv",
    table2 = "../data/stallscores/dipeptide_scores_sgd_orfs.tsv"
  output:
    "evaluate_motif_scores.nbconvert.ipynb",
    "../data/motifs/stall_motifs_controls.tsv"
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """