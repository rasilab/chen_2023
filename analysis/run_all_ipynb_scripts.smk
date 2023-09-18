rule all:
  input:
    "barcodeseq/wt_mrna_grna/scripts/plot_aggregate_effects.nbconvert.ipynb",
    "barcodeseq/wt_mrna_grna/scripts/plot_dipeptide_effects.nbconvert.ipynb",
    "barcodeseq/wt_mrna_grna/scripts/plot_supp_alignment_stats.nbconvert.ipynb",
    "barcodeseq/wt_mrna_grna/scripts/plot_supplemental_missing_data.nbconvert.ipynb",
    "barcodeseq/wt_mrna_grna/tables/destabilized_wt.csv",
    "barcodeseq/wt_mrna_grna/tables/barcode_insert_counts.tsv.gz",
    "flow_cytometry/wt_hel2_8xdicodon/scripts/plot_figure2_flow.nbconvert.ipynb",
    "barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.nbconvert.ipynb",
    "barcodeseq/hel2_syh1_mrna_grna/scripts/plot_hel2_syh1_dipeptide_effects.nbconvert.ipynb",
    "barcodeseq/hel2_syh1_mrna_grna/scripts/plot_supp_aln_qc.nbconvert.ipynb",
    "barcodeseq/wt_hel2_fk8_dms/scripts/plot_variant_effects.nbconvert.ipynb",

rule plot_aggregate_effects:
  input:
    "barcodeseq/wt_mrna_grna/scripts/plot_aggregate_effects.ipynb",
  output:
    "barcodeseq/wt_mrna_grna/scripts/plot_aggregate_effects.nbconvert.ipynb",
    wt="barcodeseq/wt_mrna_grna/tables/barcode_insert_counts.tsv.gz",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule plot_dipeptide_effects:
  input:
    notebook="barcodeseq/wt_mrna_grna/scripts/plot_dipeptide_effects.ipynb",
  output:
    "barcodeseq/wt_mrna_grna/scripts/plot_dipeptide_effects.nbconvert.ipynb",
    "barcodeseq/wt_mrna_grna/tables/destabilized_wt.csv",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0' 
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """

rule plot_wt_barcodeseq_aln_stats:
  input:
    "barcodeseq/wt_mrna_grna/scripts/plot_supp_alignment_stats.ipynb",
  output:
    "barcodeseq/wt_mrna_grna/scripts/plot_supp_alignment_stats.nbconvert.ipynb",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule plot_missing_inserts_summary:
  input:
    notebook="barcodeseq/wt_mrna_grna/scripts/plot_supplemental_missing_data.ipynb",
    wt="barcodeseq/wt_mrna_grna/tables/destabilized_wt.csv",
  output:
    "barcodeseq/wt_mrna_grna/scripts/plot_supplemental_missing_data.nbconvert.ipynb",
    "barcodeseq/wt_mrna_grna/tables/missing_dicodons_plasmid_library.csv"
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """   

rule plot_figure2_flow:
  input:
    "flow_cytometry/wt_hel2_8xdicodon/scripts/plot_figure2_flow.ipynb",
  output:
    "flow_cytometry/wt_hel2_8xdicodon/scripts/plot_figure2_flow.nbconvert.ipynb",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'  
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """

rule plot_translation_effects:
  input:
    notebook="barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb",
    table="barcodeseq/wt_mrna_grna/tables/destabilized_wt.csv",
  output:
    "barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.nbconvert.ipynb",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0' 
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """

rule plot_hel2_syh1_effects:
  input:
    notebook="barcodeseq/hel2_syh1_mrna_grna/scripts/plot_hel2_syh1_dipeptide_effects.ipynb",
    table="barcodeseq/wt_mrna_grna/tables/destabilized_wt.csv",
    wt="barcodeseq/wt_mrna_grna/tables/barcode_insert_counts.tsv.gz"
  output:
    "barcodeseq/hel2_syh1_mrna_grna/scripts/plot_hel2_syh1_dipeptide_effects.nbconvert.ipynb",
    hel2="barcodeseq/hel2_syh1_mrna_grna/tables/barcode_insert_counts.tsv.gz",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """

rule plot_hel2_syh1_aln_stats_qc:
  input:
    notebook="barcodeseq/hel2_syh1_mrna_grna/scripts/plot_supp_aln_qc.ipynb",
    wt="barcodeseq/wt_mrna_grna/tables/barcode_insert_counts.tsv.gz"
  output:
    "barcodeseq/hel2_syh1_mrna_grna/scripts/plot_supp_aln_qc.nbconvert.ipynb",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input.notebook}
    """
  
rule plot_fk8_dms_data:
  input:
    "barcodeseq/wt_hel2_fk8_dms/scripts/plot_variant_effects.ipynb",
  output:
    "barcodeseq/wt_hel2_fk8_dms/scripts/plot_variant_effects.nbconvert.ipynb",
  container: 'docker://ghcr.io/rasilab/r_python:1.1.0'
  shell:
    """
    export JUPYTER_DATA_DIR=$(pwd)
    export JUPYTER_CONFIG_DIR=$(pwd)
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=ir {input}
    """