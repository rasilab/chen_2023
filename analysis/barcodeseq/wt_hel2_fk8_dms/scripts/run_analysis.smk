"""Workflow for counting inserts and umis

  :Author: Arvind Rasi Subramaniam
  :Edited by: Katharine Chen
  :Date: June 2023
"""

# useful libraries
import os
import pandas as pd
import re
import itertools as it


# configuration specific to this analysis
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment = "#", dtype=object)
sra_annotations = pd.read_table("../../../../annotations/sra_annotations.tsv")

# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    insert_umi = expand('../data/insert_umi/{sample_name}.csv', 
      sample_name=sample_annotations['sample_name']),
    insert_umi_counts = expand('../data/insert_umi_counts/{sample_name}.csv', 
      sample_name=sample_annotations['sample_name']),
    insert_annotations = '../annotations/insert_annotations.csv',
    annotated_insert_umi_counts = expand('../data/annotated_insert_umi_counts/{sample_name}.csv', 
      sample_name=sample_annotations['sample_name']),

   
def get_fastq_file_for_sample_name(wildcards):
  """This function gets the SRR file based on the sample_id in `sample_annotations`"""
  sample_id = sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'sample_id'].item()
  srr = sra_annotations.loc[sra_annotations['sample_id'] == sample_id, 'srr'].item()
  filename = f'../../../../data/fastq/{srr}.fastq'
  return filename

rule extract_and_tabulate_all_insert_umi:
  """Extract and tabulate insert and umis
  """
  input: get_fastq_file_for_sample_name
  output: '../data/insert_umi/{sample_name}.csv'
  params:
    umi_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'umi_read'].tolist()[0],
    umi_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'umi_start'].tolist()[0],
    umi_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'umi_length'].tolist()[0],
    insert_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'insert_read'].tolist()[0],
    insert_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'insert_start'].tolist()[0],
    insert_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'insert_length'].tolist()[0],
  log: '../data/insert_umi/{sample_name}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail;
    # concatenate all input files line by line with tab separator
    paste {input} | \
    awk -v OFS="," '
      BEGIN {{print "insert","umi"}}
      NR%4 == 2  \
      {{ 
        print substr(${params.insert_read}, {params.insert_start}, {params.insert_length}), \
              substr(${params.umi_read}, {params.umi_start}, {params.umi_length}) 
      }}
    ' 1> {output}  2> {log}
    """


rule count_insert_umi_combinations:
  """Count number of reads for each combination of insert and umi
  """
  input: '../data/insert_umi/{sample_name}.csv'
  output: '../data/insert_umi_counts/{sample_name}.csv'
  log: '../data/insert_umi_counts/{sample_name}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail
    export TMPDIR=$PWD
    # skip header and skip sequences with N
      sed -n '2,${{/N/!p}}' {input} | \
      sort | \
      uniq -c | \
      sort -k1nr | \
      awk -v OFS="," 'BEGIN {{print "insert","umi","count"}} {{print $2,$1}}' \
      1> {output} 2> {log}
    """


rule get_variant_annotations:
  """Get variant annotations
  """
  input: 
    dms_variant_oligos = '../annotations/fk_dms_opool.csv',
    notebook = 'get_variant_annotations.ipynb',
  output: '../annotations/insert_annotations.csv'
  log: '../annotations/get_variant_annotations.log'
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
    notebook={input.notebook}
    script="${{notebook/.ipynb/.r}}"
    Rscript ${{script}} \
      --dms_variants_file={input.dms_variant_oligos} \
      --output_file={output} \
      &> {log}
    """


rule subset_to_annotated_inserts:
  """Subset insert-umi counts to only annotated inserts
  """
  input:
    count_file = '../data/insert_umi_counts/{sample_name}.csv',
    insert_annotations = '../annotations/insert_annotations.csv',
  output:
    annotated_insert_count_file = '../data/annotated_insert_umi_counts/{sample_name}.csv'
  params:
    # column in insert annotations file that contains the insert number and sequence
    insert_num_column = 1,
    insert_seq_column = 2,
  log: '../data/annotated_insert_umi_counts/{sample_name}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    awk -F, -v OFS=, '
      # set up parameters and header
      BEGIN {{
        insert_num_col={params.insert_num_column}; 
        insert_seq_col={params.insert_seq_column};
        print "insert_num", "umi", "count"
      }}
      # read in all insert annotations
      NR == FNR {{
        inserts[$insert_seq_col] = $insert_num_col; 
        next
      }}
      # write if insert is present
      {{
        if ($1 in inserts) {{
          print inserts[$1], $2, $3;
        }}
      }}
      ' \
      {input.insert_annotations} {input.count_file} 1>> {output.annotated_insert_count_file}
    """