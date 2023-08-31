"""Workflow for counting inserts and umis

  :Author: Arvind Rasi Subramaniam
  :Editted by: Katharine Chen
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

yeast_linkage_file = '../../8xdicodon_linkage/data/filtered_barcodes/yeast_cyto_linkage.csv'

# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    barcode_umi = expand('../data/barcode_umi/{sample_name}.csv', 
      sample_name=sample_annotations['sample_name']),
    barcode_umi_counts = expand('../data/barcode_umi_counts/{sample_name}.csv', 
      sample_name=sample_annotations['sample_name']),
    linked_barcode_counts = expand('../data/linked_barcode_counts/{sample_name}.csv', 
      sample_name=sample_annotations['sample_name']),
    linked_barcode_umi_counts = expand('../data/linked_barcode_umi_counts/{sample_name}.csv', 
      sample_name=sample_annotations['sample_name']),
   
def get_fastq_file_for_sample_name(wildcards):
  """This function gets the SRR file based on the sample_id in `sample_annotations`"""
  sample_id = sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'sample_id'].item()
  srr = sra_annotations.loc[sra_annotations['sample_id'] == sample_id, 'srr'].item()
  filename = f'../../../../data/fastq/{srr}.fastq'
  return filename

rule extract_and_tabulate_all_barcode_umi:
  """Extract and tabulate insert and umis
  """
  input: get_fastq_file_for_sample_name
  output: '../data/barcode_umi/{sample_name}.csv'
  params:
    barcode_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode_read'].tolist()[0],
    barcode_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode_start'].tolist()[0],
    barcode_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'barcode_length'].tolist()[0],
    umi_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'umi_read'].tolist()[0],
    umi_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'umi_start'].tolist()[0],
    umi_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_name'] == wildcards.sample_name, 'umi_length'].tolist()[0],
  log: '../data/barcode_umi/{sample_name}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail;
    # concatenate all input files line by line with tab separator
    paste {input} | \
    awk -v OFS="," '
      BEGIN {{print "barcode","umi"}}
      NR%4 == 2  \
      {{ 
        print substr(${params.barcode_read}, {params.barcode_start}, {params.barcode_length}), \
              substr(${params.umi_read}, {params.umi_start}, {params.umi_length}) 
      }}
    ' 1> {output}  2> {log}
    """


rule count_distinct_barcode_umi_combinations:
  """Count number of reads for each combination of barcodes and UMIs
  """
  input: '../data/barcode_umi/{sample_name}.csv'
  output: '../data/barcode_umi_counts/{sample_name}.csv'
  log: '../data/barcode_umi_counts/{sample_name}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail
    awk -F, -v OFS="," '
    # skip header and barcodes with Ns
    NR > 1 && $1 !~ /N/ {{
      array[$1","$2]++
    }}
    END {{
      print "barcode","umi","count";
      for (e in array) print e,array[e]
    }}
    ' {input} \
    1> {output} 2> {log}
    """

rule subset_to_linked_barcodes:
  """Subset barcode-UMI counts to only those barcodes identified in linkage sequencing
  """
  input:
    barcode_count_file = '../data/barcode_umi_counts/{sample_name}.csv',
    barcode_linkage_file = yeast_linkage_file,
  output:
    linked_barcode_count_file = '../data/linked_barcode_counts/{sample_name}.csv'
  params:
    # columns in annotations file
    barcode_num_column = 2,
    barcode_seq_column = 3,
  log: '../data/linked_barcode_counts/{sample_name}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    awk -F, -v OFS=, '
      # set up parameters and header
      BEGIN {{
        barcode_num_col={params.barcode_num_column}; 
        barcode_seq_col={params.barcode_seq_column};
        print "barcode_num", "umi", "count"
      }}
      # read in all barcodes from linkageseq
      NR == FNR {{
        barcodes[$barcode_seq_col] = $barcode_num_col; 
        next
      }}
      # write if barcode is present in linkageseq
      {{
        if ($1 in barcodes) {{
          print barcodes[$1], $2, $3;
        }}
      }}
      ' \
      {input.barcode_linkage_file} {input.barcode_count_file} 1>> {output.linked_barcode_count_file}
    """


rule count_umis_reads_per_barcode:
  """Count the number of UMIs and reads per barcode
  """
  input: 
    linked_barcode_count_file = '../data/linked_barcode_counts/{sample_name}.csv'
  output:
    umi_count_file = '../data/linked_barcode_umi_counts/{sample_name}.csv'
  log: '../data/linked_barcode_umi_counts/{sample_name}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    set +o pipefail
    export TMPDIR=$PWD
    awk -F, -v OFS="," ' 
      BEGIN {{print "barcode_num,umi_count,read_count"}} 
      NR > 1 {{
        array1[$1]++; 
        array2[$1] += $3
      }}
      END {{for (a in array1) print a,array1[a],array2[a]}}
      ' \
      {input.linked_barcode_count_file} > {output.umi_count_file}
    """