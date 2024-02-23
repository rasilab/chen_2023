"""Workflow for counting barcodes.

  :Date: 10 Jan 2023
"""

# useful libraries
import os
import pandas as pd
import re
import itertools as it


# configuration specific to this analysis
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment = "#", dtype=object)
print(sample_annotations)

linkage = "../../small_8xdicodon_rqcdel_linkage/data/filtered_barcodes/filtered_barcodes_all_linkage.csv"

# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    raw_barcode_counts = expand('../data/raw_barcode_counts/{sample_id}.csv', 
      sample_id=sample_annotations['sample_id']),
    linked_barcode_counts = expand('../data/linked_barcode_counts/{sample_id}.csv', 
      sample_id=sample_annotations['sample_id']),
   

def get_fastq(wildcards):
  """This function returns fastq file for sample name
  """
  sample_id = sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'sample_id'].item()
  filenames = [f'../data/fastq/{filename}' 
      for filename in filter(lambda x: re.search(f'_{sample_id}_', x) and re.search('_R1_', x), 
      os.listdir('../data/fastq/'))]
  return filenames

rule count:
  """Count each barcode
  """
  input: get_fastq
  output: '../data/raw_barcode_counts/{sample_id}.csv'
  log: '../data/raw_barcode_counts/{sample_id}.log'
  params:
    barcode_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'barcode_read'].tolist()[0],
    barcode_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'barcode_start'].tolist()[0],
    barcode_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'barcode_length'].tolist()[0],
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail
    export TMPDIR=$PWD
    awk -v OFS=, '
    BEGIN {{
      print "barcode","counts";
    }}
    NR%4 == 2 {{
      counts[substr(${params.barcode_read}, {params.barcode_start}, {params.barcode_length})]++
    }}
    END {{
      for (barcode in counts) print barcode, counts[barcode] | "sort -k2nr -t,"
    }}
    ' \
    {input} 1> {output} 2> {log}
    """

rule subset_to_linked_barcodes:
  """Subset barcode counts to only those barcodes identified in linkage sequencing
  """
  input:
    barcode_count_file = '../data/raw_barcode_counts/{sample_id}.csv',
    barcode_linkage_file = linkage
  output:
    linked_barcode_count_file = '../data/linked_barcode_counts/{sample_id}.csv'
  params:
    # column containing barcode sequence in linkage file
    barcode_col = 3
  log: '../data/linked_barcode_counts/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    awk -F, -v OFS=, '
      NR == 1 {{
        print "barcode_count",$0
      }};
      # read in linkage file
      NR == FNR && NR > 1 {{
        linkages[${params.barcode_col}] = $0;
        next
      }};
      # write barcode count along with insert-barcode linkage data
      {{
        if ($1 in linkages) print $2, linkages[$1]
      }}
      ' \
      {input.barcode_linkage_file} {input.barcode_count_file} \
      1> {output.linked_barcode_count_file} 2> {log}
    """