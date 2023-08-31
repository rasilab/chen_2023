"""Workflow for getting linkage

  :Author: Arvind Rasi Subramaniam
  :Date: 5 Jan 2023
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
sra_annotations = pd.read_table("../../../../annotations/sra_annotations.tsv")


# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    '../annotations/insert_annotations.csv',
    insert_barcodes = expand('../data/insert_barcodes/{sample_id}.csv', 
      sample_id=sample_annotations['sample_id']),
    insert_barcode_counts = expand('../data/insert_barcode_counts/{sample_id}.csv', 
      sample_id=sample_annotations['sample_id']),
    annotated_insert_barcode_counts = expand('../data/annotated_insert_barcode_counts/{sample_id}.csv', 
      sample_id=sample_annotations['sample_id']),
    ref_vs_ref_align = expand('../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode1.bam',
      sample_id=sample_annotations['sample_id']),
    filtered_barcodes = expand('../data/filtered_barcodes/{sample_id}.csv',
      sample_id=sample_annotations['sample_id']),
   

def get_fastq_file_for_sample_name(wildcards):
  """This function gets the R1 or R2 file depending on the `insert_read` column of `sample_annotations`
  """
  sample_id = sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'sample_id'].item()
  srr = sra_annotations.loc[sra_annotations['sample_id'] == sample_id, 'srr'].item()
  filename = f'../../../../data/fastq/{srr}.fastq'
  return filename


rule extract_and_tabulate_all_insert_barcodes:
  """Extract and tabulate insert and barcodes
  """
  input: get_fastq_file_for_sample_name
  output: '../data/insert_barcodes/{sample_id}.csv'
  params:
    barcode_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'barcode_read'].tolist()[0],
    barcode_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'barcode_start'].tolist()[0],
    barcode_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'barcode_length'].tolist()[0],
    insert_read = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'insert_read'].tolist()[0],
    insert_start = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'insert_start'].tolist()[0],
    insert_length = lambda wildcards: sample_annotations.loc[sample_annotations['sample_id'] == wildcards.sample_id, 'insert_length'].tolist()[0],
  log: '../data/insert_barcodes/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail;
    # concatenate all input files line by line with tab separator
    paste {input} | \
    awk -v OFS="," '
      BEGIN {{print "insert","barcode"}}
      NR%4 == 2  \
      {{ 
        print substr(${params.insert_read}, {params.insert_start}, {params.insert_length}), \
              substr(${params.barcode_read}, {params.barcode_start}, {params.barcode_length}) 
      }}
    ' 1> {output}  2> {log}
    """

rule count_insert_barcode_combinations:
  """Count number of reads for each combination of insert and barcodes
  """
  input: '../data/insert_barcodes/{sample_id}.csv'
  output: '../data/insert_barcode_counts/{sample_id}.csv'
  log: '../data/insert_barcode_counts/{sample_id}.log'
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
      awk -v OFS="," 'BEGIN {{print "insert","barcode","count"}} {{print $2,$1}}' \
      1> {output} 2> {log}
    """

rule get_variant_annotations:
  """Get variant annotations
  """
  input: 
    dms_variant_oligos = '../annotations/sdd1_fk_dms_opool.csv',
    non_dms_variant_oligos = '../annotations/stall_motifs_controls.tsv',
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
      --non_dms_variants_file={input.non_dms_variant_oligos} \
      --output_file={output} \
      &> {log}
    """

rule subset_to_annotated_inserts:
  """Subset insert-barcode counts to only annotated inserts
  """
  input:
    count_file = '../data/insert_barcode_counts/{sample_id}.csv',
    insert_annotations = '../annotations/insert_annotations.csv',
  output:
    annotated_insert_count_file = '../data/annotated_insert_barcode_counts/{sample_id}.csv'
  params:
    # column in insert annotations file that contains the insert number and sequence
    insert_num_column = 1,
    insert_seq_column = 2,
  log: '../data/annotated_insert_barcode_counts/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    awk -F, -v OFS=, '
      # set up parameters and header
      BEGIN {{
        barcode_num=1;
        insert_num_col={params.insert_num_column}; 
        insert_seq_col={params.insert_seq_column};
        print "barcode_num", "insert_num", "barcode", "count"
      }}
      # read in all insert annotations
      NR == FNR {{
        inserts[$insert_seq_col] = $insert_num_col; 
        next
      }}
      # write if insert is present
      {{
        for (insert in inserts) 
          if ($1 ~ insert) {{
            print barcode_num,inserts[insert], $2, $3;
            barcode_num++
          }}
      }}
      ' \
      {input.insert_annotations} {input.count_file} 1>> {output.annotated_insert_count_file}
    """

rule align_barcodes_against_themselves:
  """Align barcodes against themselves to find multialigners
  """
  input:
    '../data/annotated_insert_barcode_counts/{sample_name}.csv'
  output:
    sam = temp('../data/ref_vs_ref_alignments/{sample_name}/alignment_barcode1.sam'),
    bam = '../data/ref_vs_ref_alignments/{sample_name}/alignment_barcode1.bam',
    fasta = '../data/ref_vs_ref_alignments/{sample_name}/reference_barcode1.fasta',
  log:
    align = '../data/ref_vs_ref_alignments/{sample_name}/align_barcode1.log',
    build = '../data/ref_vs_ref_alignments/{sample_name}/build_barcode1.log',
  params:
    bowtie_index = '../data/ref_vs_ref_alignments/{sample_name}/reference_barcode1',
    read_count_cutoff = 25
  threads: 18
  container: 'docker://ghcr.io/rasilab/bowtie2:2.4.5'
  shell:
    """
    # write the input file to a fasta file of barcodes with name as barcode_num col from input_file
    awk -F, 'NR > 1 {{if ($4 >= {params.read_count_cutoff}) print ">" $1 "\\n" $3}}' {input} > {output.fasta}
    # create a bowtie reference of the barcodes
    bowtie2-build {output.fasta} {params.bowtie_index} 2> {log.build}
    # align against itself
    bowtie2 --threads {threads} -L 19 -N 1 --all --norc --no-unal -f -x {params.bowtie_index} -U {output.fasta} > {output.sam}  2> {log.align}
    # convert to BAM
    samtools view -@ {threads} -b {output.sam} > {output.bam}.tmp
    # sort
    samtools sort -@ {threads} {output.bam}.tmp > {output.bam}
    sleep 10
    # index
    samtools index -@ {threads} {output.bam}
    # remove unsorted bam
    rm {output.bam}.tmp
    """

rule filter_barcodes:
  """Filter barcodes to remove clashes and sequencing errors and produce a final list
  """
  input:
    bam1 = '../data/ref_vs_ref_alignments/{sample_name}/alignment_barcode1.bam',
    counts = '../data/annotated_insert_barcode_counts/{sample_name}.csv',
    notebook = 'filter_barcodes.ipynb',
  output:
    '../data/filtered_barcodes/{sample_name}.csv',
  params:
    read_count_cutoff = 25
  log:
    '../data/filtered_barcodes/{sample_name}.log',
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}

    notebook={input.notebook}
    script="${{notebook/.ipynb/.r}}"

    Rscript ${{script}} {input.bam1} {input.counts} {params.read_count_cutoff} {output} &> {log}
    """