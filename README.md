## Table of Contents
- [Abstract](#abstract)
- [Instructions for running the code repo](#running-the-code)
- [Data](data/)
  - Includes flow cytometry raw data that is also available at http://flowrepository.org/id/FR-FCM-Z6QH
- [Code for designing endogenous fragments library](analysis/library_design/endogenous_fragments/)
- [Code for processing flow cytometry data and regenerating figure panels](analysis/flow_cytometry/wt_hel2_8xdicodon/scripts)
- [Code for linking barcodes and codon pair inserts in original 8× dicodon library](analysis/barcodeseq/8xdicodon_linkage/scripts/)
- [Code for linking barcodes and codon pair inserts in frameshifted 8× dicodon library](analysis/barcodeseq/frameshifted_8xdicodon_linkage/)
- [Code for linking barcodes and codon pair inserts small-scale validation library](analysis/barcodeseq/mini_8xdicodon_linkage/scripts/)
- Code for counting barcodes and regenerating figures:
  - [codon pair library wild-type cells](analysis/barcodeseq/wt_mrna_grna/scripts/)
  - [codon pair library *hel2Δ* and *syh1Δ* cells](analysis/barcodeseq/hel2_syh1_mrna_grna/scripts/)
  - [frameshifted codon pair library](analysis/barcodeseq/wt_frameshifted_mrna_grna/scripts/)
  - [small-scale validation library](analysis/barcodeseq/wt_hel2_mini_pool/scripts/)
  - [codon pair library under glucose depletion](analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/)
- [Code for counting insert-UMI pairs in (FK)~8~ DMS library](analysis/barcodeseq/wt_hel2_fk8_dms/scripts/)
- [Code to regenerate all figure panels](analysis/run_all_ipynb_scripts.smk)

# Massively parallel identification of sequence motifs triggering ribosome-associated mRNA quality control

**Katharine Chen**<sup>1,2</sup>, **Heungwon Park**<sup>1</sup>, **Arvind Rasi Subramaniam**<sup>1,†</sup>

<sup>1</sup> Basic Sciences Division and Computational Biology Section of the Public
Health Sciences Division, Fred Hutchinson Cancer Center, Seattle, WA
98109, USA <br/>
<sup>2</sup> Molecular and Cellular Biology Program, University of Washington,
Seattle, WA 98195, USA <br/>

<sup>†</sup> Corresponding author: <rasi@fredhutch.org>

## Abstract

Decay of mRNAs can be triggered by ribosome slowdown at stretches of rare codons or positively charged amino acids.
However, the full diversity of sequence motifs that trigger co-translational mRNA decay is poorly understood.
To comprehensively identify sequence motifs that trigger mRNA decay, we use a massively parallel reporter assay to measure the effect of all possible combinations of codon pairs on mRNA levels in S. cerevisiae.
In addition to known mRNA-destabilizing sequences, we identify several dipeptide repeats whose translation reduces mRNA levels. 
These include combinations of positively charged and bulky residues, as well as proline-glycine and proline-aspartic acid dipeptide repeats.
Genetic deletion of the ribosome collision sensor Hel2 rescues the mRNA effects of these motifs, suggesting that they trigger ribosome slowdown and activate the ribosome-associated quality control (RQC) pathway.
Deep mutational scanning of an mRNA-destabilizing dipeptide repeat reveals a complex relationship between the charge, bulkiness, and location of amino acid residues in conferring mRNA instability.
Finally, we show that the mRNA effects of codon pairs are predictive of the effects of endogenous sequences.
Our work highlights the complexity of sequence motifs driving co-translational mRNA decay in eukaryotes, and presents a high-throughput approach to dissect their requirements at the codon level.

## Running the code
- To run this on a cluster with singularity containers, do:
```
module load singularity # for fred hutch cluster
conda activate snakemake # this is a minimal conda env that has snakemake-minimal and pandas for invoking snakefile
sh run_everything.sh
```

- The ```run_everything.sh``` file will:
  - Download FASTQ files from SRA
  - Run all linkage sequencing, barcode sequencing, and insert sequencing [code](analysis/barcodeseq)
  - Run [code](analysis/library_design/endogenous_fragments/scripts/run_analysis.smk) to design of the endogenous fragments library
  - Run [code](analysis/run_all_ipynb_scripts.smk) to regenerate figure panels
    - This will also run the code to process flow cytometry data

## Docker containers
- [R](https://github.com/rasilab/r/pkgs/container/r)
- [python](https://github.com/rasilab/python/pkgs/container/python)
- [R and python](https://github.com/rasilab/r_python/pkgs/container/r_python)
