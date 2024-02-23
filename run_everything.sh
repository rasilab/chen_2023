base_folder=$(pwd)
echo "Running analysis in $base_folder folder"
## For singularity on the rhino cluster only, run this from command line before running this script
# module load Singularity
# conda activate snakemake 

echo "Downloading SRA annotations from GEO"
cd $base_folder/analysis/barcodeseq/
sh submit_cluster.sh "--snakefile" download_sra_annotations.smk $@

echo "Downloading FASTQ annotations from SRA"
sh submit_cluster.sh "--snakefile" download_fastq.smk $@

echo "Running 8xdicodon linkage sequencing analysis"
cd $base_folder/analysis/barcodeseq/8xdicodon_linkage/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running WT mRNA and gDNA barcode sequencing analysis"
cd $base_folder/analysis/barcodeseq/wt_mrna_grna/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running Hel2-del/Syh1-del mRNA and gDNA barcode sequencing analysis"
cd $base_folder/analysis/barcodeseq/hel2_syh1_mrna_grna/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running WT/Hel2-del glucose deprivation mRNA and gDNA barcode sequencing analysis"
cd $base_folder/analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running WT/Hel2-del FK8 DMS mRNA and gDNA insert sequencing analysis"
cd $base_folder/analysis/barcodeseq/wt_hel2_fk8_dms/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running mini 8xdicodon library linkage sequencing analysis"
cd $base_folder/analysis/barcodeseq/mini_8xdicodon_linkage/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running mini WT/Hel2-del mRNA and gDNA barcode sequencing analysis"
cd $base_folder/analysis/barcodeseq/wt_hel2_mini_pool/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running frameshifted 8xdicodon library linkage sequencing analysis"
cd $base_folder/analysis/barcodeseq/frameshifted_8xdicodon_linkage/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running frameshifted WT mRNA and gDNA barcode sequencing analysis"
cd $base_folder/analysis/barcodeseq/wt_frameshifted_mrna_grna/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running endogenous fragments library design scripts"
cd $base_folder/analysis/library_design/endogenous_fragments/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running small 8xdicodon library linkage sequencing analysis"
cd $base_folder/analysis/barcodeseq/small_8xdicodon_rqcdel_linkage/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running small 8xdicodon library in RQC-del strains mRNA and gDNA barcode sequencing analysis"
cd $base_folder/analysis/barcodeseq/small_8xdicodon_rqcdel_mrna_grna/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running Upf1-del FK8 DMS mRNA and gDNA insert sequencing analysis"
cd $base_folder/analysis/barcodeseq/upf1_fk8_dms/scripts
sh submit_cluster.sh "--snakefile" run_analysis.smk $@

echo "Running all plotting notebooks to regenerate figures"
cd $base_folder/analysis
sh submit_cluster.sh "--snakefile" run_all_ipynb_scripts.smk $@