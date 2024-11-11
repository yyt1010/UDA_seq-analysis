# UDA-seq: universal droplet microfluidics-based combinatorial indexing for massive-scale multimodal single-cell sequencing
## UDA-seq workflow


<img width="515" alt="image" src="https://github.com/user-attachments/assets/ef334233-cf5f-41f7-aaa4-eae20e6e0160">

## Processing of UDA-seq data
### step 1: mapping of scRNA-seq&scATAC-seq multiome, scRNA-seq, VDJ-seq and CROP-seq
For n(typically n=96) round2 plate-based indexing barcodes, sequencing reads from different index were mapped to genome separately. Use script in 0_Preprocess
```
sh step1_mapping_for_paried_RNA_ATAC.sh  ##for scRNA-seq&scATAC-seq multiome
sh step1_mapping_for_RNA.sh for ##scRNA-seq
sh step1_mapping_for_TCR.sh ##for TCR in VDJ-seq
sh step1_mapping_for_BCR.sh ##for BCR in VDJ-seq
```
After that, bam files and feature-matrix grouped by n(typically n=96) index will be genrated.
### step2:decoding samples by demuxlet (https://www.nature.com/articles/nbt.4042). 
Using software demuxlet. Samples were decoding by script in 0_Preprocess
```
sh step2_demuxlet_donor.sh ## major inputs were bam files and donors VCF file.
```
### step 3: generating overall feature x cell (eg. gene x cell) matrix
by merging all round2 index barcode matrix and demuxlet, then filtering low quality cells out
```
Rscript step3_generate_overall_matrix.R
```
### step4: Scrublet doublet identification (Scrublet: Computational Identification of Cell
Doublets in Single-Cell Transcriptomic Data.https://doi.org/10.1016/j.cels.2018.11.005)
Use script in 0_Preprocess
```
python step4_Scrublet.py
```
### step5: get final feature x cell matatrix for downstream analysis.
```
Rscript step5_final_matrix.R
```
### finally, downstream analyis scripts for uda-seq paper were uploaded in code directory grouped by figures.
