for i in {1..96}  ##second round of plate-based indexing barcode, here index1-96
do
cd ~/demuxlet/Index_${i}

cat > demuxlet_${i}.sh << EOF
#!/bin/sh
#PBS -q core40
#PBS -l mem=60gb,nodes=1:ppn=1,walltime=20:10:00
#PBS -e 1_transfer.err
#PBS -o 1_transfer.txt
#STOR -s /pnas,/software,/home,/xtdisk,/p300s,/gpfs

cd ~/demuxlet/Index_${i}
gunzip ~/cellranger_outs/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
~/conda/anaconda3/envs/popscale/bin/popscle demuxlet --sam ./atac_gex.bam \  ##atac_gex bam was megered by atac_bam and gex_bam.
--vcf ~/SNPs/donor_vcf.vcf \
--field GP --out ${i}_result_GP --group-list ~/cellranger_outs/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

gzip ~/cellranger_outs/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

EOF

jsub -f demuxlet_${i}.sh

done

