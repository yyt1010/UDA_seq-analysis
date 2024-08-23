for i in {1..96} ####second round of plate-based indexing barcode, here index1-96. Running on HPC respective
do
mkdir ~/workdir/index${i}
cd ~/workdir/index${i}

cat > index${i}_cellranger.sh << EOF

#PBS -q core24
#PBS -l nodes=1:ppn=2,mem=40gb,walltime=40:00:00
#PBS -e im_${i}.err
#PBS -o im_${i}.log


cd ~/workdir/index${i}
~/software/path2cellarnger/cellranger count --id=post_pbmc${i} \
--fastqs=~/path2fastq_index${i} \
--sample=RNA_index${i} \
--localmem=30 \
--chemistry="fiveprime" \
--localcores=12 \
--include-introns=true \
--transcriptome=~/reference_for_cellranger/refdata-cellranger-GRCh38-3.0.0
EOF
dsub index${i}_cellranger.sh


