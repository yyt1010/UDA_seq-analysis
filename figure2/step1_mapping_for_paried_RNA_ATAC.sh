for i in {1..96} ##second round of plate-based indexing barcode, here index1-96
do
mkdir ~/workdir/Paired_mapping/index_${i}
cd ~/workdir/Paired_mapping/index_${i}
cat > libraries.csv << EOF
fastqs,sample,library_type
~/workdir/RNA_fastq_path,index_${i},Gene Expression
~/workdir/ATAC_fastq_path,index_${i},Chromatin Accessibility
EOF

cat > cellranger_arc_index_${i}.sh << EOF
#!/bin/sh
#PBS -q c56m256g
#PBS -l mem=60gb,nodes=1:ppn=8,walltime=100:10:00
#PBS -e 1_transfer.err
#PBS -o 1_transfer.txt
#STOR -s /work_dir,/software



cd ~/workdir/Paired_mapping/index_${i}
~/software/software/cellranger-arc-2.0.1/cellranger-arc count --id=kidney_${i} \
                       --reference=~/reference_10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=~/workdir/Paired_mapping/index_${i}/libraries.csv \
                       --localcores=8 \
                       --localmem=60
EOF
jsub -f cellranger_arc_index_${i}.sh
done

