for i in {1..96} ##second round of plate-based indexing barcode, here index1-96
do
mkdir ~/workdir/index${i}
cd ~/workdir/index${i}

cat > index${i}_TCR.sh<< EOF

#PBS -q core24
#PBS -l nodes=1:ppn=2,mem=40gb,walltime=40:00:00
#PBS -e tcr_${i}.err
#PBS -o tcr_${i}.log

cd ~/workdir/index${i}
~/software/path2cellarnger/cellranger vdj --id=tcr${i} \
                 --reference=~/reference_for_cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-4.0.0 \
                 --fastqs=~/path2tcr_fastq_index${i} \
                 --sample=index${i} \
                 --localcores=2 \
                 --localmem=20

EOF
dsub index${i}_TCR.sh
done

