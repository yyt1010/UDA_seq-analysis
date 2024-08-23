for i in {1..96} ##second round of plate-based indexing barcode, here index1-96
do
mkdir ~/workdir/index${i}
cd ~/workdir/index${i}
cat > library.csv << EOF
fastqs,sample,library_type
~/workdir/path2Crop_seq_RNA_part,Crop${i},Gene Expression
~/workdir/path2Crop_seq_sgRNA_part,CropS${i},CRISPR Guide Capture
EOF

cat > index${i}_crop.sh << EOF
#PBS -q core24
#PBS -l nodes=1:ppn=2,mem=40gb,walltime=40:00:00
#PBS -e crop_${i}.err
#PBS -o crop_${i}.log


cd  ~/workdir/index${i}
~/software/path2cellarnger/cellranger count --id=Crop${i} \
--localmem=30 \
--localcores=4 \
--include-introns=true \
--force-cells=400 \ 
--transcriptome=~/reference_for_cellranger/refdata-cellranger-GRCh38-3.0.0 \
--feature-ref=~/feature_barcode_information/library_ref.csv \
--libraries=~/workdir/index${i}/library.csv \
--chemistry='SC3Pv3'
EOF
dsub index${i}_crop.sh

done

