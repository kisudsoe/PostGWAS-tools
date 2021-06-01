#!/bin/bash

# 5/31/2021 VERSION
# This file is generated and maintained by Seungsoo Kim.

WORK_DIR="/data"
BASE_BED=$WORK_DIR"/input_bed/gwas_biomart_5892.bed"
SH_FILE=$WORK_DIR/"dist_data.sh"
ANN_PATH="/data/db_gwas"

ROAD_FILE=$ANN_PATH"/roadmap_meta.tsv"
REG_DIR=$ANN_PATH"/regulome"
ENCODE_FILE=$ANN_PATH"/wgEncodeRegTfbsClusteredWithCellsV3.bed"

### Archive code ###
#LNC_DIR=$ANN_PATH"/lncrna"
#GTEX_RDS=$ANN_PATH"/gtex_signif_5e-08_Ensgid_dt.rds"

### RUN FUNCTIONS ###
# DO NOT CHANGE BELLOW THIS CODE.

printf "\n1. Overlapping the functional annotations\n"
printf "1-1. Measuring distance by bedtools\n"
Rscript postgwas-exe.r \
  --utils bash \
  --base $BASE_BED \
  --out  $WORK_DIR \
  --ann  $ANN_PATH
#mv $SH_FILE $WORK_DIR/$SH_FILE
printf "  Run bedtools commands. This step is taking a while. "
cat $SH_FILE | parallel
printf "done\n"

printf "1-2. UCSC gene tags separation "
Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/genome_dist/ucsc_annot.tsv \
  --out    $WORK_DIR/summary \
  --infotype ucsc \
  > $WORK_DIR"/log_ucsc.txt"
printf "done\n"

printf "1-3. Roadmap merge & cell type, 1/2 "
Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/roadmap_dist/roadmap_enh_merge.tsv \
  --out    $WORK_DIR/summary \
  > $WORK_DIR"/log_roadmap_1_merge.txt"

printf "2/2 "
Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/roadmap_dist \
  --out    $WORK_DIR/roadmap_over \
  --meta   $ROAD_FILE \
  > $WORK_DIR"/log_roadmap_2_sub.txt"
printf "done\n"

printf "1-4. ENCODE Tfbs 1/2 "
Rscript src/enrich.r --split_tfbs \
  --tfbs   $WORK_DIR/genome_dist/encode_tfbs.tsv \
  --type   tsv \
  --out    $WORK_DIR/encode_dist \
  > $WORK_DIR"/log_encode_dist.txt"

printf "2/2 "
Rscript src/enrich.r --tfbs2bed \
  --tfbs   $WORK_DIR/encode_dist \
  --out    $WORK_DIR/encode_over \
  > $WORK_DIR"/log_encode_over.txt"
printf "done\n"

printf "1-5. RegulomeDB "
Rscript postgwas-exe.r \
  --dbfilt regulome \
  --base   $BASE_BED \
  --out    $WORK_DIR \
  --regulm $REG_DIR \
  > $WORK_DIR"/log_regulome.txt"
printf "done\n"

printf "1-6. lncRNASNP2 "
Rscript postgwas-exe.r \
  --dbfilt lnc_ovl \
  --base   $BASE_BED \
  --out    $WORK_DIR \
  --lncrna $ANN_PATH"/db_gwas.db" \
  > $WORK_DIR"/log_lncrnasnp.txt"
printf "done\n"

printf "1-7. GTEx eQTL genes "
Rscript postgwas-exe.r \
  --dbfilt gtex_ovl \
  --base   $BASE_BED \
  --out    $WORK_DIR \
  --gtex   $ANN_PATH"/db_gwas.db" \
  > $WORK_DIR"/log_gtex.txt"
printf "done\n\n"

printf "2. Union list to summary\n"
printf "2-1. ENCODE Tfbs "
Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/genome_dist/encode_tfbs.tsv \
  --out    $WORK_DIR/encode_dist \
  > $WORK_DIR"/log_encode.txt"
printf "done\n"

printf "2-2. Roadmap Enhancers, 1/2 "
Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/roadmap_over \
  --out    $WORK_DIR/summary \
  --sub_dir FALSE \
  --uni_save TRUE \
  > $WORK_DIR"/log_roadmap_3_summ.txt"

printf "2/2 "
Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/roadmap_over \
  --out    $WORK_DIR/summary \
  --sub_dir TRUE \
  --uni_save TRUE \
  > $WORK_DIR"/log_roadmap_4_summ_sub.txt"
printf "done\n"

printf "2-3. GTEx eQTLs union"
Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/gtex_eqtl \
  --out    $WORK_DIR/gtex_eqtl \
  --sub_dir FALSE \
  --uni_save TRUE \
  > $WORK_DIR"/log_gtex_summ.txt"
printf "done\n"

### END FUNCTIONS ###