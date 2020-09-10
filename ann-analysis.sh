#!/bin/bash

# 9/9/2020 VERSION
# This file is generated and maintained by Seungsoo Kim.

WORK_DIR="/data/20.08-PHD3"
BASE_BED=$WORK_DIR"/gwas_biomart_1569.bed"
SH_FILE="dist_20.08-PHD3.sh"
ANN_PATH="/data/db_gwas"

ROAD_FILE=$ANN_PATH"/roadmap_meta.tsv"
REG_DIR=$ANN_PATH"/regulome"

### Archive code ###
#LNC_DIR=$ANN_PATH"/lncrna"
#GTEX_RDS=$ANN_PATH"/gtex_signif_5e-08_Ensgid_dt.rds"

### RUN FUNCTIONS ###
# DO NOT CHANGE BELLOW THIS CODE.

# 1. Overlapping the functional annotations
## 1-1. Measuring distance by bedtools
Rscript postgwas-exe.r \
  --utils bash \
  --base $BASE_BED \
  --out  $WORK_DIR \
  --ann  $ANN_PATH
mv $SH_FILE /data/$SH_FILE
bash /data/$SH_FILE

## 1-2. UCSC gene tags separation
Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/genome_dist/ucsc_annot.tsv \
  --out    $WORK_DIR/summary \
  --infotype ucsc \
  > $WORK_DIR"/log_ucsc.txt"

## 1-3. Roadmap merge & cell type
Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/roadmap_dist/roadmap_enh_merge.tsv \
  --out    $WORK_DIR/summary \
  > $WORK_DIR"/log_roadmap_merge.txt"

Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/roadmap_dist \
  --out    $WORK_DIR/roadmap_over \
  --meta   $ROAD_FILE \
  > $WORK_DIR"/log_roadmap_sub.txt"

## 1-4. RegulomeDB
Rscript postgwas-exe.r \
  --dbfilt regulome \
  --base   $BASE_BED \
  --out    $WORK_DIR \
  --regulm $REG_DIR \
  > $WORK_DIR"/log_regulome.txt"

## 1-5. lncRNASNP2
Rscript postgwas-exe.r \
  --dbfilt lnc_ovl \
  --base   $BASE_BED \
  --out    $WORK_DIR \
  --lncrna $ANN_PATH"/db_gwas.db" \
  > $WORK_DIR"/log_lncrnasnp.txt"

## 1-6. GTEx eQTL genes <- Short of RAM memory..
Rscript postgwas-exe.r \
  --dbfilt gtex_ovl \
  --base   $BASE_BED \
  --out    $WORK_DIR \
  --gtex   $ANN_PATH"/db_gwas.db" \
  > $WORK_DIR"/log_gtex.txt"


# 2. Union list to summary
## 2-1. ENCODE Tfbs
Rscript postgwas-exe.r \
  --dbfilt dist \
  --base   $WORK_DIR/genome_dist/encode_tfbs.tsv \
  --out    $WORK_DIR/summary \
  > $WORK_DIR"/log_encode.txt"

## 2-2. Roadmap Enhancers
Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/roadmap_over \
  --out    $WORK_DIR/summary \
  --sub_dir FALSE \
  --uni_save TRUE \
  > $WORK_DIR"/log_roadmap_summ.txt"

Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/roadmap_over \
  --out    $WORK_DIR/summary \
  --sub_dir TRUE \
  --uni_save TRUE \
  > $WORK_DIR"/log_roadmap_summ_sub.txt"

## 2-3. GTEx eQTLs union
Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/gtex_eqtl \
  --out    $WORK_DIR/gtex_eqtl \
  --sub_dir FALSE \
  --uni_save TRUE \
  > $WORK_DIR"/log_gtex_summ.txt"

### END FUNCTIONS ###