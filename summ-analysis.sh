#!/bin/bash

WORK_DIR="/data/20.08-PHD3"
ANN_GWAS=$WORK_DIR"/gwas_biomart_fill_1569.tsv"
ANN_GTEX=$WORK_DIR"/gtex_signif_1159.tsv"

# 3. Summary table
Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/summary $WORK_DIR/gtex_eqtl \
  --out    $WORK_DIR \
  --sub_dir FALSE \
  --uni_save FALSE \
  --ann_gwas $ANN_GWAS \
  --ann_near $WORK_DIR"/genome_dist/nearest_gene.tsv" \
  --ann_encd $WORK_DIR"/genome_dist/encode_tfbs.tsv" \
  --ann_gtex $ANN_GTEX \
  > $WORK_DIR"/log_summary.txt"