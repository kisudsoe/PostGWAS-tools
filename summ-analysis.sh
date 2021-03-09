#!/bin/bash

WORK_DIR="data"
ANN_GWAS=$WORK_DIR"/input_bed/gwas_biomart_5892_uq.tsv"
ANN_GTEX=$WORK_DIR"/gtex_signif_3938.tsv"

### RUN FUNCTIONS ###
# DO NOT CHANGE BELLOW THIS CODE.

printf "1. Summary table.. "
Rscript postgwas-exe.r \
  --dbvenn summ \
  --base   $WORK_DIR/summary $WORK_DIR/gtex_eqtl $WORK_DIR/encode_over \
  --out    $WORK_DIR \
  --sub_dir FALSE \
  --uni_save FALSE \
  --ann_gwas $ANN_GWAS \
  --ann_near $WORK_DIR"/genome_dist/nearest_gene.tsv" \
  --ann_gtex $ANN_GTEX \
  > $WORK_DIR"/log_summary.txt"
#--ann_encd $WORK_DIR"/genome_dist/encode_tfbs.tsv" \
printf "done\n"

printf "2. Roadmap summary.. "
Rscript src/roadmap_summary.r \
  --pivot \
  --f_roadmap_dist $WORK_DIR/roadmap_dist_25status \
  --out            $WORK_DIR \
  > $WORK_DIR"/log_roadmap_5_summ_25status.txt"
printf "done\n"

printf "3. ENCODE summary.. "
Rscript src/encode_summary.r \
	--f_encode_dist $WORK_DIR"/genome_dist/encode_tfbs.tsv" \
  --out $WORK_DIR \
  > $WORK_DIR"/log_encode_summ.txt"
printf "done\n"

  ### END FUNCTION ###