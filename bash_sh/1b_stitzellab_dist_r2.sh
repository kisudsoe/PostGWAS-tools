bedtools sort -i	./1b_stitzellab/GSM3333912_EndoC_BH1_ATACseq_broadPeak.fdr0.05.noBlacklist.bed	| bedtools closest -d -a	gwas_hg19_biomart_5890.bed	 -b stdin >	./1b_stitzellab_dist_r2/EndoC_BH1_ATAC_broadPeak.tsv
bedtools sort -i	./1b_stitzellab/GSE118588_EndoC_BH1_ChromHMM_annotations.bed	| bedtools closest -d -a	gwas_hg19_biomart_5890.bed	 -b stdin >	./1b_stitzellab_dist_r2/EndoC_BH1_ChromHMM.tsv
bedtools sort -i	./1b_stitzellab/GSM3333898_Human_Islet_HiC_HiCCUPS_loops.bed	| bedtools closest -d -a	gwas_hg19_biomart_5890.bed	 -b stdin >	./1b_stitzellab_dist_r2/Primary_Islet_HiC_gwas.tsv
bedtools sort -i	./1b_stitzellab/GSM3333916_EndoC_BH1_HiC_HiCCUPS_loops.bed	| bedtools closest -d -a	gwas_hg19_biomart_5890.bed	 -b stdin >	./1b_stitzellab_dist_r2/EndoC_BH1_HiC_gwas.tsv
bedtools sort -i	./1b_stitzellab/GSM3333898_Human_Islet_HiC_HiCCUPS_loops.bed	| bedtools closest -d -a	./db/ucsc_annot_sort.bed	 -b stdin >	./1b_stitzellab_dist_r2/Primary_Islet_HiC_gene.tsv
bedtools sort -i	./1b_stitzellab/GSM3333916_EndoC_BH1_HiC_HiCCUPS_loops.bed	| bedtools closest -d -a	./db/ucsc_annot_sort.bed	 -b stdin >	./1b_stitzellab_dist_r2/EndoC_BH1_HiC_gene.tsv
