bedtools	sort	-i	db/ensembl_gene_hg19.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_5890.bed	-b	stdin	>	distance_r2/nearest_gene.tsv
bedtools	sort	-i	db/ucsc_tbBrowser_ensGene_CDS_hg19.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_5890.bed	-b	stdin	>	distance_r2/cds_gene.tsv
bedtools	sort	-i	db/wgEncodeRegTfbsClusteredV3.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_5890.bed	-b	stdin	>	distance_r2/encode_tfbs.tsv
bedtools	sort	-i	db/ucsc_annot.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_5890.bed	-b	stdin	>	distance_r2/ucsc_annot.tsv