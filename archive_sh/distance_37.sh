bedtools	sort	-i	db/ensembl_gene_hg19.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	distance_37/nearest_gene.tsv			
bedtools	sort	-i	db/ucsc_tbBrowser_ensGene_CDS_hg19.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	distance_37/cds_gene.tsv			
bedtools	sort	-i	db/wgEncodeRegTfbsClusteredV3.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	distance_37/encode_tfbs.tsv			
