bedtools	sort	-i	./1a_suhlab/YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled.pf_peaks.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1a_suhlab_dist_37/Scbeta_YZ_pr1_pooled.pf_peaks.tsv
bedtools	sort	-i	./1a_suhlab/YY005-beta-cell-1_R1.nodup.tn5.pr2_pooled.pf_peaks.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1a_suhlab_dist_37/Scbeta_YZ_pr2_pooled.pf_peaks.tsv
bedtools	sort	-i	./1a_suhlab/YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1a_suhlab_dist_37/SCbeta_YZ_pseudo_pool_ATAC.tsv
