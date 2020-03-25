bedtools	sort	-i	./3_tanlab/Th1_enh.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./3_tanlab_dist/Th1_enh.tsv
bedtools	sort	-i	./3_tanlab/Th1_enh_control.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./3_tanlab_dist/Th1_enh_control.tsv
bedtools	sort	-i	./3_tanlab/Th1_enh_t1d.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./3_tanlab_dist/Th1_enh_t1d.tsv
bedtools	sort	-i	./3_tanlab/Treg_enh.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./3_tanlab_dist/Treg_enh.tsv
bedtools	sort	-i	./3_tanlab/Treg_enh_control.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./3_tanlab_dist/Treg_enh_control.tsv
bedtools	sort	-i	./3_tanlab/Treg_enh_t1d.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./3_tanlab_dist/Treg_enh_t1d.tsv
