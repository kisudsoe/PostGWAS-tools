# Post-GWAS analysis: INS, SMARCE1, BACH2, RNASEH1

7 SNPs from Yousin's email (3/12/2020)

* rs689 - INS
* rs7221109 - SMARCE1
* rs10806425 - BACH2
* rs11755527 - BACH2
* rs55981318 - RNASEH1
* rs7607888 - RNASEH1
* rs1136545 - RNASEH1 (rs1136545 is not in 1000G reference panel.)

See the detailed results at `(Yousin) ATAC-seq/7 Locuses for beta cell ATAC-seq.pptx` file.

Summary figure for the 7 SNPs + 30 LD linked SNPs:

![](r2_ins_data/Fig 1.png)

# 1. From the 7 seed SNPs

## Downloading LDlink data

```CMD
Rscript postgwas-exe.r
  --ldlink down
  --base db_ins/ins_smarce1_bach2_ranseh1.tsv
  --out db_ins/ldlink
  --popul CEU TSI FIN GBR IBS
```

> ** Run function ldlink_down... 7..
> done
>   Files are moved to target folder:     db_ins/ldlink
> Job done: 2020-03-13 13:05:47 for 2.7 sec

## Filtering the LDlink data (r<sup>2</sup>>0.8)

* 1) r2 >0.6 and Dprime =1  <- The most stringent criteria.
* 2) r2 >0.8               <- Usual choice to define LD association.
* 3) Dprime =1
* 4) r2 >0.6 or Dprime =1
* 5) r2 >0.8

```CMD
Rscript postgwas-exe.r
  --ldlink filter
  --base db_ins/ins_smarce1_bach2_ranseh1.tsv db_ins/ldlink
  --out r2_ins_data
  --r2d 5
```

> ** Run function ldlink_filter...
> Read download files... 7
>   Read LDlink results           = [1] 3159   11
> Filtering by "r2 > 0.8":
>   Filtered data dimension       = [1] 35  3
>   Excluded no rsid elements     = [1] 0
> Basic summary of LDlink results:
>   SNP Tier 1                    = 7
>   SNP Tier 2                    = 30
>   SNP candidates                = 37
>   SNP source annotation table   = [1] 37  2
> Add annotations:
>   LD block annotation... [1] 5
> Search biomart for SNP coordinates:
>   Query SNPs            = [1] 37
>   Hg19 result table     = [1] 37  4
>   Hg38 result table     = [1] 37  4
>   Cytoband annotation... 37.. done
>   Merged table          = [1] 37 11
> Write file: r2_ins_data/gwas_biomart.tsv
> Job done: 2020-03-13 13:37:16 for 23.6 sec

## Generating BED files (hg19 and hg38)

```CMD
Rscript postgwas-exe.r
  --ldlink bed
  --base r2_ins_data/gwas_biomart.tsv
  --out r2_ins_data
```

> ** Run function ldlink_bed... [1] 37 11
> Write file:     r2_ins_data/gwas_hg19_biomart_37.bed
> Write file:     r2_ins_data/gwas_hg38_biomart_37.bed
> Job done: 2020-03-13 13:41:59 for 0.3 sec

# 2. Downloading annotation data

This step is already done.



# 3. Filtering the annotations

This step is already done.



# 4. Calculating distance from the annotations

## For Roadmap each cell type

See details in `0_roadmap_dist_37.sh` file:

```bash
bedtools	sort	-i	0_roadmap/roadmap_001_enh.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	0_roadmap_dist_37/roadmap_001_enh.tsv
...
bedtools	sort	-i	0_roadmap/roadmap_129_enh.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	0_roadmap_dist_37/roadmap_129_enh.tsv
```

## For Melton lab's SC-β-cell multiomic data

```bash
bedtools	sort	-i	./1_meltonlab/GSM4171636_PP2_ATAC_rep1.peaks.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1_meltonlab_dist_37/PP2_ATAC_rep1.tsv
...
bedtools	sort	-i	./1_meltonlab/GSE140500_SCbeta_TE.byH3K27ac.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1_meltonlab_dist_37/SCbeta_TE_enh.tsv
```

## For Suh lab (Yizhou)'s SC-β-cell ATAC-seq data

```bash
bedtools	sort	-i	./1a_suhlab/YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled.pf_peaks.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1a_suhlab_dist_37/Scbeta_YZ_pr1_pooled.pf_peaks.tsv
bedtools	sort	-i	./1a_suhlab/YY005-beta-cell-1_R1.nodup.tn5.pr2_pooled.pf_peaks.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1a_suhlab_dist_37/Scbeta_YZ_pr2_pooled.pf_peaks.tsv
bedtools	sort	-i	./1a_suhlab/YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./1a_suhlab_dist_37/SCbeta_YZ_pseudo_pool_ATAC.tsv
```

## For Pritchard lab's immune cell ATAC-seq atlas

```bahs
bedtools	sort	-i	./2_pritchardlab/Pritchard_X1001.Bulk_B.S.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./2_pritchardlab_dist_37/1001.Bulk_B.S.tsv
...
bedtools	sort	-i	./2_pritchardlab/Pritchard_X1011.Naive_Teffs.S.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_2003.bed	-b	stdin	>	./2_pritchardlab_dist_37/1011.Naive_Teffs.S.tsv
```

## For Tan lab's Th1 and Treg enhancer data

```bash
bedtools	sort	-i	./3_tanlab/Th1_enh.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./3_tanlab_dist_37/Th1_enh.tsv
bedtools	sort	-i	./3_tanlab/Th1_enh_control.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./3_tanlab_dist_37/Th1_enh_control.tsv
bedtools	sort	-i	./3_tanlab/Th1_enh_t1d.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./3_tanlab_dist_37/Th1_enh_t1d.tsv
bedtools	sort	-i	./3_tanlab/Treg_enh.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./3_tanlab_dist_37/Treg_enh.tsv
bedtools	sort	-i	./3_tanlab/Treg_enh_control.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./3_tanlab_dist_37/Treg_enh_control.tsv
bedtools	sort	-i	./3_tanlab/Treg_enh_t1d.bed	|	bedtools	closest	-d	-a	gwas_hg19_biomart_37.bed	-b	stdin	>	./3_tanlab_dist_37/Treg_enh_t1d.tsv
```



# 5. Overlapping the annotations

## For Roadmap each cell type

```CMD
Rscript postgwas-exe.r
  --dbfilt dist
  --base r2_ins_data/0_roadmap_dist_37
  --out r2_ins_data/0_roadmap_over_37
  --meta db_gwas/roadmap_meta.tsv
```

> ...
>
> File roadmap_129_enh... nrow= 37.. done
>   Annotations occupied by SNPs  = [1] 2
>   SNPs in annotations           = [1] 2
>   Write file: r2_ins_data/0_roadmap_over_37/BONE/snp_roadmap_129_enh_2.bed
> Job process: 0 sec
>
> Job done: 2020-03-13 18:29:32 for 4.8 sec

## For ENCODE Tfbs data

```CMD
Rscript postgwas-exe.r ^
	--dbfilt dist ^
	--base r2_ins_data/distance_37/encode_tfbs.tsv ^
	--out r2_ins_data/summary
```

> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 1
> File encode_tfbs... nrow= 146.. done
>   Annotations occupied by SNPs  = [1] 118
>   SNPs in annotations           = [1] 12
>   Write file: r2_ins_data/summary/snp_encode_tfbs_12.bed
>   1/1 r2_ins_data/distance_37/encode_tfbs.tsv
> Job process: 0.3 sec
>
> Job done: 2020-03-17 19:32:39 for 0.4 sec

## For Melton lab's data

```CMD
Rscript postgwas-exe.r
  --dbfilt dist
  --base r2_ins_data/1_meltonlab_dist_37
  --out r2_ins_data/1_meltonlab_over_37
  --meta db_gwas/meltonlab_meta.tsv
```

> ...
>
> File SCbeta_TE_enh... nrow= 37.. done
>   Annotations occupied by SNPs  = [1] 0
>   SNPs in annotations           = [1] 0
>   Write file: r2_ins_data/1_meltonlab_over_37/beta_enh/snp_SCbeta_TE_enh_0.bed
> Job process: 0 sec
>
> Job done: 2020-03-13 18:55:13 for 1.3 sec

### For Suh lab's data

```CMD
Rscript postgwas-exe.r
  --dbfilt dist
  --base r2_ins_data/1a_suhlab_dist_37
  --out r2_ins_data/1a_suhlab_over_37
```

> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 3
> File Scbeta_YZ_pr1_pooled.pf_peaks... nrow= 64.. done
>   Annotations occupied by SNPs  = [1] 2
>   SNPs in annotations           = [1] 2
>   Write file: r2_ins_data/1a_suhlab_over_37/snp_Scbeta_YZ_pr1_pooled.pf_peaks_2.bed
>   1/3 r2_ins_data/1a_suhlab_dist_37/Scbeta_YZ_pr1_pooled.pf_peaks.tsv
> Job process: 0.2 sec
>
> File Scbeta_YZ_pr2_pooled.pf_peaks... nrow= 53.. done
>   Annotations occupied by SNPs  = [1] 1
>   SNPs in annotations           = [1] 1
>   Write file: r2_ins_data/1a_suhlab_over_37/snp_Scbeta_YZ_pr2_pooled.pf_peaks_1.bed
>   2/3 r2_ins_data/1a_suhlab_dist_37/Scbeta_YZ_pr2_pooled.pf_peaks.tsv
> Job process: 0 sec
>
> File SCbeta_YZ_pseudo_pool_ATAC... nrow= 64.. done
>   Annotations occupied by SNPs  = [1] 2
>   SNPs in annotations           = [1] 2
>   Write file: r2_ins_data/1a_suhlab_over_37/snp_SCbeta_YZ_pseudo_pool_ATAC_2.bed
>   3/3 r2_ins_data/1a_suhlab_dist_37/SCbeta_YZ_pseudo_pool_ATAC.tsv
> Job process: 0.1 sec
>
> Job done: 2020-03-13 18:58:46 for 0.4 sec

## For Pritchard lab's data

```CMD
Rscript postgwas-exe.r
  --dbfilt dist
  --base r2_ins_data/2_pritchardlab_dist_37
  --out r2_ins_data/2_pritchardlab_over_37
  --meta db_gwas/pritchardlab_meta.tsv
```

> ...
>
> File 1011.Naive_Teffs.S... nrow= 37.. done
>   Annotations occupied by SNPs  = [1] 14
>   SNPs in annotations           = [1] 17
>   Write file: r2_ins_data/2_pritchardlab_over_37/T-CD4-cells_stim/snp_1011.Naive_Teffs.S_17.bed
> Job process: 0 sec
>
> Job done: 2020-03-13 21:30:01 for 7.3 sec

## For Tan lab's data

```CMD
Rscript postgwas-exe.r
  --dbfilt dist
  --base r2_ins_data/3_tanlab_dist_37
  --out r2_ins_data/3_tanlab_over_37
  --meta db_gwas/tanlab_meta.tsv
```

> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 6
>   Read metadata file dim        = [1] 6 2
> File Th1_enh... nrow= 37.. done
>   Annotations occupied by SNPs  = [1] 2
>   SNPs in annotations           = [1] 4
>   Write file: r2_ins_data/3_tanlab_over_37/Th1_enh/snp_Th1_enh_4.bed
>   1/6 r2_ins_data/3_tanlab_dist_37/Th1_enh.tsv
> Job process: 0.2 sec
>
> ...
>
> File Treg_enh_t1d... nrow= 37.. done
>   Annotations occupied by SNPs  = [1] 1
>   SNPs in annotations           = [1] 3
>   Write file: r2_ins_data/3_tanlab_over_37/Treg_enh/snp_Treg_enh_t1d_3.bed
>   6/6 r2_ins_data/3_tanlab_dist_37/Treg_enh_t1d.tsv
> Job process: 0 sec
>
> Job done: 2020-03-13 19:21:44 for 0.6 sec

# 6. Generating BED summary

## ENCODE Tfbs list

```CMD
Rscript postgwas-exe.r
  --dbfilt dist
  --base r2_ins_data/distance_37/encode_tfbs.tsv
  --out r2_ins_data/summary
```

> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 1
> File encode_tfbs... nrow= 146.. done
>   Annotations occupied by SNPs  = [1] 118
>   SNPs in annotations           = [1] 12
>   Write file: r2_ins_data/summary/snp_encode_tfbs_12.bed
>   1/1 r2_ins_data/distance_37/encode_tfbs.tsv
> Job process: 0.2 sec
>
> Job done: 2020-03-13 22:37:35 for 0.3 sec

## Roadmap union list

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/0_roadmap_over_37
  --out r2_ins_data/summary
  --sub_dir FALSE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 30 Files/folders input.
>   1 1 files in the ADRENAL
>   2 27 files in the BLOOD
>   3 1 files in the BONE
>   4 13 files in the BRAIN
>   5 3 files in the BREAST
>   6 1 files in the CERVIX
>   7 8 files in the ESC
>   8 9 files in the ESC_DERIVED
>   9 3 files in the FAT
>   10 3 files in the GI_COLON
>   11 2 files in the GI_DUODENUM
>   12 1 files in the GI_ESOPHAGUS
>   13 3 files in the GI_INTESTINE
>   14 3 files in the GI_RECTUM
>   15 4 files in the GI_STOMACH
>   16 4 files in the HEART
>   17 5 files in the IPSC
>   18 1 files in the KIDNEY
>   19 2 files in the LIVER
>   20 5 files in the LUNG
>   21 7 files in the MUSCLE
>   22 1 files in the MUSCLE_LEG
>   23 1 files in the OVARY
>   24 2 files in the PANCREAS
>   25 2 files in the PLACENTA
>   26 8 files in the SKIN
>   27 1 files in the SPLEEN
>   28 2 files in the STROMAL_CONNECTIVE
>   29 2 files in the THYMUS
>   30 2 files in the VASCULAR
> Total 127 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
> In addition: Warning message:
> package 'eulerr' was built under R version 3.6.2
>   [ERROR] r2_ins_data/0_roadmap_over_37/BREAST/snp_roadmap_119_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/GI_ESOPHAGUS/snp_roadmap_079_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/HEART/snp_roadmap_105_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/PANCREAS/snp_roadmap_098_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/SKIN/snp_roadmap_057_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/SKIN/snp_roadmap_058_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/SKIN/snp_roadmap_059_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/SKIN/snp_roadmap_061_enh_0.bed
>   Read 127 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1]  23 123
>   Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_23.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-13 20:53:04 for 3.4 sec

### union list by groups

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/0_roadmap_over_37
  --out r2_ins_data/summary
  --sub_dir TRUE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 30 Files/folders input.
>   1 sub_dir 1: 1 file(s) in the ADRENAL folder
>   2 sub_dir 2: 27 file(s) in the BLOOD folder
>   3 sub_dir 3: 1 file(s) in the BONE folder
>   4 sub_dir 4: 13 file(s) in the BRAIN folder
>   5 sub_dir 5: 3 file(s) in the BREAST folder
>   6 sub_dir 6: 1 file(s) in the CERVIX folder
>   7 sub_dir 7: 8 file(s) in the ESC folder
>   8 sub_dir 8: 9 file(s) in the ESC_DERIVED folder
>   9 sub_dir 9: 3 file(s) in the FAT folder
>   10 sub_dir 10: 3 file(s) in the GI_COLON folder
>   11 sub_dir 11: 2 file(s) in the GI_DUODENUM folder
>   12 sub_dir 12: 1 file(s) in the GI_ESOPHAGUS folder
>   13 sub_dir 13: 3 file(s) in the GI_INTESTINE folder
>   14 sub_dir 14: 3 file(s) in the GI_RECTUM folder
>   15 sub_dir 15: 4 file(s) in the GI_STOMACH folder
>   16 sub_dir 16: 4 file(s) in the HEART folder
>   17 sub_dir 17: 5 file(s) in the IPSC folder
>   18 sub_dir 18: 1 file(s) in the KIDNEY folder
>   19 sub_dir 19: 2 file(s) in the LIVER folder
>   20 sub_dir 20: 5 file(s) in the LUNG folder
>   21 sub_dir 21: 7 file(s) in the MUSCLE folder
>   22 sub_dir 22: 1 file(s) in the MUSCLE_LEG folder
>   23 sub_dir 23: 1 file(s) in the OVARY folder
>   24 sub_dir 24: 2 file(s) in the PANCREAS folder
>   25 sub_dir 25: 2 file(s) in the PLACENTA folder
>   26 sub_dir 26: 8 file(s) in the SKIN folder
>   27 sub_dir 27: 1 file(s) in the SPLEEN folder
>   28 sub_dir 28: 2 file(s) in the STROMAL_CONNECTIVE folder
>   29 sub_dir 29: 2 file(s) in the THYMUS folder
>   30 sub_dir 30: 2 file(s) in the VASCULAR folder
> Total 30 sub-folder(s) is/are input
> Total 0 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
>   1 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_ADRENAL_2.bed
>   2 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_BLOOD_21.bed
>   3 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_BONE_2.bed
>   4 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_BRAIN_4.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
> In addition: Warning message:
> package 'eulerr' was built under R version 3.6.2
>   [ERROR] r2_ins_data/0_roadmap_over_37/BREAST/snp_roadmap_119_enh_0.bed
>   5 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_BREAST_3.bed
>   6 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_CERVIX_2.bed
>   7 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_ESC_5.bed
>   8 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_ESC_DERIVED_9.bed
>   9 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_FAT_3.bed
>   10 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_COLON_7.bed
>   11 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_DUODENUM_5.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/GI_ESOPHAGUS/snp_roadmap_079_enh_0.bed
>   12 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_ESOPHAGUS_.bed
>   13 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_INTESTINE_6.bed
>   14 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_RECTUM_5.bed
>   15 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_STOMACH_5.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/HEART/snp_roadmap_105_enh_0.bed
>   16 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_HEART_2.bed
>   17 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_IPSC_4.bed
>   18 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_KIDNEY_1.bed
>   19 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_LIVER_7.bed
>   20 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_LUNG_6.bed
>   21 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_MUSCLE_3.bed
>   22 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_MUSCLE_LEG_2.bed
>   23 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_OVARY_1.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/PANCREAS/snp_roadmap_098_enh_0.bed
>   24 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_PANCREAS_2.bed
>   25 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_PLACENTA_8.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/0_roadmap_over_37/SKIN/snp_roadmap_057_enh_0.bed
>   26 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_SKIN_1.bed
>   27 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_SPLEEN_3.bed
>   28 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_STROMAL_CONNECTIVE_2.bed
>   29 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_THYMUS_13.bed
>   30 Write a BED file: r2_ins_data/summary/snp_union_0_roadmap_over_37_VASCULAR_4.bed

## Melton lab union list

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/1_meltonlab_over_37
  --out r2_ins_data/summary
  --sub_dir FALSE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 4 Files/folders input.
>   1 6 files in the beta_ATAC
>   2 6 files in the beta_enh
>   3 4 files in the other_ATAC
>   4 9 files in the other_enh
> Total 25 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
> In addition: Warning message:
> package 'eulerr' was built under R version 3.6.2
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_ATAC/snp_Beta_ATAC_rep1_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_DE_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_DE_TE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_SCbeta_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_SCbeta_TE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_Alpha_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_EN_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_PH_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_PH_TE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_PP_TE_enh_0.bed
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1]  8 19
>   Write a BED file: r2_ins_data/summary/snp_union_1_meltonlab_over_37_8.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-13 20:52:18 for 3.3 sec

### union list by groups

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/1_meltonlab_over_37
  --out r2_ins_data/summary
  --sub_dir TRUE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 4 Files/folders input.
>   1 sub_dir 1: 6 file(s) in the beta_ATAC folder
>   2 sub_dir 2: 6 file(s) in the beta_enh folder
>   3 sub_dir 3: 4 file(s) in the other_ATAC folder
>   4 sub_dir 4: 9 file(s) in the other_enh folder
> Total 4 sub-folder(s) is/are input
> Total 0 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
> In addition: Warning message:
> package 'eulerr' was built under R version 3.6.2
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_ATAC/snp_Beta_ATAC_rep1_0.bed
>   1 Write a BED file: r2_ins_data/summary/snp_union_1_meltonlab_over_37_beta_ATAC_4.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_DE_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_DE_TE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_SCbeta_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/beta_enh/snp_SCbeta_TE_enh_0.bed
>   2 Write a BED file: r2_ins_data/summary/snp_union_1_meltonlab_over_37_beta_enh_4.bed
>   3 Write a BED file: r2_ins_data/summary/snp_union_1_meltonlab_over_37_other_ATAC_2.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_Alpha_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_EN_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_PH_SE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_PH_TE_enh_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/1_meltonlab_over_37/other_enh/snp_PP_TE_enh_0.bed
>   4 Write a BED file: r2_ins_data/summary/snp_union_1_meltonlab_over_37_other_enh_5.bed

### Suh lab union list

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/1a_suhlab_over_37
  --out r2_ins_data/summary
  --sub_dir FALSE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 3 Files/folders input.
>   1 r2_ins_data/1a_suhlab_over_37/snp_Scbeta_YZ_pr1_pooled.pf_peaks_2.bed
>   2 r2_ins_data/1a_suhlab_over_37/snp_Scbeta_YZ_pr2_pooled.pf_peaks_1.bed
>   3 r2_ins_data/1a_suhlab_over_37/snp_SCbeta_YZ_pseudo_pool_ATAC_2.bed
> Total 3 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Scbeta_YZ_pr1_pooled.pf_peaks_2
>   Read 2: snp_Scbeta_YZ_pr2_pooled.pf_peaks_1
>   Read 3: snp_SCbeta_YZ_pseudo_pool_ATAC_2
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 2 7
>   Write a BED file: r2_ins_data/summary/snp_union_1a_suhlab_over_37_2.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-13 21:01:07 for 2.7 sec

## Pritchard lab union list

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/2_pritchardlab_over_37
  --out r2_ins_data/summary
  --sub_dir FALSE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 15 Files/folders input.
>   1 12 files in the B-cells_rest
>   2 10 files in the B-cells_stim
>   3 4 files in the Gamma_delta_T_rest
>   4 3 files in the Gamma_delta_T_stim
>   5 3 files in the Monocytes_rest
>   6 6 files in the Monocytes_stim
>   7 3 files in the Myeloid_DCs_rest
>   8 15 files in the NK-cells_rest
>   9 6 files in the NK-cells_stim
>   10 3 files in the pDCs_rest
>   11 3 files in the Plasmablasts_rest
>   12 38 files in the T-CD4-cells_rest
>   13 38 files in the T-CD4-cells_stim
>   14 16 files in the T-CD8-cells_rest
>   15 15 files in the T-CD8-cells_stim
> Total 175 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 175 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1]  17 179
>   Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_17.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-13 21:30:22 for 3.4 sec

### union list by groups

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/2_pritchardlab_over_37
  --out r2_ins_data/summary
  --sub_dir TRUE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 15 Files/folders input.
>   1 sub_dir 1: 12 file(s) in the B-cells_rest folder
>   2 sub_dir 2: 10 file(s) in the B-cells_stim folder
>   3 sub_dir 3: 4 file(s) in the Gamma_delta_T_rest folder
>   4 sub_dir 4: 3 file(s) in the Gamma_delta_T_stim folder
>   5 sub_dir 5: 3 file(s) in the Monocytes_rest folder
>   6 sub_dir 6: 6 file(s) in the Monocytes_stim folder
>   7 sub_dir 7: 3 file(s) in the Myeloid_DCs_rest folder
>   8 sub_dir 8: 15 file(s) in the NK-cells_rest folder
>   9 sub_dir 9: 6 file(s) in the NK-cells_stim folder
>   10 sub_dir 10: 3 file(s) in the pDCs_rest folder
>   11 sub_dir 11: 3 file(s) in the Plasmablasts_rest folder
>   12 sub_dir 12: 38 file(s) in the T-CD4-cells_rest folder
>   13 sub_dir 13: 38 file(s) in the T-CD4-cells_stim folder
>   14 sub_dir 14: 16 file(s) in the T-CD8-cells_rest folder
>   15 sub_dir 15: 15 file(s) in the T-CD8-cells_stim folder
> Total 15 sub-folder(s) is/are input
> Total 0 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
>   1 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_B-cells_rest_17.bed
>   2 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_B-cells_stim_17.bed
>   3 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Gamma_delta_T_rest_17.bed
>   4 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Gamma_delta_T_stim_17.bed
>   5 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Monocytes_rest_17.bed
>   6 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Monocytes_stim_17.bed
>   7 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Myeloid_DCs_rest_17.bed
>   8 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_NK-cells_rest_17.bed
>   9 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_NK-cells_stim_17.bed
>   10 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_pDCs_rest_17.bed
>   11 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Plasmablasts_rest_17.bed
>   12 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD4-cells_rest_17.bed
>   13 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD4-cells_stim_17.bed
>   14 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD8-cells_rest_17.bed
>   15 Write a BED file: r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD8-cells_stim_17.bed

## Tan lab union list

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/3_tanlab_over_37
  --out r2_ins_data/summary
  --sub_dir FALSE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 2 Files/folders input.
>   1 3 files in the Th1_enh
>   2 3 files in the Treg_enh
> Total 6 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Th1_enh_4
>   Read 2: snp_Th1_enh_control_1
>   Read 3: snp_Th1_enh_t1d_4
>   Read 4: snp_Treg_enh_4
>   Read 5: snp_Treg_enh_control_1
>   Read 6: snp_Treg_enh_t1d_3
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1]  4 10
>   Write a BED file: r2_ins_data/summary/snp_union_3_tanlab_over_37_4.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-13 21:31:21 for 2.8 sec

### union list by groups

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/3_tanlab_over_37
  --out r2_ins_data/summary
  --sub_dir TRUE
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 2 Files/folders input.
>   1 sub_dir 1: 3 file(s) in the Th1_enh folder
>   2 sub_dir 2: 3 file(s) in the Treg_enh folder
> Total 2 sub-folder(s) is/are input
> Total 0 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
>   1 Write a BED file: r2_ins_data/summary/snp_union_3_tanlab_over_37_Th1_enh_4.bed
>   2 Write a BED file: r2_ins_data/summary/snp_union_3_tanlab_over_37_Treg_enh_4.bed

# 7. Generating summary table

```CMD
Rscript postgwas-exe.r
  --dbvenn summ
  --base r2_ins_data/summary
  --out r2_ins_data/
  --uni_save FALSE
  --ann_gwas r2_ins_data/gwas_biomart.tsv
  --ann_encd r2_ins_data/distance_37/encode_tfbs.tsv
  --ann_near r2_ins_data/distance_37/nearest_gene.tsv
  --ann_cds r2_ins_data/distance_37/cds_gene.tsv
```

> ** Run function: db_venn.r/summ...
> 58 Files/folders input.
>   1 r2_ins_data/summary/snp_union_0_roadmap_over_37_23.bed
>   2 r2_ins_data/summary/snp_union_0_roadmap_over_37_ADRENAL_2.bed
>   3 r2_ins_data/summary/snp_union_0_roadmap_over_37_BLOOD_21.bed
>   4 r2_ins_data/summary/snp_union_0_roadmap_over_37_BONE_2.bed
>   5 r2_ins_data/summary/snp_union_0_roadmap_over_37_BRAIN_4.bed
>   6 r2_ins_data/summary/snp_union_0_roadmap_over_37_BREAST_3.bed
>   7 r2_ins_data/summary/snp_union_0_roadmap_over_37_CERVIX_2.bed
>   8 r2_ins_data/summary/snp_union_0_roadmap_over_37_ESC_5.bed
>   9 r2_ins_data/summary/snp_union_0_roadmap_over_37_ESC_DERIVED_9.bed
>   10 r2_ins_data/summary/snp_union_0_roadmap_over_37_FAT_3.bed
>   11 r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_COLON_7.bed
>   12 r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_DUODENUM_5.bed
>   13 r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_ESOPHAGUS_0.bed
>   14 r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_INTESTINE_6.bed
>   15 r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_RECTUM_5.bed
>   16 r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_STOMACH_5.bed
>   17 r2_ins_data/summary/snp_union_0_roadmap_over_37_HEART_2.bed
>   18 r2_ins_data/summary/snp_union_0_roadmap_over_37_IPSC_4.bed
>   19 r2_ins_data/summary/snp_union_0_roadmap_over_37_KIDNEY_1.bed
>   20 r2_ins_data/summary/snp_union_0_roadmap_over_37_LIVER_7.bed
>   21 r2_ins_data/summary/snp_union_0_roadmap_over_37_LUNG_6.bed
>   22 r2_ins_data/summary/snp_union_0_roadmap_over_37_MUSCLE_3.bed
>   23 r2_ins_data/summary/snp_union_0_roadmap_over_37_MUSCLE_LEG_2.bed
>   24 r2_ins_data/summary/snp_union_0_roadmap_over_37_OVARY_1.bed
>   25 r2_ins_data/summary/snp_union_0_roadmap_over_37_PANCREAS_2.bed
>   26 r2_ins_data/summary/snp_union_0_roadmap_over_37_PLACENTA_8.bed
>   27 r2_ins_data/summary/snp_union_0_roadmap_over_37_SKIN_1.bed
>   28 r2_ins_data/summary/snp_union_0_roadmap_over_37_SPLEEN_3.bed
>   29 r2_ins_data/summary/snp_union_0_roadmap_over_37_STROMAL_CONNECTIVE_2.bed
>   30 r2_ins_data/summary/snp_union_0_roadmap_over_37_THYMUS_13.bed
>   31 r2_ins_data/summary/snp_union_0_roadmap_over_37_VASCULAR_4.bed
>   32 r2_ins_data/summary/snp_union_1_meltonlab_over_37_8.bed
>   33 r2_ins_data/summary/snp_union_1_meltonlab_over_37_beta_ATAC_0.bed
>   34 r2_ins_data/summary/snp_union_1_meltonlab_over_37_beta_ATAC_4.bed
>   35 r2_ins_data/summary/snp_union_1_meltonlab_over_37_beta_enh_4.bed
>   36 r2_ins_data/summary/snp_union_1_meltonlab_over_37_other_ATAC_2.bed
>   37 r2_ins_data/summary/snp_union_1_meltonlab_over_37_other_enh_0.bed
>   38 r2_ins_data/summary/snp_union_1_meltonlab_over_37_other_enh_5.bed
>   39 r2_ins_data/summary/snp_union_1a_suhlab_over_37_2.bed
>   40 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_17.bed
>   41 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_B-cells_rest_17.bed
>   42 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_B-cells_stim_17.bed
>   43 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Gamma_delta_T_rest_17.bed
>   44 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Gamma_delta_T_stim_17.bed
>   45 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Monocytes_rest_17.bed
>   46 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Monocytes_stim_17.bed
>   47 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Myeloid_DCs_rest_17.bed
>   48 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_NK-cells_rest_17.bed
>   49 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_NK-cells_stim_17.bed
>   50 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_pDCs_rest_17.bed
>   51 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_Plasmablasts_rest_17.bed
>   52 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD4-cells_rest_17.bed
>   53 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD4-cells_stim_17.bed
>   54 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD8-cells_rest_17.bed
>   55 r2_ins_data/summary/snp_union_2_pritchardlab_over_37_T-CD8-cells_stim_17.bed
>   56 r2_ins_data/summary/snp_union_3_tanlab_over_37_4.bed
>   57 r2_ins_data/summary/snp_union_3_tanlab_over_37_Th1_enh_4.bed
>   58 r2_ins_data/summary/snp_union_3_tanlab_over_37_Treg_enh_4.bed
> Total 58 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
> In addition: Warning message:
> package 'eulerr' was built under R version 3.6.2
>   [ERROR] r2_ins_data/summary/snp_union_0_roadmap_over_37_GI_ESOPHAGUS_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/summary/snp_union_1_meltonlab_over_37_beta_ATAC_0.bed
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
>   [ERROR] r2_ins_data/summary/snp_union_1_meltonlab_over_37_other_enh_0.bed
>   Read 58 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 30 59
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 37 11
>   Merge dim     = [1] 37 63
>   Write a CSV file: r2_ins_data//summary_gwas.csv
>
>   ENCODE dim    = [1] 146  13
>   Merge dim     = [1] 121  58
>   Write a CSV file: r2_ins_data//summary_encode.csv
>
>   Nearest gene dim      = [1] 40  9
>   Search biomaRt... 7.. 6.. [1] 40  5
>   CDS dim               = [1] 166  11
>   Merge dim             = [1] 42 63
>   Write a CSV file: r2_ins_data//summary_nearest.csv
>
> Job done: 2020-03-13 22:23:43 for 25.3 sec