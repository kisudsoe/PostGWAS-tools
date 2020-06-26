# Type 1 Diabetes (r2 >0.6)

This is the result from the LDlink threshold as r<sup>2</sup> >0.6.

![](r2_data/Fig 1.png)



# 1. From the seed SNPs

## Filtering the LDlink data

Options: --r2d

* 1 = r2>0.6 and D'=1
* 2 = r2>0.6
* 3 = D'=1
* 4 = r2>0.6 or D'=1

```CMD
Rscript postgwas-exe.r ^
  --ldlink 	filter ^
  --base 	db_gwas/gwas_5e-08_152.tsv db_gwas/ldlink ^
  --out 	r2_data ^
  --r2d 	2
```

> ** Run function ldlink_filter...
> Read download files... 152
> Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
>   line 1 did not have 11 elements
>   [ERROR] rs75793288
>   Read LDlink results           = [1] 303054     12
> Filtering by "r2 > 0.6":
>   Filtered data dimension       = [1] 7630    3
>   Excluded no rsid elements     = [1] 117
> Basic summary of LDlink results:
>   **SNP Tier 1                    = 152**
>   **SNP Tier 2                    = 5738**
>   **SNP candidates                = 5890**
>   SNP source annotation table   = [1] 5890    2
> Add annotations:
>   LD block annotation... [1] 119
> Search biomart for SNP coordinates:
>   Query SNPs            = [1] 5890
>   Hg19 result table     = [1] 5835    4
>   Hg38 result table     = [1] 5886    4
>   Cytoband annotation... 5891.. done
>   Merged table          = [1] 5891   11
> Write file: r2_data/gwas_biomart.tsv
> Job done: 2020-03-10 18:59:28 for 38.8 sec

## Generating BED files for hg19 coordinates

File name change: `r2_data/gwas_biomart.tsv` -> `r2_data/gwas_biomart.fill.tsv`

```CMD
Rscript postgwas-exe.r ^
  --ldlink 	bed ^
  --base 	r2_data/gwas_biomart_fill.tsv ^
  --out 	r2_data	
```

> ** Run function ldlink_bed... [1] 5891   11
> Write file:     r2_data/gwas_hg19_biomart_5890.bed
> Write file:     r2_data/gwas_hg38_biomart_5886.bed
> Job done: 2020-04-09 19:36:41 for 3.9 sec

# 2. Downloading annotation data

Roadmap, ENCODE, RegulomeDB, GTEx, lncRNASNP2, and Ensembl Gene annotation were already downloaded in `r2d1_report.md`.

## ATAC data compare

1. Generating merged bed file

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 bedmerge ^
  --base	 db_1_meltonlab/atac_beta db_1a_Suhlab/atac db_1b_Stitzellab/atac db_2_Pritchardlab/bed_ATAC_T_mono db_4_Ferrerlab/atac ^
  --out		 pca/atac_merge_input.bed
```

> ** Run function: db_compare.r/bedmerge...
> 1 db_1_meltonlab/atac_beta: [1] 6
>   db_1_meltonlab/atac_beta/GSM4171638_SCbeta_ATAC_rep1.peaks.bed
>    db_1_meltonlab/atac_beta/GSM4171640_Beta_ATAC_rep1.peaks.bed
>    db_1_meltonlab/atac_beta/GSM4171642_in_vivo-matured_SCbeta_ATAC_rep1.peaks.bed
>    db_1_meltonlab/atac_beta/GSM4171643_in_vivo-matured_SCbeta_ATAC_rep2.peaks.bed
>    db_1_meltonlab/atac_beta/GSM4171644_in_vitro-matured_SCbeta_ATAC_12h.peaks.bed
>    db_1_meltonlab/atac_beta/GSM4171645_in_vitro-matured_SCbeta_ATAC_72h.peaks.bed
> 2 db_1a_Suhlab/atac: [1] 2
>   db_1a_Suhlab/atac/YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled.pf_peaks.bed
>    db_1a_Suhlab/atac/YY005-beta-cell-1_R1.nodup.tn5.pr2_pooled.pf_peaks.bed
> 3 db_1b_Stitzellab/atac: [1] 1
>   db_1b_Stitzellab/atac/GSM3333912_EndoC_BH1_ATACseq_broadPeak.fdr0.05.noBlacklist.bed
> 4 db_2_Pritchardlab/bed_ATAC_T_mono: [1] 123
> 5 db_4_Ferrerlab/atac: [1] 1
>   db_4_Ferrerlab/atac/atac_consistent_peaks.bed
> ** 133 BED file paths are ready.
>
>   Data merge...
>     20/133 process done.
>     40/133 process done.
>     60/133 process done.
>     80/133 process done.
>     100/133 process done.
>     120/133 process done. 133..   done: [1] 98083348        4
>   Write file: pca/atac_merge_input.bed
> Job done: 2020-04-29 21:05:02 for 16.7 min

2. Using bedtools merge

```bash
sortBed -i pca/atac_merge_input.bed | mergeBed -i stdin > pca/atac_merge.bed
```

3. Using bedtools distance

```bash
bash bash_sh/dbcomp_atac_dist.sh
```

4. Get union table for PCA analysis

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 uniontb ^
  --base	 pca/atac_dist ^
  --out		 pca/atac_uniontb.csv
```

> ** Run function: db_compare.r/uniontb...
> 1 pca/atac_dist: [1] 135
> ** 135 dist file paths are ready.
>   
>      Read data...
>        20/135 process done.
>        40/135 process done.
>        60/135 process done.
>        80/135 process done.
>        100/135 process done.
>        120/135 process done. 135..
>      Data union... [1] 1038289     136
>      Write file: pca/atac_uniontb.csv
>   Job done: 2020-04-29 22:14:08 for 21.4 min

5. Draw PCA plot

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 pcaplot ^
  --base	 pca/atac_uniontb.csv ^
  --out		 pca ^
  --meta	 pca/atac_meta.tsv ^
  --filt	 Melton Suh Stitzel Pritchard Ferrer
```

> ** Run function: db_compare.r/pcaplot...
>   Read file: pca/atac_uniontb.csv... [1] 1038289     136
>   Convert T/F to 1/0... done
>   Read meta file: pca/atac_meta.tsv... [1] 133   3
>   Filter option: Melton, Suh, Stitzel, Pritchard, Ferrer >> [1] 133
>   PCA analysis... done
>   Write file: pca/atac_uniontb_pca.csv
>   Draw PCA plot... hull.. done
> Save plot: pca/atac_uniontb_pca.png
> Job done: 2020-04-29 22:18:21 for 3.4 min

Change file name: `atac_uniontb_pca.png` -> `atac_uniontb-all_pca.png`

![](pca\atac_uniontb-all_pca.png)

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 pcaplot ^
  --base	 pca/atac_uniontb.csv ^
  --out		 pca ^
  --meta	 pca/atac_meta.tsv ^
  --filt	 Melton Stitzel Pritchard Ferrer
```

> ** Run function: db_compare.r/pcaplot...
>     Read file: pca/atac_uniontb.csv... [1] 1038289     136
>     Convert T/F to 1/0... done
>     Read meta file: pca/atac_meta.tsv... [1] 133   3
>     Filter option: Melton, Stitzel, Pritchard, Ferrer >> [1] 131
>     PCA analysis... done
>     Write file: pca/atac_uniontb_pca.csv
>     Draw PCA plot... hull.. done
> Save plot: pca/atac_uniontb_pca.png
> Job done: 2020-04-29 22:29:51 for 3.3 min

File name change: `atac_uniontb_pca.png` -> `atac_uniontb-suh_pca.png`

![](pca\atac_uniontb-suh_pca.png)



## ChIP-seq Enhancer data compare

1.  Generating merged bed file

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 bedmerge ^
  --base	 db_1_meltonlab/enh_beta db_1b_Stitzellab/enh db_3_Tanlab/enh_bed_files db_4_Ferrerlab/enh_active db_gwas/roadmap_enh_panc,T,mono ^
  --out		 pca/enh_merge_input.bed
```

> ** Run function: db_compare.r/bedmerge...
> 1 db_1_meltonlab/enh_beta: [1] 4
>   db_1_meltonlab/enh_beta/GSE140500_Beta_SE.byH3K27ac.bed
>      db_1_meltonlab/enh_beta/GSE140500_Beta_TE.byH3K27ac.bed
>    db_1_meltonlab/enh_beta/GSE140500_SCbeta_SE.byH3K27ac.bed
>      db_1_meltonlab/enh_beta/GSE140500_SCbeta_TE.byH3K27ac.bed
>    2 db_1b_Stitzellab/enh: [1] 1
>      db_1b_Stitzellab/enh/EndoC_BH1_enh.bed
>    3 db_3_Tanlab/enh_bed_files: [1] 6
>      db_3_Tanlab/enh_bed_files/Th1_enh.bed
>       db_3_Tanlab/enh_bed_files/Th1_enh_control.bed
>    db_3_Tanlab/enh_bed_files/Th1_enh_t1d.bed
>    db_3_Tanlab/enh_bed_files/Treg_enh.bed
>   db_3_Tanlab/enh_bed_files/Treg_enh_control.bed
>      db_3_Tanlab/enh_bed_files/Treg_enh_t1d.bed
>    4 db_4_Ferrerlab/enh_active: [1] 1
>      db_4_Ferrerlab/enh_active/islet_regulome_active_enh.bed
>    5 db_gwas/roadmap_enh_panc,T,mono: [1] 20
>    ** 32 BED file paths are ready.
>    
>      Data merge...
>        20/32 process done. 32..   done: [1] 2965364       4
>     Write file: pca/enh_merge_input.bed
> Job done: 2020-04-30 11:05:14 for 13.3 sec

2. Using bed tools merge

```bash
sortBed -i pca/enh_merge_input.bed | mergeBed -i stdin > pca/enh_merge.bed
```

3. Using bedtools distance

```bash
bash bash_sh/dbcomp_enh_dist.sh
```

4. Get union table for PCA analysis

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 uniontb ^
  --base	 pca/enh_dist ^
  --out		 pca/enh_uniontb.csv
```

> ** Run function: db_compare.r/uniontb...
> 1 pca/enh_dist: [1] 30
> ** 30 dist file paths are ready.
>
>     Read data...
>        10/30 process done.
>        20/30 process done.
>        30/30 process done. 30..
>      Data union... [1] 291127     31
>      Write file: pca/enh_uniontb.csv
>    Job done: 2020-04-30 11:09:29 for 23.2 sec

5. Draw PCA plot

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 pcaplot ^
  --base	 pca/enh_uniontb.csv ^
  --out		 pca ^
  --meta	 pca/enh_meta.tsv ^
  --filt	 roadmap_immn roadmap_panc Melton Suh Stitzel Tan Ferrer
```

> ** Run function: db_compare.r/pcaplot...
>     Read file: pca/enh_uniontb.csv... [1] 291127     31
>     Convert T/F to 1/0... done
>     Read meta file: pca/enh_meta.tsv... [1] 30  3
>     Filter option: roadmap, Melton, Suh, Stitzel, Tan, Ferrer >> [1] 30
>     PCA analysis... done
>     Write file: pca/enh_uniontb_pca.csv
>     Draw PCA plot... hull.. done
> Save plot: pca/enh_uniontb_pca.png
> Job done: 2020-04-30 11:23:19 for 29.1 sec

File name change: `enh_uniontb_pca.png` -> `enh_uniontb-all_pca.png`

![](pca\enh_uniontb-all_pca.png)

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 pcaplot ^
  --base	 pca/enh_uniontb.csv ^
  --out		 pca ^
  --meta	 pca/enh_meta.tsv ^
  --filt	 roadmap_immn Suh Tan
```

> ** Run function: db_compare.r/pcaplot...
>   Read file: pca/enh_uniontb.csv... [1] 291127     31
>   Convert T/F to 1/0... done
>   Read meta file: pca/enh_meta.tsv... [1] 30  3
>   Filter option: roadmap_immn, Suh, Tan >> [1] 22
>   PCA analysis... done
>   Write file: pca/enh_uniontb_pca.csv
>   Draw PCA plot... hull.. done
> Save plot: pca/enh_uniontb_pca.png
> Job done: 2020-04-30 11:25:52 for 27.8 sec

File name change: `enh_uniontb_pca.csv` -> `enh_uniontb-immn_pca.csv`

## Hi-C data compare

1. Generating merged bed file

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 bedmerge ^
  --base	 db_1a_Suhlab/hic_bed_lftovr_hg19 db_1b_Stitzellab/hic db_3a_Allanlab/hic_bed db_3b_Aidenlab/hic_bed db_4_Ferrerlab/hic ^
  --out		 pca/hic_merge_input.bed
```

> ** Run function: db_compare.r/bedmerge...
> 1 db_1a_Suhlab/hic_bed_lftovr_hg19: [1] 8
>   db_1a_Suhlab/hic_bed_lftovr_hg19/hic_mono_1018_hiccups_merged_loops_hg19.bed
>    db_1a_Suhlab/hic_bed_lftovr_hg19/hic_mono_1027_hiccups_merged_loops_hg19.bed
>    db_1a_Suhlab/hic_bed_lftovr_hg19/hic_mono_1035_hiccups_merged_loops_hg19.bed
>    db_1a_Suhlab/hic_bed_lftovr_hg19/hic_mono_1046_hiccups_merged_loops_hg19.bed
>    db_1a_Suhlab/hic_bed_lftovr_hg19/hic_T_1018_hiccups_merged_loops_hg19.bed
>    db_1a_Suhlab/hic_bed_lftovr_hg19/hic_T_1027_hiccups_merged_loops_hg19.bed
>    db_1a_Suhlab/hic_bed_lftovr_hg19/hic_T_1035_hiccups_merged_loops_hg19.bed
>    db_1a_Suhlab/hic_bed_lftovr_hg19/hic_T_1046_hiccups_merged_loops_hg19.bed
> 2 db_1b_Stitzellab/hic: [1] 2
>   db_1b_Stitzellab/hic/GSM3333898_Human_Islet_HiC_HiCCUPS_loops.bed
>    db_1b_Stitzellab/hic/GSM3333916_EndoC_BH1_HiC_HiCCUPS_loops.bed
> 3 db_3a_Allanlab/hic_bed: [1] 6
>   db_3a_Allanlab/hic_bed/GSM2827786_CD4T1.bed
>    db_3a_Allanlab/hic_bed/GSM2827787_CD4T2.bed
>    db_3a_Allanlab/hic_bed/GSM2827788_CD8T1.bed
>    db_3a_Allanlab/hic_bed/GSM2827789_CD8T2.bed
>    db_3a_Allanlab/hic_bed/GSM2827790_HB1.bed
>    db_3a_Allanlab/hic_bed/GSM2827791_HB2.bed
> 4 db_3b_Aidenlab/hic_bed: [1] 3
>   db_3b_Aidenlab/hic_bed/GSE63525_GM12878_primary_HiCCUPS_looplist.bed
>    db_3b_Aidenlab/hic_bed/GSE63525_GM12878_replicate_HiCCUPS_looplist.bed
>    db_3b_Aidenlab/hic_bed/GSE63525_K562_HiCCUPS_looplist.bed
> 5 db_4_Ferrerlab/hic: [1] 1
>   db_4_Ferrerlab/hic/pi_merged_washu_text.bed
> ** 20 BED file paths are ready.
>
>   Data merge... 20..   done: [1] 1465395       4
>   Write file: pca/hic_merge_input.bed
> Job done: 2020-04-29 19:12:46 for 14.1 sec

2. Using bed tools merge

```bash
sortBed -i pca/hic_merge_input.bed | mergeBed -i stdin > pca/hic_merge.bed
```

3. Using bedtools distance

```bash
bash bash_sh/dbcomp_hic_dist.sh
```

4. Get union table for PCA analysis

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 uniontb ^
  --base	 pca/hic_dist ^
  --out		 pca/hic_uniontb.csv ^
  --hic		 TRUE
```

> ** Run function: db_compare.r/uniontb...
> 1 pca/hic_dist: [1] 20
> ** 20 dist file paths are ready.
>
>   Read data...
>     5/20 process done.
>     10/20 process done.
>     15/20 process done.
>     20/20 process done. 20..
>   Data union... [1] 13631    21
>   Write file: pca/hic_uniontb.csv
> Job done: 2020-04-29 19:52:42 for 28 min

5. Draw PCA plot

```CMD
Rscript postgwas-exe.r ^
  --dbcomp	 pcaplot ^
  --base	 pca/hic_uniontb.csv ^
  --out		 pca ^
  --meta	 pca/hic_meta.tsv ^
  --filt	 Suh Stitzel Allan Aiden Ferrer
```

> ** Run function: db_compare.r/pcaplot...
>   Read file: pca/hic_uniontb.csv... [1] 13631    21
>   Convert T/F to 1/0... done
>   Read meta file: pca/hic_meta.tsv... [1] 20  3
>   Filter option: Suh, Stitzel, Allan, Aiden, Ferrer >> [1] 18
>   PCA analysis... done
>   Write file: pca/hic_uniontb_pca.csv
>   Draw PCA plot... hull.. done
> Save plot: pca/hic_uniontb_pca.png
> Job done: 2020-04-29 20:07:19 for 2.8 sec

File name change: `hic_uniontb_pca.png` -> `hic_uniontb-all_pca.png`

![](pca\hic_uniontb-all_pca.png)

```CMD
# WITHOUT beta/islet related cells
Rscript postgwas-exe.r ^
  --dbcomp	 pcaplot ^
  --base	 pca/hic_uniontb.csv ^
  --out		 pca ^
  --meta	 pca/hic_meta.tsv ^
  --filt	 Suh Allan Aiden
```

> ** Run function: db_compare.r/pcaplot...
>   Read file: pca/hic_uniontb.csv... [1] 13631    21
>   Convert T/F to 1/0... done
>   Read meta file: pca/hic_meta.tsv... [1] 20  3
>   Filter option: Suh, Allan, Aiden >> [1] 15
>   PCA analysis... done
>   Write file: pca/hic_uniontb_pca.csv
>   Draw PCA plot... hull.. done
> Save plot: pca/hic_uniontb_pca.png
> Job done: 2020-04-30 10:47:56 for 2.5 sec

File name change: `hic_uniontb_pca.png` -> ``

# 3. Filtering/converting the annotations

UCSC annotations, Roadmap, GTEx, Stitzel lab's Hi-C data were already filtered in `r2d1_report.md`.

## Roadmap enhancer merge/closest

```CMD
Rscript postgwas-exe.r ^
  --dbfilt	 roadmap ^
  --base	 db_gwas/roadmap ^
  --out		 db_gwas
```

> ** Run function: db_filter.r/roadmap_filt...
>   Reading files..
>     10/129 being processed.
>     20/129 being processed.
>     30/129 being processed.
>     40/129 being processed.
>     50/129 being processed.
> Error in gzfile(file, "rb") : cannot open the connection
> In addition: Warning message:
> In gzfile(file, "rb") :
>   cannot open compressed file 'db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
>   db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds - file not found.
>     60/129 being processed.
> Error in gzfile(file, "rb") : cannot open the connection
> In addition: Warning message:
> In gzfile(file, "rb") :
>   cannot open compressed file 'db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
>   db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds - file not found.
>     70/129 being processed.
>     80/129 being processed.
>     90/129 being processed.
>     100/129 being processed.
>     110/129 being processed.
>     120/129 being processed.
>   Finished reading and filtering 129 files.
>
> Write file: db_gwas/roadmap_enh.bed
> Job done: 2020-04-11 19:17:57 for 6 min

## GTEx filter

Filtering GTEx eQTL data by p-value < 5e-8 is already done. See details in `r2d1_report.md`.

## Converting additional data

### Stitzel lab's Hi-C data

This step is already done. See details in `r2d1_report.md`.

### Ferror lab's Enhancer, pcHi-C data

This step is simply done by excel. See details for pcHi-C data conversion in `r2d1_report.md`.



# 4. Distances from the annotations

## For general annotations: CDS_gene, encode_tfbs, nearest_gene

```bash
bash distance_r2.sh
```

## For roadmap annotations

```bash
bash 0_roadmap_dist_r2.sh
```

```CMD
# generate and distance bedtools merge file
bedtools sort -i db_gwas/roadmap_enh.bed | bedtools merge -i stdin -c 1 -o count > db_gwas/roadmap_enh_merge.bed
bedtools sort -i r2_data/gwas_hg19_biomart_5890.bed | bedtools closest -d -a stdin -b db_gwas/roadmap_enh_merge.bed > r2_data/distance_r2/roadmap_enh_merge.tsv
```

## Preparing additional data

Using AWS server, see details in `scp.txt` file.

bash files are located in `bash_sh` folder.

### Melton lab's β-cell ATAC-seq, enhancer data

```bash
bash 1_meltonlab_dist_r2.sh
```

### Suh lab Yizhou's β-cell ATAC-seq, Tri-Hic data

```bash
bash 1a_suhlab_dist_r2.sh
```

```bash
bash 1a_suhlab_dist_hic_r2.sh
```

### Stitzel lab's EndoC-βH1 cell ATAC-seq, enhancer, Hi-C data

```bash
bash 1b_stitzellab_dist_r2.sh
```

### Pritchard lab's immune cell ATAC-seq data

```bash
bash 2_pritchardlab_dist_r2.sh
```

### Tan lab's Th1 and Treg cell enhancer data

```bash
bash 3_3_tan_dist_r2.sh
```

### Allan lab's CD4 T, CD8 T, and B cell Hi-C data

```bash
bash 3a_allanlab_dist_r2.sh
```

### Aiden lab's GM12878, K562 cells Hi-C data

```bash
bash 3b_aidenlab_dist_hic_r2.sh
```

### Ferrer lab's islet ATAC, enhancer, pcHi-C data

```bash
bash 4_ferrerlab_dist_r2.sh
```



# 5. Overlapping the annotations

## UCSC gene annotations

```CMD
Rscript postgwas-exe.r ^
  --dbfilt 		dist ^
  --base 		r2_data/distance_r2/ucsc_annot.tsv ^
  --out 		r2_data/summary_r2 ^
  --infotype 	ucsc
```

> Input file/folder N     = [1] 1
>
> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 1
> File ucsc_annot... nrow= 23459.. done
>   Annotations occupied by SNPs  = [1] 2668
>   SNPs in annotations           = [1] 3811
>
>   UCSC annotations: 0 tags
>     1 NA:       0..     Save at: r2_data/summary/snp_ucsc_NA_0.bed
>     0 : 0..     Save at: r2_data/summary/snp_ucsc__0.bed
> Job done: 2020-04-10 00:59:10 for 0.8 sec

## Roadmap each cell type

```CMD
Rscript postgwas-exe.r ^
  --dbfilt	dist ^
  --base	r2_data/0_roadmap_dist_r2 ^
  --out 	r2_data/0_roadmap_over_r2 ^
  --meta 	db_gwas/roadmap_meta.tsv
```

> ...
>
> File roadmap_128_enh... nrow= 5891.. done
>   Annotations occupied by SNPs  = [1] 184
>   SNPs in annotations           = [1] 240
>   Write file: r2_data/0_roadmap_over_r2/LUNG_encode/snp_roadmap_128_enh_240.bed
> Job process: 0.1 sec
>
> File roadmap_129_enh... nrow= 5891.. done
>   Annotations occupied by SNPs  = [1] 189
>   SNPs in annotations           = [1] 240
>   Write file: r2_data/0_roadmap_over_r2/BONE_encode/snp_roadmap_129_enh_240.bed
> Job process: 0.1 sec
>
> Job done: 2020-04-10 01:11:36 for 12.8 sec

## RegulomeDB annotations

```CMD
Rscript postgwas-exe.r ^
  --dbfilt	regulome ^
  --base	r2_data/gwas_hg19_biomart_5890.bed ^
  --regulm 	db_gwas/regulome ^
  --out 	r2_data/summary_r2
```

> ** Run function: db_filter.r/regulome_filt...
> Input GWAS SNPs N       = [1] 5890
> 2 Regulome data load...
>   Read: db_gwas/regulome/dbSNP132.Category1.txt.gz.rds; dim = [1] 39432     5
>   Read: db_gwas/regulome/dbSNP132.Category2.txt.gz.rds; dim = [1] 407796      5
>
>   Regulome score >=2b, SNPs             = [1] 430528
>   Functional motifs (1a~2b - 1f only)   = [1] 34705
>
>   Regulome >=2b, GWAS SNPs              = [1] 332
>   GWAS SNPs occupied in
>     functional motifs (1a~2b - 1f only) = [1] 177
>
> Write file: r2_data/summary/regulome_332.tsv
> Write file: r2_data/summary/snp_regulome2b_332.bed
>
> Job done: 2020-04-10 01:15:01 for 8.3 sec

Move file: `data/summary/regulome_332.tsv` -> `data/distance_r2/regulome_332.tsv`

## GTEx eQTL data

```CMD
Rscript postgwas-exe.r ^
  --dbfilt gtex_ovl ^
  --base r2_data/gwas_hg19_biomart_5890.bed ^
  --gtex db_gwas/gtex_signif_5e-08.rds ^
  --out r2_data/gtex_eqtl
```

> ** Run function: db_filter.r/gtex_overlap...
> Input GWAS SNPs N       = 5890
>   gtex_signif_5e-08.rds, dim    = [1] 30613850        9
>   Overlapped eQTL-gene pairs    = [1] 270976
>   eQTLs N               = [1] 3921
>   Associated eGenes     = [1] 286
>
> Write file: r2_data/gtex_eqtl/gtex_signif_3921.tsv
>
> Generating BED files for 49 tissues.. done
>
> Job done: 2020-04-21 10:47:30 for 2.3 min

Move file: `r2_data/gtex_eqtl/gtex_signif_3921.tsv` -> `r2_data/gtex_signif_3921.tsv`

## lncRNASNP2 data

```CMD
Rscript postgwas-exe.r ^
  --dbfilt lnc_ovl ^
  --base r2_data/gwas_hg19_biomart_5890.bed ^
  --lncrna db_gwas/lncrna ^
  --out r2_data
```

> ** Run function: db_filter.r/lncrna_overlap...
> Input GWAS SNPs N = 5890
> 3 lncRNASNP2 data load...
>   Read: db_gwas/lncrna/lncRNASNP2_snplist.txt.rds;              dim = [1] 10205295        3
>   Read: db_gwas/lncrna/lncrnas.txt.rds;                         dim = [1] 141271      4
>   Read: db_gwas/lncrna/lncrna-diseases_experiment.txt.rds;      dim = [1] 753   3
>
> Summary =
>   lncRNA SNPs
> 1    166  316
>
>   Write file: r2_data/snp_lncrnasnp_316.bed
>   Write file: r2_data/lncrnasnp_316.tsv
> Job done: 2020-04-21 10:49:23 for 28.6 sec

Move file: `r2_data/snp_lncrnasnp_316.bed` -> `r2_data/summary_r2/snp_lncrnasnp_316.bed`

## Additional data

### Melton lab's data - ATAC, Enhancer

```CMD
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/1_meltonlab_dist_r2 ^
  --out r2_data/1_meltonlab_over_r2 ^
  --meta db_gwas/meltonlab_meta.tsv
```

> Input file/folder N     = [1] 1
>
> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 25
>   Read metadata file dim        = [1] 25  2
> File Alpha_ATAC_rep1... nrow= 5890.. done
>   Annotations occupied by SNPs  = [1] 31
>   SNPs in annotations           = [1] 41
>   Write file: r2_data/1_meltonlab_over_r2/other_ATAC/snp_Alpha_ATAC_rep1_41.bed
> Job process: 0.3 sec
>
> ...
>
> File SCbeta_TE_enh... nrow= 5890.. done
>   Annotations occupied by SNPs  = [1] 17
>   SNPs in annotations           = [1] 37
>   Write file: r2_data/1_meltonlab_over_r2/beta_enh/snp_SCbeta_TE_enh_37.bed
> Job process: 0.2 sec
>
> Job done: 2020-04-21 10:52:18 for 3.2 sec

### Suh lab Yizhou's data - ATAC, Tri-HiC

```CMD
# ATAC-seq data
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/1a_suhlab_dist_r2
  --out r2_data/1a_suhlab_over_r2
```

> Input file/folder N     = [1] 1
>
> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 3
> File Scbeta_YZ_pr1_pooled.pf_peaks... nrow= 9306.. done
>   Annotations occupied by SNPs  = [1] 115
>   SNPs in annotations           = [1] 168
>   Write file: r2_data/1a_suhlab_over_r2/snp_Scbeta_YZ_pr1_pooled.pf_peaks_168.bed
>   1/3 done: r2_data/1a_suhlab_dist_r2/Scbeta_YZ_pr1_pooled.pf_peaks.tsv
> Job process: 0.3 sec
>
> File Scbeta_YZ_pr2_pooled.pf_peaks... nrow= 10166.. done
>   Annotations occupied by SNPs  = [1] 123
>   SNPs in annotations           = [1] 177
>   Write file: r2_data/1a_suhlab_over_r2/snp_Scbeta_YZ_pr2_pooled.pf_peaks_177.bed
>   2/3 done: r2_data/1a_suhlab_dist_r2/Scbeta_YZ_pr2_pooled.pf_peaks.tsv
> Job process: 0.2 sec
>
> File SCbeta_YZ_pseudo_pool_ATAC... nrow= 9743.. done
>   Annotations occupied by SNPs  = [1] 96
>   SNPs in annotations           = [1] 140
>   Write file: r2_data/1a_suhlab_over_r2/snp_SCbeta_YZ_pseudo_pool_ATAC_140.bed
>   3/3 done: r2_data/1a_suhlab_dist_r2/SCbeta_YZ_pseudo_pool_ATAC.tsv
> Job process: 0.2 sec
>
> Job done: 2020-04-21 10:53:37 for 0.7 sec

```CMD
# Tri-HiC data for mono_1018 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/mono_1018_gwas.tsv r2_data/1a_suhlab_dist_r2/mono_1018_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/mono_1018_gwas.tsv, length= 12968,  overlap= 787
>   r2_data/1a_suhlab_dist_r2/mono_1018_gene.tsv, length= 3187877,        overlap= 828970
>
>   Process gwas_loop.. 787.. done
>   Process gene_loop.. 828970.. done
>   Process merge for TSV.. 787.. [1] 15839     5
>   Write file: r2_data/summary_gene_r2/hic_mono_1018_gwas_1440.tsv
>
>   Process extract for BED.. 787.. [1] 288   4
>   Write file: r2_data/summary_r2/snp_hic_mono_1018_gwas_288.bed
> Job done: 2020-04-28 22:48:36 for 24.6 sec

```CMD
# Tri-HiC data for mono_1027 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/mono_1027_gwas.tsv r2_data/1a_suhlab_dist_r2/mono_1027_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/mono_1027_gwas.tsv, length= 12143,  overlap= 798
>   r2_data/1a_suhlab_dist_r2/mono_1027_gene.tsv, length= 3195753,        overlap= 832957
>
>   Process gwas_loop.. 798.. done
>   Process gene_loop.. 832957.. done
>   Process merge for TSV.. 798.. [1] 14104     5
>   Write file: r2_data/summary_gene_r2/hic_mono_1027_gwas_1345.tsv
>
>   Process extract for BED.. 798.. [1] 325   4
>   Write file: r2_data/summary_r2/snp_hic_mono_1027_gwas_325.bed
> Job done: 2020-04-28 22:50:15 for 24.8 sec

```CMD
# Tri-HiC data for mono_1035 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/mono_1035_gwas.tsv r2_data/1a_suhlab_dist_r2/mono_1035_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/mono_1035_gwas.tsv, length= 14052,  overlap= 937
>   r2_data/1a_suhlab_dist_r2/mono_1035_gene.tsv, length= 3656912,        overlap= 957084
>
>   Process gwas_loop.. 937.. done
>   Process gene_loop.. 957084.. done
>   Process merge for TSV.. 937.. [1] 19948     5
>   Write file: r2_data/summary_gene_r2/hic_mono_1035_gwas_1699.tsv
>
>   Process extract for BED.. 937.. [1] 311   4
>   Write file: r2_data/summary_r2/snp_hic_mono_1035_gwas_311.bed
> Job done: 2020-04-28 22:51:40 for 30.2 sec

```CMD
# Tri-HiC data for mono_1046 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/mono_1046_gwas.tsv r2_data/1a_suhlab_dist_r2/mono_1046_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/mono_1046_gwas.tsv, length= 14856,  overlap= 1047
>   r2_data/1a_suhlab_dist_r2/mono_1046_gene.tsv, length= 3831503,        overlap= 1048554
>
>   Process gwas_loop.. 1047.. done
>   Process gene_loop.. 1048554.. done
>   Process merge for TSV.. 1047.. [1] 20809     5
>   Write file: r2_data/summary_gene_r2/hic_mono_1046_gwas_1871.tsv
>
>   Process extract for BED.. 1047.. [1] 372   4
>   Write file: r2_data/summary_r2/snp_hic_mono_1046_gwas_372.bed
> Job done: 2020-04-28 22:52:47 for 33.9 sec

```CMD
# Tri-HiC data for T_1018 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/T_1018_gwas.tsv r2_data/1a_suhlab_dist_r2/T_1018_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/T_1018_gwas.tsv, length= 12942,     overlap= 528
>   r2_data/1a_suhlab_dist_r2/T_1018_gene.tsv, length= 2864793,   overlap= 536919
>
>   Process gwas_loop.. 528.. done
>   Process gene_loop.. 536919.. done
>   Process merge for TSV.. 528.. [1] 11233     5
>   Write file: r2_data/summary_gene_r2/hic_T_1018_gwas_1233.tsv
>
>   Process extract for BED.. 528.. [1] 222   4
>   Write file: r2_data/summary_r2/snp_hic_T_1018_gwas_222.bed
> Job done: 2020-04-28 22:53:39 for 15.7 sec

```CMD
# Tri-HiC data for T_1027 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/T_1027_gwas.tsv r2_data/1a_suhlab_dist_r2/T_1027_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/T_1027_gwas.tsv, length= 11986,     overlap= 738
>   r2_data/1a_suhlab_dist_r2/T_1027_gene.tsv, length= 3552381,   overlap= 755083
>
>   Process gwas_loop.. 738.. done
>   Process gene_loop.. 755083.. done
>   Process merge for TSV.. 738.. [1] 15000     5
>   Write file: r2_data/summary_gene_r2/hic_T_1027_gwas_1462.tsv
>
>   Process extract for BED.. 738.. [1] 244   4
>   Write file: r2_data/summary_r2/snp_hic_T_1027_gwas_244.bed
> Job done: 2020-04-28 22:54:19 for 22.4 sec

```CMD
# Tri-HiC data for T_1035 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/T_1035_gwas.tsv r2_data/1a_suhlab_dist_r2/T_1035_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/T_1035_gwas.tsv, length= 12129,     overlap= 822
>   r2_data/1a_suhlab_dist_r2/T_1035_gene.tsv, length= 3781085,   overlap= 798756
>
>   Process gwas_loop.. 822.. done
>   Process gene_loop.. 798756.. done
>   Process merge for TSV.. 822.. [1] 16740     5
>   Write file: r2_data/summary_gene_r2/hic_T_1035_gwas_1626.tsv
>
>   Process extract for BED.. 822.. [1] 258   4
>   Write file: r2_data/summary_r2/snp_hic_T_1035_gwas_258.bed
> Job done: 2020-04-28 22:59:22 for 24.8 sec

```CMD
# Tri-HiC data for T_1046 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1a_suhlab_dist_r2/T_1046_gwas.tsv r2_data/1a_suhlab_dist_r2/T_1046_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1a_suhlab_dist_r2/T_1046_gwas.tsv, length= 12267,     overlap= 682
>   r2_data/1a_suhlab_dist_r2/T_1046_gene.tsv, length= 3312327,   overlap= 667108
>
>   Process gwas_loop.. 682.. done
>   Process gene_loop.. 667108.. done
>   Process merge for TSV.. 682.. [1] 15410     5
>   Write file: r2_data/summary_gene_r2/hic_T_1046_gwas_1522.tsv
>
>   Process extract for BED.. 682.. [1] 244   4
>   Write file: r2_data/summary_r2/snp_hic_T_1046_gwas_244.bed
> Job done: 2020-04-28 22:59:55 for 19.8 sec

### Stitzel lab's data - ATAC, Enhancer, Hi-C

```CMD
# ATAC-seq data
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/1b_stitzellab_dist_r2/EndoC_BH1_ATAC_broadPeak.tsv ^
  --out r2_data/1b_stitzellab_over_r2
```

> Input file/folder N     = [1] 1
>
> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 1
> File EndoC_BH1_ATAC_broadPeak... nrow= 5890.. done
>   Annotations occupied by SNPs  = [1] 100
>   SNPs in annotations           = [1] 151
>   Write file: r2_data/1b_stitzellab_over_r2/snp_EndoC_BH1_ATAC_broadPeak_151.bed
> Job done: 2020-04-21 10:55:14 for 0.9 sec

```CMD
# ChromHMM annotation
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/1b_stitzellab_dist_r2/EndoC_BH1_ChromHMM.tsv ^
  --out r2_data/1b_stitzellab_over_r2 ^
  --infotype tags
```

> Input file N    = [1] 1
> File EndoC_BH1_ChromHMM... nrow= 5890.. done
>   Annotations occupied by SNPs  = [1] 1099
>   SNPs in annotations           = [1] 5790
>
>  Annotations: 13 tags
>    1 Active.enhancer.1:        155..   Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Active.enhancer.1_155.bed
>    2 Active.enhancer.2:        42..    Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Active.enhancer.2_42.bed
>    3 Active.TSS:       73..    Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Active.TSS_73.bed
>    4 Bivalent.poised.TSS:      16..    Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Bivalent.poised.TSS_16.bed
>    5 Flanking.TSS:     21..    Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Flanking.TSS_21.bed
>    6 Genic.enhancer:   19..    Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Genic.enhancer_19.bed
>    7 Quiescent.low.signal:     1155..  Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Quiescent.low.signal_1155.bed
>    8 Repressed.polycomb:       74..    Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Repressed.polycomb_74.bed
>
>    9 Strong.transcription:     241..   Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Strong.transcription_241.bed
>    10 Weak.enhancer:   151..   Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Weak.enhancer_151.bed
>    11 Weak.repressed.polycomb: 2376..  Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Weak.repressed.polycomb_2376.bed
>    12 Weak.transcription:      1316..  Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Weak.transcription_1316.bed
>    13 Weak.TSS:        151..   Save at: r2_data/1b_stitzellab_over_r2/snp_tags_Weak.TSS_151.bed
> Job done: 2020-04-21 10:56:15 for 0.8 sec

```CMD
# Hi-C data for primary islet (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1b_stitzellab_dist_r2/Primary_Islet_HiC_gwas.tsv r2_data/1b_stitzellab_dist_r2/Primary_Islet_HiC_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1b_stitzellab_dist_r2/Primary_Islet_HiC_gwas.tsv, length= 6491,       overlap= 240
>   r2_data/1b_stitzellab_dist_r2/Primary_Islet_HiC_gene.tsv, length= 1217288,    overlap= 43158
>
>   Process gwas_loop.. 240.. done
>   Process gene_loop.. 43158.. done
>   Process merge for TSV.. 240.. [1] 1352    5
>   Write file: r2_data/summary_gene_r2/hic_Primary_Islet_HiC_gwas_179.tsv
>
>   Process merge for TSV.. 240.. [1] 103   4
>   Write file: r2_data/summary_r2/snp_hic_Primary_Islet_HiC_gwas_103.bed
> Job done: 2020-04-22 17:52:52 for 4.6 sec

```CMD
# Hi-C data for EndoC-BH1 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/1b_stitzellab_dist_r2/EndoC_BH1_HiC_gwas.tsv r2_data/1b_stitzellab_dist_r2/EndoC_BH1_HiC_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/1b_stitzellab_dist_r2/EndoC_BH1_HiC_gwas.tsv, length= 6638,   overlap= 426
>   r2_data/1b_stitzellab_dist_r2/EndoC_BH1_HiC_gene.tsv, length= 1372782,        overlap= 151476
>
>   Process gwas_loop.. 426.. done
>   Process gene_loop.. 151476.. done
>   Process merge for TSV.. 426.. [1] 5351    5
>   Write file: r2_data/summary_gene_r2/hic_EndoC_BH1_HiC_gwas_453.tsv
>
>   Process merge for TSV.. 426.. [1] 265   4
>   Write file: r2_data/summary_r2/snp_hic_EndoC_BH1_HiC_gwas_265.bed
> Job done: 2020-04-22 18:00:07 for 7.9 sec

### Pritchard lab's data - ATAC-seq

```CMD
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/2_pritchardlab_dist_r2 ^
  --out r2_data/2_pritchardlab_over_r2 ^
  --meta db_gwas/pritchardlab_meta.tsv
```

> ...
>
> File 1011.Naive_Teffs.S... nrow= 5891.. done
>   Annotations occupied by SNPs  = [1] 706
>   SNPs in annotations           = [1] 1153
>   Write file: r2_data/2_pritchardlab_over_r2/T-CD4-cells_stim/snp_1011.Naive_Teffs.S_1153.bed
> Job process: 0.3 sec
>
> Job done: 2020-04-21 11:23:42 for 1.2 min

### Tan lab's data - Enhancer

```CMD
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/3_tanlab_dist_r2 ^
  --out r2_data/3_tanlab_over_r2 ^
  --meta db_gwas/tanlab_meta.tsv
```

> ...
>
> File Treg_enh_t1d... nrow= 5891.. done
>   Annotations occupied by SNPs  = [1] 91
>   SNPs in annotations           = [1] 242
>   Write file: r2_data/3_tanlab_over_r2/Treg_enh/snp_Treg_enh_t1d_242.bed
>   6/6 done: r2_data/3_tanlab_dist_r2/Treg_enh_t1d.tsv
> Job process: 0.3 sec
>
> Job done: 2020-04-21 11:24:47 for 3 sec

### Allan lab's data - Hi-C

Memory error was occurred in 16 GB RAM system.

```CMD
# Hi-C data for B cells 1 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3a_allanlab_dist_r2/B1_gwas.tsv r2_data/3a_allanlab_dist_r2/B1_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3a_allanlab_dist_r2/B1_gwas.tsv, length= 13783,       overlap= 9628
>   r2_data/3a_allanlab_dist_r2/B1_gene.tsv, length= 3572972,     overlap= 2866785
>
>   Process gwas_loop.. 9628.. done
>   Process gene_loop.. 2866785.. done
>   Process merge for TSV.. 9628.. [1] 435077      5
>   Write file: r2_data/summary_gene_r2/hic_B1_gwas_4071.tsv
>
>   Process merge for TSV.. 9628.. [1] 3215    4
>   Write file: r2_data/summary_r2/snp_hic_B1_gwas_3214.bed
> Job done: 2020-04-22 18:09:10 for 7.6 min

```CMD
# Hi-C data for B cells 2 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3a_allanlab_dist_r2/B2_gwas.tsv r2_data/3a_allanlab_dist_r2/B2_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3a_allanlab_dist_r2/B2_gwas.tsv, length= 13680,       overlap= 9493
>   r2_data/3a_allanlab_dist_r2/B2_gene.tsv, length= 3597090,     overlap= 2878821
>
>   Process gwas_loop.. 9493.. done
>   Process gene_loop.. 2878821.. done
>   Process merge for TSV.. 9493.. [1] 444364      5
>   Write file: r2_data/summary_gene_r2/hic_B2_gwas_4342.tsv
>
>   Process extract for BED.. 9493.. [1] 3323    4
>   Write file: r2_data/summary_r2/snp_hic_B2_gwas_3322.bed
> Job done: 2020-04-22 18:18:17 for 7.5 min

```CMD
# Hi-C data for CD4 T cells 1 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3a_allanlab_dist_r2/CD4_T1_gwas.tsv r2_data/3a_allanlab_dist_r2/CD4_T1_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3a_allanlab_dist_r2/CD4_T1_gwas.tsv, length= 13279,   overlap= 8460
>   r2_data/3a_allanlab_dist_r2/CD4_T1_gene.tsv, length= 2968818, overlap= 2158049
>
>   Process gwas_loop.. 8460.. done
>   Process gene_loop.. 2158049.. done
>   Process merge for TSV.. 8460.. [1] 453403      5
>   Write file: r2_data/summary_gene_r2/hic_CD4_T1_gwas_3209.tsv
>
>   Process extract for BED.. 8460.. [1] 2950    4
>   Write file: r2_data/summary_r2/snp_hic_CD4_T1_gwas_2949.bed
> Job done: 2020-04-22 18:28:51 for 5.3 min

```CMD
# Hi-C data for CD4 T cells 2 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene hic_pair ^
  --base r2_data/3a_allanlab_dist_r2/CD4_T2_gwas.tsv r2_data/3a_allanlab_dist_r2/CD4_T2_gene.tsv ^
  --out r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3a_allanlab_dist_r2/CD4_T2_gwas.tsv, length= 16189,   overlap= 9289
>   r2_data/3a_allanlab_dist_r2/CD4_T2_gene.tsv, length= 3595761, overlap= 2848415
>
>   Process gwas_loop.. 9289.. done
>   Process gene_loop.. 2848415.. done
>   Process merge for TSV.. 9289.. [1] 534803      5
>   Write file: r2_data/summary_gene_r2/hic_CD4_T2_gwas_4339.tsv
>
>   Process extract for BED.. 9289.. [1] 3439    4
>   Write file: r2_data/summary_r2/snp_hic_CD4_T2_gwas_3438.bed
> Job done: 2020-04-22 18:37:05 for 7.3 min

```CMD
# Hi-C data for CD8 T cells 1 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3a_allanlab_dist_r2/CD8_T1_gwas.tsv r2_data/3a_allanlab_dist_r2/CD8_T1_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3a_allanlab_dist_r2/CD8_T1_gwas.tsv, length= 11041,   overlap= 7583
>   r2_data/3a_allanlab_dist_r2/CD8_T1_gene.tsv, length= 2846663, overlap= 2018663
>
>   Process gwas_loop.. 7583.. done
>   Process gene_loop.. 2018663.. done
>   Process merge for TSV.. 7583.. [1] 375074      5
>   Write file: r2_data/summary_gene_r2/hic_CD8_T1_gwas_3021.tsv
>
>   Process extract for BED.. 7583.. [1] 2906    4
>   Write file: r2_data/summary_r2/snp_hic_CD8_T1_gwas_2905.bed
> Job done: 2020-04-22 18:42:23 for 4.6 min

```CMD
# Hi-C data for CD8 T cells 2 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3a_allanlab_dist_r2/CD8_T2_gwas.tsv r2_data/3a_allanlab_dist_r2/CD8_T2_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3a_allanlab_dist_r2/CD8_T2_gwas.tsv, length= 11294,   overlap= 7864
>   r2_data/3a_allanlab_dist_r2/CD8_T2_gene.tsv, length= 2929292, overlap= 2125239
>
>   Process gwas_loop.. 7864.. done
>   Process gene_loop.. 2125239.. done
>   Process merge for TSV.. 7864.. [1] 338456      5
>   Write file: r2_data/summary_gene_r2/hic_CD8_T2_gwas_3386.tsv
>
>   Process extract for BED.. 7864.. [1] 3721    4
>   Write file: r2_data/summary_r2/snp_hic_CD8_T2_gwas_3720.bed
> Job done: 2020-04-22 18:47:33 for 4.9 min

### Aiden lab's data - Hi-C

```CMD
# Hi-C data for GM12878_primary (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3b_aidenlab_dist_r2/GM12878_primary_gwas.tsv r2_data/3b_aidenlab_dist_r2/GM12878_primary_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3b_aidenlab_dist_r2/GM12878_primary_gwas.tsv, length= 7506,   overlap= 614
>   r2_data/3b_aidenlab_dist_r2/GM12878_primary_gene.tsv, length= 1496697,        overlap= 162699
>
>   Process gwas_loop.. 614.. done
>   Process gene_loop.. 162699.. done
>   Process merge for TSV.. 614.. [1] 6923    5
>   Write file: r2_data/summary_gene_r2/hic_GM12878_primary_gwas_485.tsv
>
>   Process extract for BED.. 614.. [1] 347   4
>   Write file: r2_data/summary_r2/snp_hic_GM12878_primary_gwas_347.bed
> Job done: 2020-04-28 23:03:55 for 9.1 sec

```CMD
# Hi-C data for GM12878_replicate (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3b_aidenlab_dist_r2/GM12878_replicate_gwas.tsv r2_data/3b_aidenlab_dist_r2/GM12878_replicate_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3b_aidenlab_dist_r2/GM12878_replicate_gwas.tsv, length= 7312, overlap= 553
>   r2_data/3b_aidenlab_dist_r2/GM12878_replicate_gene.tsv, length= 1440751,      overlap= 159468
>
>   Process gwas_loop.. 553.. done
>   Process gene_loop.. 159468.. done
>   Process merge for TSV.. 553.. [1] 6055    5
>   Write file: r2_data/summary_gene_r2/hic_GM12878_replicate_gwas_562.tsv
>
>   Process extract for BED.. 553.. [1] 279   4
>   Write file: r2_data/summary_r2/snp_hic_GM12878_replicate_gwas_279.bed
> Job done: 2020-04-28 23:04:37 for 8.7 sec

```CMD
# Hi-C data for K562 (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/3b_aidenlab_dist_r2/K562_gwas.tsv r2_data/3b_aidenlab_dist_r2/K562_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/3b_aidenlab_dist_r2/K562_gwas.tsv, length= 6777,      overlap= 586
>   r2_data/3b_aidenlab_dist_r2/K562_gene.tsv, length= 1332121,   overlap= 148040
>
>   Process gwas_loop.. 586.. done
>   Process gene_loop.. 148040.. done
>   Process merge for TSV.. 586.. [1] 11288     5
>   Write file: r2_data/summary_gene_r2/hic_K562_gwas_433.tsv
>
>   Process extract for BED.. 586.. [1] 366   4
>   Write file: r2_data/summary_r2/snp_hic_K562_gwas_366.bed
> Job done: 2020-04-28 23:05:39 for 9.5 sec

### Ferrer lab's data - ATAC, Enhancer, pcHi-C

```CMD
# ATAC-seq data
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/4_ferrerlab_dist_r2/islet_atac_peaks.tsv ^
  --out r2_data/4_ferrerlab_over_r2
```

> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 1
> File islet_atac_peaks... nrow= 5896.. done
>   Annotations occupied by SNPs  = [1] 209
>   SNPs in annotations           = [1] 325
>   Write file: r2_data/4_ferrerlab_over_r2/snp_islet_atac_peaks_325.bed
> Job done: 2020-04-29 13:11:50 for 0.3 sec

```CMD
# ChromHMM annotation
Rscript postgwas-exe.r ^
  --dbfilt dist ^
  --base r2_data/4_ferrerlab_dist_r2/islet_active_enh.tsv ^
  --out r2_data/4_ferrerlab_over_r2
```

> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 1
> File islet_active_enh... nrow= 5891.. done
>   Annotations occupied by SNPs  = [1] 51
>   SNPs in annotations           = [1] 73
>   Write file: r2_data/4_ferrerlab_over_r2/snp_islet_active_enh_73.bed
> Job done: 2020-04-29 13:19:47 for 0.3 sec

```CMD
# pcHi-C data for human islet (TSV, BED)
Rscript postgwas-exe.r ^
  --dbgene	 hic_pair ^
  --base	 r2_data/4_ferrerlab_dist_r2/islet_pcHiC_gwas.tsv r2_data/4_ferrerlab_dist_r2/islet_pcHiC_gene.tsv ^
  --out		 r2_data/summary_gene_r2 r2_data/summary_r2 ^
  --bed		 TRUE
```

> ** Run function: db_gene.r/hic... ready
>   r2_data/4_ferrerlab_dist_r2/islet_pcHiC_gwas.tsv, length= 5890,       overlap= 0
>   r2_data/4_ferrerlab_dist_r2/islet_pcHiC_gene.tsv, length= 1163547,    overlap= 0
>
>   Process gwas_loop.. 0.. done
>   Process gene_loop.. 0.. done
>   Process merge for TSV.. 0.. [1] 0 0
>   Write file: r2_data/summary_gene_r2/hic_islet_pcHiC_gwas_0.tsv
>
>   Process extract for BED.. 0.. [1] 0 0
>   Write file: r2_data/summary_r2/snp_hic_islet_pcHiC_gwas_0.bed
> Job done: 2020-04-29 13:20:20 for 3.9 sec

# 6. BED union list

## General annotations: ENCODE Tfbs, UCSC gene regions

```CMD
Rscript postgwas-exe.r ^
  --dbfilt	dist ^
  --base	r2_data/distance_r2/encode_tfbs.tsv r2_data/distance_r2/cds_gene.tsv \^
  --out		r2_data/summary_r2
```

> ** Run function: db_filter.r/distance_filt_multi...
> Input file N    = [1] 2
> File encode_tfbs... nrow= 12201.. done
>   Annotations occupied by SNPs  = [1] 5230
>   SNPs in annotations           = [1] 1532
>   Write file: r2_data/summary/snp_encode_tfbs_1532.bed
>   1/2 done: r2_data/distance_r2/encode_tfbs.tsv
> Job process: 0.4 sec
>
> File cds_gene... nrow= 18361.. done
>   Annotations occupied by SNPs  = [1] 148
>   SNPs in annotations           = [1] 160
>   Write file: r2_data/summary/snp_cds_gene_160.bed
>   2/2 done: r2_data/distance_r2/cds_gene.tsv
> Job process: 0.2 sec
>
> Job done: 2020-04-10 14:34:09 for 0.6 sec

## Roadmap union list

```CMD
Rscript postgwas-exe.r ^
  --dbvenn		summ ^
  --base		r2_data/0_roadmap_over_r2 ^
  --out			r2_data/summary_r2 ^
  --sub_dir		FALSE ^
  --uni_save	TRUE
```

> ** Run function: db_venn.r/summ... ready
> 56 Files/folders input.
>   1 1 files in the ADRENAL
>   2 2 files in the Blood_B-cell_CD19p
>   3 4 files in the Blood_HSC_CD34p
>   4 2 files in the Blood_Leukemia_encode
>   5 1 files in the Blood_Lymphoblastoid_encode
>   6 1 files in the Blood_Monocytes_CD14
>   7 1 files in the Blood_Monocytes_CD14_encode
>   8 1 files in the Blood_Mononuclear_cell
>   9 1 files in the Blood_NK_cell_CD56
>   10 1 files in the Blood_Nutrophils_CD15
>   11 2 files in the Blood_T-cell_CD3
>   12 2 files in the Blood_T-cell_CD8
>   13 2 files in the Blood_Th_CD4
>   14 1 files in the Blood_Th_CD4p_CD25m
>   15 1 files in the Blood_Th_memory_CD4p_CD25m
>   16 1 files in the Blood_Th_naive_CD4p_CD25m
>   17 1 files in the Blood_Th_PMA-I_stim_CD4p_CD25m_IL17m
>   18 1 files in the Blood_Th17_PMA-I_stim_CD4p_CD25m_IL17p
>   19 1 files in the Blood_Tmem_CD4p_CD25int_CD127p
>   20 1 files in the Blood_Treg_CD4p_CD25p
>   21 1 files in the BONE_encode
>   22 12 files in the BRAIN
>   23 1 files in the BRAIN_encode
>   24 2 files in the BREAST
>   25 1 files in the BREAST_encode
>   26 1 files in the CERVIX_encode
>   27 8 files in the ESC
>   28 9 files in the ESC_DERIVED
>   29 3 files in the FAT
>   30 3 files in the GI_COLON
>   31 2 files in the GI_DUODENUM
>   32 1 files in the GI_ESOPHAGUS
>   33 3 files in the GI_INTESTINE
>   34 3 files in the GI_RECTUM
>   35 4 files in the GI_STOMACH
>   36 4 files in the HEART
>   37 5 files in the IPSC
>   38 1 files in the KIDNEY
>   39 1 files in the LIVER
>   40 1 files in the LIVER_encode
>   41 3 files in the LUNG
>   42 2 files in the LUNG_encode
>   43 5 files in the MUSCLE
>   44 2 files in the MUSCLE_encode
>   45 1 files in the MUSCLE_LEG
>   46 1 files in the OVARY
>   47 1 files in the PANCREAS
>   48 1 files in the PANCREAS_Islets
>   49 2 files in the PLACENTA
>   50 6 files in the SKIN
>   51 2 files in the SKIN_encode
>   52 1 files in the SPLEEN
>   53 2 files in the STROMAL_CONNECTIVE
>   54 2 files in the THYMUS
>   55 1 files in the VASCULAR
>   56 1 files in the VASCULAR_encode
> Total 127 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 127 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 2233  131
>   Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_2233.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-04-10 11:49:34 for 3.7 sec

```CMD
# Roadmap list by groups
Rscript postgwas-exe.r ^
  --dbvenn		summ ^
  --base		r2_data/0_roadmap_over_r2 ^
  --out			r2_data/summary_r2 ^
  --sub_dir		TRUE ^
  --uni_save	TRUE
```

> Total 56 sub-folder(s) is/are input
> Total 0 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
>   1 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_ADRENAL_285.bed
>   2 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_B-cell_CD19p_1030.bed
>   3 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_HSC_CD34p_835.bed
>   4 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Leukemia_encode_779.bed
>   5 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Lymphoblastoid_encode_860.bed
>   6 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Monocytes_CD14_615.bed
>   7 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Monocytes_CD14_encode_615.bed
>   8 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Mononuclear_cell_638.bed
>   9 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_NK_cell_CD56_697.bed
>   10 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Nutrophils_CD15_752.bed
>   11 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_T-cell_CD3_946.bed
>   12 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_T-cell_CD8_738.bed
>   13 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Th_CD4_701.bed
>   14 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Th_CD4p_CD25m_607.bed
>   15 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Th_memory_CD4p_CD25m_641.bed
>   16 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Th_naive_CD4p_CD25m_609.bed
>   17 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Th_PMA-I_stim_CD4p_CD25m_IL17m_613.bed
>   18 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Th17_PMA-I_stim_CD4p_CD25m_IL17p_615.bed
>
>   19 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Tmem_CD4p_CD25int_CD127p_619.bed
>   20 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_Blood_Treg_CD4p_CD25p_658.bed
>   21 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_BONE_encode_240.bed
>   22 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_BRAIN_442.bed
>   23 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_BRAIN_encode_223.bed
>   24 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_BREAST_383.bed
>   25 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_BREAST_encode_239.bed
>   26 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_CERVIX_encode_263.bed
>   27 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_ESC_335.bed
>   28 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_ESC_DERIVED_551.bed
>   29 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_FAT_483.bed
>   30 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_GI_COLON_559.bed
>   31 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_GI_DUODENUM_471.bed
>   32 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_GI_ESOPHAGUS_266.bed
>   33 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_GI_INTESTINE_481.bed
>   34 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_GI_RECTUM_517.bed
>   35 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_GI_STOMACH_577.bed
>   36 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_HEART_365.bed
>   37 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_IPSC_306.bed
>   38 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_KIDNEY_241.bed
>   39 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_LIVER_356.bed
>   40 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_LIVER_encode_276.bed
>   41 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_LUNG_481.bed
>   42 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_LUNG_encode_332.bed
>   43 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_MUSCLE_521.bed
>   44 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_MUSCLE_encode_275.bed
>   45 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_MUSCLE_LEG_204.bed
>   46 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_OVARY_184.bed
>   47 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_PANCREAS_236.bed
>   48 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_PANCREAS_Islets_185.bed
>   49 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_PLACENTA_430.bed
>   50 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_SKIN_424.bed
>   51 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_SKIN_encode_335.bed
>   52 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_SPLEEN_438.bed
>   53 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_STROMAL_CONNECTIVE_256.bed
>   54 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_THYMUS_795.bed
>   55 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_VASCULAR_161.bed
>   56 Write a BED file: r2_data/summary/snp_union_0_roadmap_over_r2_VASCULAR_encode_180.bed
> Job done: 2020-04-10 11:51:51 for 4.2 sec

## GTEx eQTLs

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/gtex_eqtl ^
  --out r2_data/summary_r2 ^
  --sub_dir FALSE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
> 49 Files/folders input.
>   1 r2_data/gtex_eqtl/snp_gtex_Adipose_Subcutaneous_2444.bed
>   2 r2_data/gtex_eqtl/snp_gtex_Adipose_Visceral_Omentum_1982.bed
>   3 r2_data/gtex_eqtl/snp_gtex_Adrenal_Gland_1554.bed
>   4 r2_data/gtex_eqtl/snp_gtex_Artery_Aorta_1906.bed
>   5 r2_data/gtex_eqtl/snp_gtex_Artery_Coronary_1197.bed
>   6 r2_data/gtex_eqtl/snp_gtex_Artery_Tibial_2086.bed
>   7 r2_data/gtex_eqtl/snp_gtex_Brain_Amygdala_1021.bed
>   8 r2_data/gtex_eqtl/snp_gtex_Brain_Anterior_cingulate_cortex_BA24_1027.bed
>   9 r2_data/gtex_eqtl/snp_gtex_Brain_Caudate_basal_ganglia_1265.bed
>   10 r2_data/gtex_eqtl/snp_gtex_Brain_Cerebellar_Hemisphere_1621.bed
>   11 r2_data/gtex_eqtl/snp_gtex_Brain_Cerebellum_1985.bed
>   12 r2_data/gtex_eqtl/snp_gtex_Brain_Cortex_1338.bed
>   13 r2_data/gtex_eqtl/snp_gtex_Brain_Frontal_Cortex_BA9_1028.bed
>   14 r2_data/gtex_eqtl/snp_gtex_Brain_Hippocampus_1173.bed
>   15 r2_data/gtex_eqtl/snp_gtex_Brain_Hypothalamus_1184.bed
>   16 r2_data/gtex_eqtl/snp_gtex_Brain_Nucleus_accumbens_basal_ganglia_1196.bed
>   17 r2_data/gtex_eqtl/snp_gtex_Brain_Putamen_basal_ganglia_1071.bed
>   18 r2_data/gtex_eqtl/snp_gtex_Brain_Spinal_cord_cervical_c-1_678.bed
>   19 r2_data/gtex_eqtl/snp_gtex_Brain_Substantia_nigra_528.bed
>   20 r2_data/gtex_eqtl/snp_gtex_Breast_Mammary_Tissue_1876.bed
>   21 r2_data/gtex_eqtl/snp_gtex_Cells_Cultured_fibroblasts_1866.bed
>   22 r2_data/gtex_eqtl/snp_gtex_Cells_EBV-transformed_lymphocytes_1282.bed
>   23 r2_data/gtex_eqtl/snp_gtex_Colon_Sigmoid_1953.bed
>   24 r2_data/gtex_eqtl/snp_gtex_Colon_Transverse_2038.bed
>   25 r2_data/gtex_eqtl/snp_gtex_Esophagus_Gastroesophageal_Junction_1982.bed
>   26 r2_data/gtex_eqtl/snp_gtex_Esophagus_Mucosa_2367.bed
>   27 r2_data/gtex_eqtl/snp_gtex_Esophagus_Muscularis_2149.bed
>   28 r2_data/gtex_eqtl/snp_gtex_Heart_Atrial_Appendage_2096.bed
>   29 r2_data/gtex_eqtl/snp_gtex_Heart_Left_Ventricle_1743.bed
>   30 r2_data/gtex_eqtl/snp_gtex_Kidney_Cortex_323.bed
>   31 r2_data/gtex_eqtl/snp_gtex_Liver_1457.bed
>   32 r2_data/gtex_eqtl/snp_gtex_Lung_2095.bed
>   33 r2_data/gtex_eqtl/snp_gtex_Minor_Salivary_Gland_868.bed
>   34 r2_data/gtex_eqtl/snp_gtex_Muscle_Skeletal_2073.bed
>   35 r2_data/gtex_eqtl/snp_gtex_Nerve_Tibial_2720.bed
>   36 r2_data/gtex_eqtl/snp_gtex_Ovary_1068.bed
>   37 r2_data/gtex_eqtl/snp_gtex_Pancreas_1927.bed
>   38 r2_data/gtex_eqtl/snp_gtex_Pituitary_1815.bed
>   39 r2_data/gtex_eqtl/snp_gtex_Prostate_1338.bed
>   40 r2_data/gtex_eqtl/snp_gtex_Skin_Not_Sun_Exposed_Suprapubic_2277.bed
>   41 r2_data/gtex_eqtl/snp_gtex_Skin_Sun_Exposed_Lower_leg_2463.bed
>   42 r2_data/gtex_eqtl/snp_gtex_Small_Intestine_Terminal_Ileum_1403.bed
>   43 r2_data/gtex_eqtl/snp_gtex_Spleen_1849.bed
>   44 r2_data/gtex_eqtl/snp_gtex_Stomach_1635.bed
>   45 r2_data/gtex_eqtl/snp_gtex_Testis_2366.bed
>   46 r2_data/gtex_eqtl/snp_gtex_Thyroid_2527.bed
>   47 r2_data/gtex_eqtl/snp_gtex_Uterus_776.bed
>   48 r2_data/gtex_eqtl/snp_gtex_Vagina_925.bed
>   49 r2_data/gtex_eqtl/snp_gtex_Whole_Blood_2296.bed
> Total 49 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 49 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 3921   53
>   Write a BED file: r2_data/summary_r2/snp_union_gtex_eqtl_3921.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-04-21 12:05:11 for 3.9 sec

## Additional data

### Melton lab union list

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/1_meltonlab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir FALSE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
> 4 Files/folders input.
>   1 6 files in the beta_ATAC
>   2 6 files in the beta_enh
>   3 4 files in the other_ATAC
>   4 9 files in the other_enh
> Total 25 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 25 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 739  29
>   Write a BED file: r2_data/summary_r2/snp_union_1_meltonlab_over_r2_739.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-04-21 12:07:29 for 12.7 sec

```CMD
# union list by groups
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/1_meltonlab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir TRUE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
> 4 Files/folders input.
>   1 sub_dir 1: 6 file(s) in the beta_ATAC folder
>   2 sub_dir 2: 6 file(s) in the beta_enh folder
>   3 sub_dir 3: 4 file(s) in the other_ATAC folder
>   4 sub_dir 4: 9 file(s) in the other_enh folder
> Total 4 sub-folder(s) is/are input
> Total 0 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
>   1 Write a BED file: r2_data/summary_r2/snp_union_1_meltonlab_over_r2_beta_ATAC_234.bed
>   2 Write a BED file: r2_data/summary_r2/snp_union_1_meltonlab_over_r2_beta_enh_279.bed
>   3 Write a BED file: r2_data/summary_r2/snp_union_1_meltonlab_over_r2_other_ATAC_195.bed
>   4 Write a BED file: r2_data/summary_r2/snp_union_1_meltonlab_over_r2_other_enh_512.bed
> Job done: 2020-04-21 12:08:27 for 3.4 sec

### Suh lab union list

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/1a_suhlab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir FALSE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
> 3 Files/folders input.
>   1 r2_data/1a_suhlab_over_r2/snp_Scbeta_YZ_pr1_pooled.pf_peaks_168.bed
>   2 r2_data/1a_suhlab_over_r2/snp_Scbeta_YZ_pr2_pooled.pf_peaks_177.bed
>   3 r2_data/1a_suhlab_over_r2/snp_SCbeta_YZ_pseudo_pool_ATAC_140.bed
> Total 3 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Scbeta_YZ_pr1_pooled.pf_peaks_168
>   Read 2: snp_Scbeta_YZ_pr2_pooled.pf_peaks_177
>   Read 3: snp_SCbeta_YZ_pseudo_pool_ATAC_140
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 228   7
>   Write a BED file: r2_data/summary_r2/snp_union_1a_suhlab_over_r2_228.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-04-21 12:09:30 for 6.4 sec

### Stitzel lab union list

```CMD
# union list by groups
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/1b_stitzellab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir TRUE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
> 11 Files/folders input.
>   1 sub_dir 1: 1 file(s) in the ATAC folder
>   2 sub_dir 2: 4 file(s) in the ChromHMM_enh folder
>   3 r2_data/1b_stitzellab_over_r2/snp_tags_Active.TSS_73.bed
>   4 r2_data/1b_stitzellab_over_r2/snp_tags_Bivalent.poised.TSS_16.bed
>   5 r2_data/1b_stitzellab_over_r2/snp_tags_Flanking.TSS_21.bed
>   6 r2_data/1b_stitzellab_over_r2/snp_tags_Quiescent.low.signal_1155.bed
>   7 r2_data/1b_stitzellab_over_r2/snp_tags_Repressed.polycomb_74.bed
>   8 r2_data/1b_stitzellab_over_r2/snp_tags_Strong.transcription_241.bed
>   9 r2_data/1b_stitzellab_over_r2/snp_tags_Weak.repressed.polycomb_2376.bed
>   10 r2_data/1b_stitzellab_over_r2/snp_tags_Weak.transcription_1316.bed
>   11 r2_data/1b_stitzellab_over_r2/snp_tags_Weak.TSS_151.bed
> Total 2 sub-folder(s) is/are input
> Total 9 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
>   1 Write a BED file: r2_data/summary_r2/snp_union_1b_stitzellab_over_r2_ATAC_151.bed
>   2 Write a BED file: r2_data/summary_r2/snp_union_1b_stitzellab_over_r2_ChromHMM_enh_367.bed
> Job done: 2020-04-21 12:47:20 for 3 sec

### Pritchard lab union list

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/2_pritchardlab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir FALSE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
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
>   Returned union list dim       = [1] 1312  179
>   Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_1312.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-04-21 12:14:41 for 36.9 sec

```CMD
# union list by groups
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/2_pritchardlab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir TRUE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
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
>   1 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_B-cells_rest_1309.bed
>   2 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_B-cells_stim_1312.bed
>   3 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_Gamma_delta_T_rest_1283.bed
>   4 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_Gamma_delta_T_stim_1309.bed
>   5 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_Monocytes_rest_1309.bed
>   6 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_Monocytes_stim_1312.bed
>   7 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_Myeloid_DCs_rest_1196.bed
>   8 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_NK-cells_rest_1312.bed
>   9 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_NK-cells_stim_1312.bed
>   10 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_pDCs_rest_1198.bed
>   11 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_Plasmablasts_rest_1302.bed
>   12 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_T-CD4-cells_rest_1312.bed
>   13 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_T-CD4-cells_stim_1312.bed
>   14 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_T-CD8-cells_rest_1312.bed
>   15 Write a BED file: r2_data/summary_r2/snp_union_2_pritchardlab_over_r2_T-CD8-cells_stim_1312.bed
> Job done: 2020-04-21 12:15:28 for 5.4 sec

### Tan lab union list

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/3_tanlab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir FALSE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
> 2 Files/folders input.
>   1 3 files in the Th1_enh
>   2 3 files in the Treg_enh
> Total 6 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Th1_enh_251
>   Read 2: snp_Th1_enh_control_222
>   Read 3: snp_Th1_enh_t1d_231
>   Read 4: snp_Treg_enh_285
>   Read 5: snp_Treg_enh_control_278
>   Read 6: snp_Treg_enh_t1d_242
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 321  10
>   Write a BED file: r2_data/summary_r2/snp_union_3_tanlab_over_r2_321.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-04-21 12:16:38 for 7.5 sec

```CMD
# union list by groups
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2_data/3_tanlab_over_r2 ^
  --out r2_data/summary_r2 ^
  --sub_dir TRUE ^
  --uni_save TRUE
```

> ** Run function: db_venn.r/summ... ready
> 2 Files/folders input.
>   1 sub_dir 1: 3 file(s) in the Th1_enh folder
>   2 sub_dir 2: 3 file(s) in the Treg_enh folder
> Total 2 sub-folder(s) is/are input
> Total 0 file(s) is/are input.
>
> Option sub_dir = TRUE, summary table are not going to be generated.
>   1 Write a BED file: r2_data/summary_r2/snp_union_3_tanlab_over_r2_Th1_enh_251.bed
>   2 Write a BED file: r2_data/summary_r2/snp_union_3_tanlab_over_r2_Treg_enh_285.bed
> Job done: 2020-04-21 12:17:13 for 3.2 sec

### Ferrer lab union list

Copy the files in the `r2_data/4_ferrerlab_over_r2` to `r2_data/summary_r2`.

* `snp_islet_active_enh_73.bed` -> `snp_union_4_ferrerlab_over_r2_islet_active_enh_73.bed`
* `snp_islet_atac_peaks_325` -> `snp_union_4_ferrerlab_over_r2_islet_atac_peaks_325.bed`

# 7. Summary annotations

## SNP level summary

```CMD
Rscript postgwas-exe.r ^
  --dbvenn		 summ ^
  --base		 r2_data/summary_r2 r2_data/gtex_eqtl ^
  --out			 r2_data ^
  --sub_dir		 FALSE ^
  --uni_save	 FALSE ^
  --ann_gwas	 r2_data/gwas_biomart_fill.tsv ^
  --ann_encd	 r2_data/distance_r2/encode_tfbs.tsv ^
  --ann_near	 r2_data/distance_r2/nearest_gene.tsv ^
  --ann_cds		 r2_data/distance_r2/cds_gene.tsv
```

> ** Run function: db_venn.r/summ... ready
> 164 Files/folders input.
>   1 r2_data/gtex_eqtl/snp_gtex_Adipose_Subcutaneous_2444.bed
>   2 r2_data/gtex_eqtl/snp_gtex_Adipose_Visceral_Omentum_1982.bed
>
>   ...
>
>   163 r2_data/summary_r2/snp_union_4_ferrerlab_over_r2_islet_atac_peaks_325.bed
>   164 r2_data/summary_r2/snp_union_gtex_eqtl_3921.bed
> Total 164 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
> Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
>   no lines available in input
> In addition: Warning message:
> package 'eulerr' was built under R version 3.6.2
>   [ERROR] r2_data/summary_r2/snp_hic_islet_pcHiC_gwas_0.bed
>   Read 164 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 10472   167
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 5891   11
>   Merge dim     = [1] 5890  171
>   Write a CSV file: r2_data/summary_r2_gwas.csv
>
>   ENCODE dim    = [1] 12201    13
>   Merge dim     = [1] 7828  166
>   Write a CSV file: r2_data/summary_r2_encode.csv
>
>   Nearest gene dim      = [1] 6511    9
>   Search biomaRt... 424.. 407.. [1] 6511    5
>   CDS dim               = [1] 18361    11
>   Merge dim             = [1] 6890  171
>   Write a CSV file: r2_data/summary_r2_nearest.csv
>
> Job done: 2020-04-29 15:04:57 for 18.2 sec

## Gene lavel summary

### GTEx data

```CMD
Rscript postgwas-exe.r ^
  --dbgene gtex_pair ^
  --base r2_data/gtex_signif_3921.tsv ^
  --out r2_data/summary_gene_r2
```

>** Run function: db_gene.r/gtex_pair... ready
>  GTEx dim      = [1] 270976      9
>  Extracting Ensgid.. done
>  Spreading GTEx data.. [1] 30218    51
>  Write file: r2_data/summary_gene_r2/gtex_eqtl_286.tsv
>Job done: 2020-04-22 20:57:38 for 3.9 sec

### Summary Hi-C, GTEx, Nearest genes

```CMD
Rscript postgwas-exe.r
  --dbgene summary
  --base r2_data/summary_gene_r2
  --nearest r2_data/summary_r2_nearest.csv
  --out r2_data
```

> ** Run function: db_gene.r/summary_gene... ready
> 21 Files/folders input.
>   1 r2_data/summary_gene_r2/gtex_eqtl_286.tsv
>   2 r2_data/summary_gene_r2/hic_B1_gwas_4071.tsv
>   3 r2_data/summary_gene_r2/hic_B2_gwas_4342.tsv
>   4 r2_data/summary_gene_r2/hic_CD4_T1_gwas_3209.tsv
>   5 r2_data/summary_gene_r2/hic_CD4_T2_gwas_4339.tsv
>   6 r2_data/summary_gene_r2/hic_CD8_T1_gwas_3021.tsv
>   7 r2_data/summary_gene_r2/hic_CD8_T2_gwas_3386.tsv
>   8 r2_data/summary_gene_r2/hic_EndoC_BH1_HiC_gwas_453.tsv
>   9 r2_data/summary_gene_r2/hic_GM12878_primary_gwas_485.tsv
>   10 r2_data/summary_gene_r2/hic_GM12878_replicate_gwas_562.tsv
>   11 r2_data/summary_gene_r2/hic_islet_pcHiC_gwas_2131.tsv
>   12 r2_data/summary_gene_r2/hic_K562_gwas_433.tsv
>   13 r2_data/summary_gene_r2/hic_mono_1018_gwas_1440.tsv
>   14 r2_data/summary_gene_r2/hic_mono_1027_gwas_1345.tsv
>   15 r2_data/summary_gene_r2/hic_mono_1035_gwas_1699.tsv
>   16 r2_data/summary_gene_r2/hic_mono_1046_gwas_1871.tsv
>   17 r2_data/summary_gene_r2/hic_Primary_Islet_HiC_gwas_179.tsv
>   18 r2_data/summary_gene_r2/hic_T_1018_gwas_1233.tsv
>   19 r2_data/summary_gene_r2/hic_T_1027_gwas_1462.tsv
>   20 r2_data/summary_gene_r2/hic_T_1035_gwas_1626.tsv
>   21 r2_data/summary_gene_r2/hic_T_1046_gwas_1522.tsv
> Total 21 files are input.
> Nearest file: r2_data/summary_r2_nearest.csv
>   Read gene-snp pair files... done
>   Read nearest gene file... [1] 6890  172
> [1] 3416    3
>
>   Parsing gene name... 3415.. [1] 3416    3
>   ENSGid-Rsid list... 173286.. [1] 173286      3
>   Merging biomaRt annotations.. [1] 173292      5
>   Merging eQTL, Hi-C, nearest genes... [1] 173292     27
>   Merging GTEx eQTL SNP slopes... [1] 173292     78
>
>   Write file: r2_data/summary_gene_r2_pairs.csv
> Job done: 2020-04-29 15:08:22 for 1.5 min

### Pivoting Hi-C, GTEx, Nearest genes

```CMD
# Nearest genes
Rscript postgwas-exe.r
  --dbgene	 pivot_gene
  --base	 r2_data/summary_gene_r2_nearest_pivot.tsv
  --out		 r2_data
```

> ** Run function: db_gene.r/pivot_gene... ready
>   Read summary_gene_pair TSV file... [1] 2379   49
>   Extract ensgids... union table.. [1] 424  49
>   Write file: r2_data/summary_gene_pivot.csv
> Job done: 2020-04-29 16:12:46 for 0.1 sec

File name change: `summary_gene_pivot.csv` -> `summary_gene_nearest_pivot.csv`

```CMD
# Nearest, GTEx, Hi-C
Rscript postgwas-exe.r
  --dbgene	 pivot_gene
  --base	 r2_data/summary_gene_r2_total_pivot.tsv
  --out		 r2_data
```

> ** Run function: db_gene.r/pivot_gene... ready
>     Read summary_gene_pair TSV file... [1] 57375    54
>     Extract ensgids... union table.. [1] 8294   54
>     Write file: r2_data/summary_gene_pivot.csv
> Job done: 2020-04-29 18:11:44 for 1 sec

File name change: `summary_gene_pivot.tsv` -> `summary_gene_total_pivot.tsv`

## DAVID GO/KEGG analysis

```CMD
Rscript postgwas-exe.r ^
  --dbgene		 david_go ^
  --base		 r2_data/go_analysis_r2.tsv ^
  --out			 r2_data ^
  --criteria	 0.001 ^
  --stat		 pval ^
  --dataset		 Tissue_6
```

> ** Run function: db_gene.r/david_go... ready
>   Read DAVID result file        [1] 788  14
>   Filtering dataset = Tissue_6... [1] 143  14
>
>   Pval criteria = 0.001
>   Filtered GO terms     = [1] 23 14
>
>   Extract the gene list... 34.. done
>   Venn analysis... [1] 34 24
>   Write file: r2_data/go_analysis.Tissue_6.pval-0.001.csv
> Job done: 2020-04-12 16:52:00 for 0.1 sec



# #. Comparing R2 and Dprime effects

## Total LDlink data

```CMD
Rscript postgwas-exe.r
  --ldlink filter
  --base db_gwas/gwas_5e-08_152.tsv db_gwas/ldlink
  --out r2_data
  --r2d 6
```

> ** Run function ldlink_filter...
> Read download files... 152
> Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
>   line 1 did not have 11 elements
>   [ERROR] rs75793288
>   Read LDlink results           = [1] 303054     12
> No filtering:
>   Filtered data dimension       = [1] 303040      3
>   Excluded no rsid elements     = [1] 16522
>
>   Merged table          = [1] 286518      3
>   Write file: r2_data/gwas_ldlink_r0.2.csv
> Basic summary of LDlink results:
>   SNP Tier 1                    = 152
>   SNP Tier 2                    = 135669
>   SNP candidates                = 135821
>   SNP source annotation table   = [1] 135821      2
> Add annotations:
>   LD block annotation... [1] 88
> Search biomart for SNP coordinates:
>   Query SNPs            = [1] 135821
> [1] 134092      4
> Batch submitting query [=====>-------------------------]  18% eta:  2hError in getBM(attributes = snp_attr1, filters = "snp_filter", values = snp_cand,  :
>   The query to the BioMart webservice returned an invalid result: biomaRt expected a character string of length 1.
> Please report this on the support site at http://support.bioconductor.org
> Calls: gwas_ldlink -> ldlink_filter -> %>% -> eval -> eval -> getBM
> Execution halted

# #. To identify ethnic effect

## Downloading LDlink data

```CMD
# Download LD SNPs from 152 query SNPs takes about ~1.2 hrs
# --popul CEU TSI FIN GBR IBS 
Rscript postgwas-exe.r
  --ldlink down
  --base db_gwas/gwas_5e-08_152.tsv
  --out db_gwas/ldlink2
  --popul ALL
```

> ** Run function ldlink_down... 152..
>   error: rs75793288 is not in 1000G reference panel.,
> done
>   Files are moved to target folder:     db_gwas/ldlink2
> Job done: 2020-04-07 01:50:28 for 3.5 hr

## Filtering the LDlink data

```CMD
Rscript postgwas-exe.r
  --ldlink filter
  --base db_gwas/gwas_5e-08_152.tsv db_gwas/ldlink2
  --out r2_data
  --r2d 6
```

> ** Run function ldlink_filter...
> Read download files... 152
> Error in file(file, "rt") : cannot open the connection
> In addition: Warning message:
> In file(file, "rt") :
>   cannot open file 'db_gwas/ldlink2/rs75793288.txt': No such file or directory
>   [ERROR] rs75793288
>   Read LDlink results           = [1] 277840     12
> Minimum filtering by r2 > 0.2:
>   Filtered data dimension       = [1] 23484     3
>   Excluded no rsid elements     = [1] 603
>
>   Merged table          = [1] 22881     3
>   Write file: r2_data/gwas_ldlink_r0.2.csv
> Basic summary of LDlink results:
>   SNP Tier 1                    = 152
>   SNP Tier 2                    = 16698
>   SNP candidates                = 16850
>   SNP source annotation table   = [1] 16850     2
> Add annotations:
>   LD block annotation... [1] 107
> Search biomart for SNP coordinates:
>   Query SNPs            = [1] 16850
> [1] 16620     4
> Batch submitting query [=========================>-----]  82% eta:  2m
> Error in curl::curl_fetch_memory(url, handle = handle) :
>   Operation was aborted by an application callback
> Calls: gwas_ldlink ... request_fetch -> request_fetch.write_memory -> <Anonymous>
> Execution halted