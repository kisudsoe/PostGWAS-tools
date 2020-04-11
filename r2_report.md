# Type 1 Diabetes (r2 >0.6)

This is the result from the LDlink threshold as r<sup>2</sup> >0.6.

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

# 3. Filtering the annotations

UCSC annotations, Roadmap, GTEx, Stitzel lab's Hi-C data were already filtered in `r2d1_report.md`.

# 4. Distances from the annotations

## For general annotations: CDS_gene, encode_tfbs, nearest_gene

```bash
bash distance_r2.sh
```

## For roadmap annotations

```bash
bash 0_roadmap_dist_r2.sh
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

# 6. BED union list

## General annotations: ENCODE Tfbs, UCSC gene regions

```CMD
Rscript postgwas-exe.r \
  --dbfilt	dist \
  --base	r2_data/distance_r2/encode_tfbs.tsv r2_data/distance_r2/cds_gene.tsv \
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

### Roadmap list by groups

```CMD
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

# 7. Summary annotations

## SNP level summary

```CMD
Rscript postgwas-exe.r ^
  --dbvenn		summ ^
  --base		r2_data/summary_r2 ^
  --out			r2_data ^
  --sub_dir		FALSE ^
  --uni_save	FALSE ^
  --ann_gwas	r2_data/gwas_biomart_fill.tsv ^
  --ann_encd	r2_data/distance_r2/encode_tfbs.tsv ^
  --ann_near	r2_data/distance_r2/nearest_gene.tsv ^
  --ann_cds		r2_data/distance_r2/cds_gene.tsv
```

> ** Run function: db_venn.r/summ... ready
> 64 Files/folders input.
>   1 r2_data/summary/snp_cds_gene_160.bed
>
>   ...
>
>   64 r2_data/summary/snp_union_0_roadmap_over_r2_VASCULAR_encode_180.bed
> Total 64 file(s) is/are input.
>
> Total 64 file(s) is/are input.
>
> ** Run function: db_venn.r/venn_bed...
>   Read 64 files
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 4671   68
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 5891   11
>   Merge dim     = [1] 5890   72
>   Write a CSV file: r2_data/summary_gwas.csv
>
>   ENCODE dim    = [1] 12201    13
>   Merge dim     = [1] 7828   67
>   Write a CSV file: r2_data/summary_encode.csv
>
>   Nearest gene dim      = [1] 6511    9
>   Search biomaRt... 424.. 407.. [1] 6511    5
>   CDS dim               = [1] 18361    11
>   Merge dim             = [1] 6890   72
>   Write a CSV file: r2_data/summary_nearest.csv
>
> Job done: 2020-04-10 14:03:52 for 31.2 sec

## DAVID GO/KEGG analysis

```CMD
Rscript postgwas-exe.r ^
  --dbgene	 david_go ^
  --base	 r2_data/go_analysis_r2.tsv ^
  --out		 r2_data ^
  --fdr		 0.05
```

> ** Run function: db_gene.r/david_go... ready
>   Read DAVID result file        [1] 686  14
>   FDR criteria  = 0.05
>   Filtered GO terms     = [1] 54 14
>   Extract the gene list... 17.. done
>   Venn analysis... [1] 17 55
>   Write file: r2_data/go_analysis.csv
> Job done: 2020-04-11 17:32:23 for 0.1 sec



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