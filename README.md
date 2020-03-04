# Post-GWAS analysis : Type 1 Diabetes

This is a log file for analyzing "Type I Diabetes Mellitus" (EFO0001359) from GWAS Catalog.

As of 2020-01-20, the GWAS Catalog contains 4,390 publications and 171,674 associations.
GWAS Catalog data is currently mapped to Genome Assembly GRCh38.p13 and dbSNP Build 153.

And LDlink results are returned with **hg19** coordinates.

## Brief result description

![](C:\Users\sk4748\OneDrive\Suh's Lab\2020-03 T1D v2.0\Fig1a.png)

![](C:\Users\sk4748\OneDrive\Suh's Lab\2020-03 T1D v2.0\Fig1b.png)



![](C:\Users\sk4748\OneDrive\Suh's Lab\2020-03 T1D v2.0\Fig2.png)



![](C:\Users\sk4748\OneDrive\Suh's Lab\2020-03 T1D v2.0\Fig3.png)



**Additional tools to be added:**

* De novo motif: TrawlerWeb [http://trawler.erc.monash.edu.au](http://trawler.erc.monash.edu.au/)
  * Dang et al., 2018, BMC Genomics, [10.1186/s12864-018-4630-0](https://doi.org/10.1186/s12864-018-4630-0)
  * BED format motif list entry (hg19/hg38)
* Core transcriptional Regulatory Circuitries: dbCoRC http://dbcorc.cam-su.org/
  * Huang et al., 2018, Nucleic Acids Res, [10.1093/nar/gkx796](https://dx.doi.org/10.1093%2Fnar%2Fgkx796), pmid [28977473](https://www.ncbi.nlm.nih.gov/pubmed/28977473)
  * Gene symbol list entry (gh19)

# 1. From the seed SNPs (EFO0001359) to candidate risk SNPs

## Generating pivot tables

From the downloaded GWAS Catalog data, basic statistics of the trait are need to know. This function will generate three summary pivot tables.

Usage: `Rscript postgwas-exe.r --gwas <options> --base <input file> --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
	--gwas trait gene study ^
    --base db_gwas/EFO0001359_2020-02-19.tsv ^
    --out db_gwas
```

> ** Run function trait_pivot:
> Write TRAITS pivot:     db_gwas/EFO0001359_2020-02-19_snps.tsv
>
> ** Run function gene_pivot:
> Write gene pivot:       db_gwas/EFO0001359_2020-02-19_genes.tsv
>
> ** Run function study_pivot:
> Write study pivot:      db_gwas/EFO0001359_2020-02-19_studies.tsv

## Filtering the GWAS SNPs

Using the downloaded GWAS Catalog file, SNPs are needed to filter by the P-value criteria.

Usage: `Rscript postgwas-exe.r --gwas <option> --base <input file> --out <out folder> --p.criteria <number>`

```CMD
Rscript postgwas-exe.r ^
	--gwas filter ^
    --base db_gwas/EFO0001359_2020-02-19.tsv ^
    --out db_gwas ^
    --p.criteria 5e-8
```

> ** Run function gwas_filt:
>     gwas dim              = [1] 330  38
>     gwas (5e-08)          = [1] 192  38
> Write gwas filter:      db_gwas/gwas_5e-08_152.tsv



## Downloading LDlink data

To expand the seed SNPs to the neighbor SNPs within the LD block, LS associated SNPs from the seed SNPs are downloaded from the LDlink DB (Machiela et al, 2015, Bioinformatics, pmid 26139635).

Usage: `Rscript postgwas-exe.r --ldlink down --base <input file> --out <out folder> --popul <options>`

```CMD
# Download LD SNPs from 152 query SNPs takes about ~1.2 hrs
Rscript postgwas-exe.r ^
    --ldlink down ^
    --base db_db/gwas_5e-08_152.tsv ^
    --out db_gwas/ldlink ^
    --popul CEU TSI FIN GBR IBS
```

> ** Run function ldlink_dn... 152.. done
>   Files are moved to target folder: db_gwas/ldlink
> Job done: 2020-02-21 15:05:45 for 1.2 hr
> Warning message:
> package 'LDlinkR' was built under R version 3.6.2

## Filtering the LDlink data

Filtering the LDlink data by the criteria, r<sup>2</sup> > 0.6 and D' = 1.

Usage: `Rscript postgwas-exe.r --ldlink filter --base <input file> <input folder> --out <out folder> --r2d <option>`

```CMD
Rscript postgwas-exe.r ^
    --ldlink filter ^
    --base db_gwas/gwas_5e-08_152.tsv db_gwas/ldlink ^
    --out r2d1_data_gwas ^
    --r2d 1
```

> ** Run function ldlink_ln...
> Read download files... 152
> Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
>   line 1 did not have 11 elements
>   rs75793288
>   Read LDlink results           = [1] 303054     12
> Filtering by "r2 > 0.6 and Dprime = 1":
>   Filtered data dimension       = [1] 2147    3
>   Excluded no rsid elements     = [1] 25
> Basic summary of LDlink results:
>   SNP Tier 1                    = 152
>   SNP Tier 2                    = 1851
>   SNP candidates        = 2003
> Search biomart for SNP coordinates:
>   Query SNPs            = [1] 2003
>   Hg19 result table     = [1] 1992    4
>   Hg38 result table     = [1] 2003    4
>   Merged table          = [1] 2126    9
> Write file:     r2d1_data_gwas/gwas_biomart.tsv
> Job done: 2020-02-25 19:14:14 for 11.7 sec

## Generating BED files for coordinates from hg19 and hg38

To prepare candidate SNPs, you have to fill the "NA" coordinate information from the **hg19** and **hg38** columns of the "--ldlink fl" function result file (e.g. `gwas_biomart.tsv`) by searching from the Ensembl biomart. As a note, I found hg19 SNP `rs71080526` have two hg38 SNPs, `rs67544786` and `rs1553647310`. Among those two SNPs, I selected `rs67544786` as same for the Ensemble homepage searching result. And I found the 11 SNPs have no coordinate from the Ensembl, I filled their coordinates from LDlink results (e.g. ld_coord column). Among the 11 hg19 SNP, `rs78898539` have no information of hg38. Therefore, I found coordinate information of the total 2,003 hg19 SNPs and 2,002 hg38 SNPs. This result was saved as `gwas_snp_2003.bed`.

File name change: `db_gwas/gwas_biomart.bed` -> `db_gwas/gwas_biomart_fill.bed`

For further analysis, you have to generate **BED** format files for both **hg19** and **hg38**.

Usage: `Rscript postgwas-exe.r --ldlink bed --base <input file> <input folder> --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
    --ldlink bed ^
    --base r2d1_data_gwas/gwas_biomart_fill.tsv ^
    --out r2d1_data_gwas
```

> ** Run function ldlink_bed... [1] 2003    8
> Write file:     r2d1_data_gwas/gwas_hg19_biomart_2003.bed
> Write file:     r2d1_data_gwas/gwas_hg38_biomart_2002.bed
> Job done: 2020-02-25 23:56:09 for 1.9 sec

# 2. Downloading annotation data

## Roadmap download

Downloading the [127 cell type](https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types)-specific Roadmap **ChmmModels** BED files (**hg19**) to designated out folder. This process takes about ~35 min depending on your internet speed.

Usage: `Rscript postgwas-exe.r --dbdown roadmap --out <out folder>`

```cmd
# Downloading roadmap data takes ~35 min
Rscript postgwas-exe.r ^
	--dbdown roadmap ^
	--out db_gwas/roadmap
```

> ** Run function: db_download.r/roadmap_down...
>   Directory generated: db_gwas/roadmap
> trying URL 'https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E001_25_imputed12marks_dense.bed.gz'
> Content type 'application/x-gzip' length 8930181 bytes (8.5 MB)
> downloaded 8.5 MB
>
> ...
>
> db_gwas/roadmap/E129_25_imputed12marks_dense.bed.gz
>   Convert: db_gwas/roadmap/E129_25_imputed12marks_dense.bed
>
> File write: db_gwas/roadmap/E129_25_imputed12marks_dense.bed.rds
> Job done: 2020-02-27 16:46:38 for 34.3 min

## ENCODE download

Downloading `wgEncodeRegTfbsClusteredV3.bed.gz` file (81 MB) which is regulatory transcription factor binding site (Reg-TFBS) cluster data from ENCODE.

Usage: `Rscript postgwas-exe.r --dbdown encode --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
	--dbdown encode ^
	--out db_gwas/encode
```

> ** Run function: db_download.r/encode_down...
>   Directory generated: db_gwas/encode
> trying URL 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz'
> Content type 'application/x-gzip' length 84986946 bytes (81.0 MB)
> downloaded 81.0 MB
>
> db_gwas/encode/wgEncodeRegTfbsClusteredV3.bed.gz
> Job done: 2020-02-27 17:01:53 for 9.2 sec

## RegulomeDB download

Downloading category scores for SNPs by evidences such as eQTL, TF binding, matched TF motif, matched DNase Footprint, and DNase peak from the RegulomeDB. Default filtering criteria is **score ≥2b**.

Downloading files:

* `RegulomeDB.dbSNP132.Category1.txt.gz` (2 MB)
* `RegulomeDB.dbSNP132.Category2.txt.gz` (39.3 MB)
* (optional) total dataset `RegulomeDB.dbSNP141.txt.gz` (2.8 GB)

Usage: `Rscript postgwas-exe.r --dbdown regulome --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
	--dbdown regulome ^
	--out db_gwas/regulome
```

> ** Run function: db_download.r/regulome_down...
>   Directory exists: db_gwas/regulome
>   Download RegulomeDB data
>
> trying URL 'http://legacy.regulomedb.org/downloads/RegulomeDB.dbSNP132.Category1.txt.gz'
> Content type 'application/gzip' length 2096454 bytes (2.0 MB)
> downloaded 2.0 MB
>
> [1] 39432     5
> File write: db_gwas/regulome/dbSNP132.Category1.txt.gz.rds
> trying URL 'http://legacy.regulomedb.org/downloads/RegulomeDB.dbSNP132.Category2.txt.gz'
> Content type 'application/gzip' length 41253483 bytes (39.3 MB)
> downloaded 39.3 MB
>
> [1] 407796      5
> File write: db_gwas/regulome/dbSNP132.Category2.txt.gz.rds
> Job done: 2020-02-27 17:07:17 for 34.7 sec

## GTEx download

GTEx v8 includes 17,382 samples of 54 tissues from 948 donors. See [README_eQTL_v8.txt](https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/README_eQTL_v8.txt) for description of GTEx file format.

Usage: `Rscript postgwas-exe.r --dbdown gtex --out <out folder>`

```CMD
# Converting gtex data takes ~1 hr
Rscript postgwas-exe.r ^
	--dbdown gtex ^
	--out db_gwas/gtex
```

> ** Run function: db_download.r/gtex_down...
>   Directory exists: db_gwas/gtex
>   Download GTEx data
>
>   db_gwas/gtex/GTEx_Analysis_v8_eQTL.tar
>     File write: db_gwas/gtex/gtex_files.txt
>   db_gwas/gtex/GTEx_Analysis_2017-06-05_v8_lookup_table.txt.gz
>
>   Loading GTEx BED files
>   File reading...
>   (1/49) Adipose_Subcutaneous
>   (2/49) Adipose_Visceral_Omentum
>   ...
>   (48/49) Vagina
>   (49/49) Whole_Blood
>  gte_df.pval_nominal
>  Min.   :0.000e+00
>  1st Qu.:0.000e+00
>  Median :3.168e-07
>  Mean   :2.399e-05
>  3rd Qu.:1.633e-05
>  Max.   :1.759e-03
>   GTEx table, rows= 71478479 cols= 13
>   BED file read complete. Job process: 13.9 min
>   Loading annotation file
>   Annotation file read complete. Job process: 51 min
>   Annotation file, rows= 46569704 cols= 8
>   GTEx annotation, rows= 71478479 cols= 9
>   Saving a compiled RDS file..  db_gwas/gtex/Gtex_Analysis_v8_eQTL_rsid.rds
> Job done: 2020-02-27 19:07:32 for 52.9 min

## lncRNASNP2 download

In lncRNASNP2, 141,353 human lncRNAs and 10,205,295 SNPs were archived.

```CMD
Rscript postgwas-exe.r ^
	--dbdown lncrna ^
	--out db_gwas/lncrna
```

> ** Run function: db_download.r/regulome_down...
>   Directory exists: db_gwas/lncrna
> trying URL 'http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/snps_mod.txt'
> Content type 'text/plain; charset=GBK' length 477785336 bytes (455.7 MB)
> downloaded 455.7 MB
>
> trying URL 'http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncrnas.txt'
> Content type 'text/plain; charset=GBK' length 7005411 bytes (6.7 MB)
> downloaded 6.7 MB
>
> trying URL 'http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncRNA_associated_disease_experiment.txt'
> Content type 'text/plain; charset=GBK' length 31542 bytes (30 KB)
> downloaded 30 KB
>
> [1] 10205295        3
> File write: db_gwas/lncrna/lncRNASNP2_snplist.txt.rds
> [1] 141271      4
> File write: db_gwas/lncrna/lncrnas.txt.rds
> [1] 753   3
> File write: db_gwas/lncrna/lncrna-diseases_experiment.txt.rds
> Job done: 2020-02-27 21:10:59 for 3.2 min

## Ensembl Genes

Downloading Ensembl gene coordinates through biomaRt (hg19)

```CMD
Rscript postgwas-exe.r ^
	--dbdown gene ^
	--out db_gwas ^
	--hg hg19
```

> ** Run function: db_download.r/biomart_gene...
>   BiomaRt table, dim    = [1] 63677     5
>   File write: db_gwas/ensembl_gene_ann_hg19.tsv
>
>   Filtered table, dim   = [1] 57736     5
>   File write: db_gwas/ensembl_gene_hg19.bed
>
> Job done: 2020-02-28 00:38:13 for 6 sec

## Coding exons from UCSC table browser

Go to UCSC browser hompage https://genome.ucsc.edu

1. For downloading hg38 coding exons,
   * Go to **Table Browser** page
   * assembly: **hg38**
   * group: Genes and Genes Predictions
   * track: **GENCODE v32**
   * table: knownGene
   * output format: BED - browser extensible data

2. Then click the "**get output**" button. In the new page,
   * visibility= hide
   * Create one BED record per: **Coding Exons**

3. Then download the file by click the "**get BED**" button.
   * Set the file name as `ucsc_tbBrowser_Gencode_v32_CDS_hg38.bed`.

4. For downloading hg19 coding exons:
   * Go to **Table Browser** page
   * assembly: **hg19**
   * group: Genes and Genes Predictions
   * track: **Ensembl Genes**
   * table: ensGene
   * output format: BED - browser extensible data
5. Then clisk the "**get output**" button. In the new page,
   * visibility= hid
   * Create one BED record per: **Coding Exons**
6. Then download the file by click the "**get BED**" button.
   * Set the file name as `ucsc_tbBrowser_ensGene_CDS_hg19.bed`.

# 3. Filtering the annotations

## Roadmap filter

### Enhancers

Filtering the Roadmap data by "Enhancers" (e.g., 13_EnhA1, 14_EnhA2, 15_EnhAF, 16_EnhW1, 17_EnhW2, 18_EnhAc). This process might take ~3 min.

The roadmap annotation code information is [here](https://egg2.wustl.edu/roadmap/web_portal/imputed.html).

Usage: `Rscript postgwas-exe.r --dbfilt roadmap --base <base folder> --out <out folder> <...>`

```CMD
Rscript postgwas-exe.r ^
	--dbfilt roadmap ^
	--base db_gwas/roadmap ^
	--out db_gwas ^
	--enh TRUE
```

> ** Run function: db_filter.r/roadmap_filt...
>   Reading files..
>     10/129 being processed.
>     ...
>    Error in gzfile(file, "rb") : cannot open the connection
>    In addition: Warning messages:
>    1: package 'plyr' was built under R version 3.6.2
>    2: package 'data.table' was built under R version 3.6.2
> 3: package 'numbers' was built under R version 3.6.2
> 4: In gzfile(file, "rb") :
> cannot open compressed file 'db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
> db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds - file not found.
> Error in gzfile(file, "rb") : cannot open the connection
> In addition: Warning message:
>   In gzfile(file, "rb") :
>   cannot open compressed file 'db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
> db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds - file not found.
>  ...
>  120/129 being processed.
>   Finished reading and filtering 129 files.
>   
>    Write file: db_gwas/roadmap_enh.bed
>    Job done: 2020-02-27 23:57:08 for 2.5 min

### Total codes

Total roadmap annotations are achieved by this code.

```CMD
# This step takes a lot of system memomry (~20 GB) and about ~6 min.
Rscript postgwas-exe.r ^
	--dbfilt roadmap ^
	--base db_gwas/roadmap ^
	--out db_gwas ^
	--enh FALSE
```

> ** Run function: db_filter.r/roadmap_filt...
>   Reading files..
>     10/129 being processed.
>     ...
> Error in gzfile(file, "rb") : cannot open the connection
> In addition: Warning message:
> In gzfile(file, "rb") :
>   cannot open compressed file 'db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
>   db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds - file not found.
> Error in gzfile(file, "rb") : cannot open the connection
> In addition: Warning message:
> In gzfile(file, "rb") :
>   cannot open compressed file 'db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
>   db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds - file not found.
>     ...
>     120/129 being processed.
>   Finished reading and filtering 129 files.
>
> Write file: db_gwas/roadmap_total.bed
> Job done: 2020-02-28 12:00:35 for 5.7 min

### Enhancer in each cell type

To identify cell type-specific enhancers, filtering the enhancer tags by each cell type. See details in `db_gwas/roadmap_metadata.tsv`:

| EID  | GROUP          | STD_NAME                                                     | EDACC_NAME                                                   |
| ---- | -------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| E033 | Blood & T-cell | Primary T cells from cord blood                              | CD3_Primary_Cells_Cord_BI                                    |
| E034 | Blood & T-cell | Primary T cells from peripheral blood                        | CD3_Primary_Cells_Peripheral_UW                              |
| E037 | Blood & T-cell | Primary T helper memory cells from   peripheral blood 2      | CD4_Memory_Primary_Cells                                     |
| E038 | Blood & T-cell | Primary T helper naive cells from peripheral blood           | CD4_Naive_Primary_Cells                                      |
| E039 | Blood & T-cell | Primary T helper naive cells from peripheral blood           | CD4+_CD25-_CD45RA+_Naive_Primary_Cells                       |
| E040 | Blood & T-cell | Primary T helper memory cells from   peripheral blood 1      | CD4+_CD25-_CD45RO+_Memory_Primary_Cells                      |
| E041 | Blood & T-cell | Primary T helper cells PMA-I stimulated                      | CD4+_CD25-_IL17-_PMA-Ionomycin_stimulated_MACS_purified_Th_Primary_Cells |
| E042 | Blood & T-cell | Primary T helper 17 cells PMA-I stimulated                   | CD4+_CD25-_IL17+_PMA-Ionomcyin_stimulated_Th17_Primary_Cells |
| E043 | Blood & T-cell | Primary T helper cells   from peripheral blood               | CD4+_CD25-_Th_Primary_Cells                                  |
| E044 | Blood & T-cell | Primary T regulatory cells   from peripheral blood           | CD4+_CD25+_CD127-_Treg_Primary_Cells                         |
| E045 | Blood & T-cell | Primary T cells effector/memory enriched   from peripheral blood | CD4+_CD25int_CD127+_Tmem_Primary_Cells                       |
| E047 | Blood & T-cell | Primary T CD8+ naive cells from peripheral   blood           | CD8_Naive_Primary_Cells                                      |
| E048 | Blood & T-cell | Primary T CD8+ memory cells from peripheral   blood          | CD8_Memory_Primary_Cells                                     |
| E062 | Blood & T-cell | Primary mononuclear cells   from peripheral blood            | Peripheral_Blood_Mononuclear_Primary_Cells                   |
| E029 | HSC & B-cell   | Primary monocytes from peripheral blood                      | CD14_Primary_Cells                                           |
| E030 | HSC & B-cell   | Primary neutrophils from peripheral blood                    | CD15_Primary_Cells                                           |
| E031 | HSC & B-cell   | Primary B cells from cord blood                              | CD19_Primary_Cells_Cord_BI                                   |
| E032 | HSC & B-cell   | Primary B cells from peripheral blood                        | CD19_Primary_Cells_Peripheral_UW                             |
| E035 | HSC & B-cell   | Primary hematopoietic stem cells                             | CD34_Primary_Cells                                           |
| E036 | HSC & B-cell   | Primary hematopoietic stem cells short term   culture        | CD34_Cultured_Cells                                          |
| E046 | HSC & B-cell   | Primary Natural Killer cells   from peripheral blood         | CD56_Primary_Cells                                           |
| E050 | HSC & B-cell   | Primary hematopoietic stem cells   G-CSF-mobilized Female    | Mobilized_CD34_Primary_Cells_Female                          |
| E051 | HSC & B-cell   | Primary hematopoietic stem cells   G-CSF-mobilized Male      | Mobilized_CD34_Primary_Cells_Male                            |
| E087 | Other          | Pancreatic Islets                                            | Pancreatic_Islets                                            |
| E098 | Other          | Pancreas                                                     | Pancreas                                                     |

```CMD
Rscript postgwas-exe.r ^
	--dbfilt roadmap ^
	--base db_gwas/roadmap ^
	--out db_gwas/roadmap_enh ^
	--enh TRUE ^
	--sep TRUE
```

> ** Run function: db_filter.r/roadmap_filt...
>   Directory generated: db_gwas/roadmap_enh
>   Reading files..
>     10/129 being processed.
>     ...
> Error in gzfile(file, "rb") : cannot open the connection
> In addition: Warning message:
> In gzfile(file, "rb") :
>   cannot open compressed file 'db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
>   db_gwas/roadmap/E060_25_imputed12marks_dense.bed.rds - file not found.
> Error in gzfile(file, "rb") : cannot open the connection
> In addition: Warning message:
> In gzfile(file, "rb") :
>   cannot open compressed file 'db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds', probable reason 'No such file or directory'
>   db_gwas/roadmap/E064_25_imputed12marks_dense.bed.rds - file not found.
>     ...
>     120/129 being processed.
>   Finished processing 129 files.
>
> Job done: 2020-03-03 22:04:51 for 3.1 min

## GTEx filter

Filtering GTEx eQTL data by p-value <5e-8.

```CMD
Rscript postgwas-exe.r ^
	--dbfilt gtex ^
	--base db_gwas/gtex/Gtex_Analysis_v8_eQTL_rsid.rds ^
	--out db_gwas
	--pval 5e-8
```

> ** Run function: db_filter.r/gtex_filt...
>   P-value threshold     = [1] 5e-08
>   GTEx data, dim        = [1] 71478479        9
>  gtex_sig.pval_nominal
>  Min.   :0.000e+00
>  1st Qu.:0.000e+00
>  Median :8.100e-13
>  Mean   :3.518e-09
>  3rd Qu.:9.148e-10
>  Max.   :5.000e-08
>   GTEx <5e-08, dim      =[1] 30613850        9
> Write file: db_gwas/gtex_signif_5e-08.rds
> Job done: 2020-02-28 00:15:35 for 2.9 min

# 4. Distances from the annotations

## For General annotations

Using `bedtools closest` in bash, this process required a lot of resouce (e.g., ~19GB of RAM). Therefore, files were moved to AWS to BNL server which is highly secured. Then I ran the bellow processes:

```bash
./bin/bedtools sort -i 2020_t1d/db/roadmap_enh.bed | ./bin/bedtools closest -d -a 2020_t1d/gwas_hg19_biomart_2003.bed -b stdin > 2020_t1d/data/roadmap_enh.tsv
```

```bash
./bin/bedtools sort -i 2020_t1d/db/roadmap_total.bed | ./bin/bedtools closest -d -a 2020_t1d/gwas_hg19_biomart_2003.bed -b stdin > 2020_t1d/data/roadmap_total.tsv
```

```bash
#./bin/bedtools sort -i 2020_t1d/db/roadmap_087_enh.bed | ./bin/bedtools closest -d -a 2020_t1d/gwas_hg19_biomart_2003.bed -b stdin > 2020_t1d/data/roadmap_087_enh.tsv
```

```bash
./bin/bedtools sort -i 2020_t1d/db/wgEncodeRegTfbsClusteredV3.bed | ./bin/bedtools closest -d -a 2020_t1d/gwas_hg19_biomart_2003.bed -b stdin > 2020_t1d/data/encode_tfbs.tsv
```

```bash
./bin/bedtools sort -i 2020_t1d/db/ensembl_gene_hg19.bed | ./bin/bedtools closest -d -a 2020_t1d/gwas_hg19_biomart_2003.bed -b stdin > 2020_t1d/data/ensGene_hg19.tsv
```

```bash
./bin/bedtools sort -i 2020_t1d/db/ucsc_tbBrowser_ensGene_CDS_hg19.bed | ./bin/bedtools closest -d -a 2020_t1d/gwas_hg19_biomart_2003.bed -b stdin > 2020_t1d/data/ensGene_cds_hg19.tsv
```

Then result files were downloaded from the BNL server to AWS EC2 to local `data_gwas/distance/` folder.

## For roadmap each cell type

Full code was wrote in `db_gwas/roadmap_dist.sh`:

```bash
bedtools sort -i roadmap_enh/roadmap_001_enh.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > roadmap_dist/roadmap_001_enh.tsv
bedtools sort -i roadmap_enh/roadmap_002_enh.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > roadmap_dist/roadmap_002_enh.tsv
...
bedtools sort -i roadmap_enh/roadmap_002_enh.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > roadmap_dist/roadmap_002_enh.tsv
```

## Preparing Meltonlab's β-cell ATAC-seq data

See details in `db_Meltonlab/README.md` file. Original download files are in `db_Meltonlab/UCSC/` folder.

GEO access: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139817

peaks.bed file (mapped on **hg19**) list:

* GSM4171638_SCbeta_ATAC_rep1.peaks.bed
* GSM4171640_Beta_ATAC_rep1.peaks.bed
* GSM4171642_in_vivo-matured_SCbeta_ATAC_rep1.peaks.bed
* GSM4171643_in_vivo-matured_SCbeta_ATAC_rep2.peaks.bed
* GSM4171644_in_vitro-matured_SCbeta_ATAC_12h.peaks.bed
* GSM4171645_in_vitro-matured_SCbeta_ATAC_72h.peaks.bed

Using `bedtools closest` in bash:

```bash
#bedtools closest -d \
#	-a r2d1_data_gwas/gwas_hg19_biomart_2003.bed \
#	-b db_Meltonlab/GSM4171638_SCbeta_ATAC_rep1.peaks.bed \
#	> r2d1_data_gwas/distance/SCbeta_ATAC-2.tsv
```

```bash
#bedtools closest -d\
#	-a r2d1_data_gwas/gwas_hg19_biomart_2003.bed\
#	-b db_Meltonlab/GSM4171640_Beta_ATAC_rep1.peaks.bed\
#	> r2d1_data_gwas/distance/Beta_ATAC.tsv
```

```bash
#bedtools closest -d\
#	-a r2d1_data_gwas/gwas_hg19_biomart_2003.bed\
#	-b db_Meltonlab/GSM4171642_in_vivo-matured_SCbeta_ATAC_rep1.peaks.bed\
#	> r2d1_data_gwas/distance/SCbeta_in_vivo_ATAC_rep1.tsv
```

```bash
#bedtools closest -d\
#	-a r2d1_data_gwas/gwas_hg19_biomart_2003.bed\
#	-b db_Meltonlab/GSM4171642_in_vivo-matured_SCbeta_ATAC_rep2.peaks.bed\
#	> r2d1_data_gwas/distance/SCbeta_in_vivo_ATAC_rep2.tsv
```

```bash
bedtools closest -d\
	-a r2d1_data_gwas/gwas_hg19_biomart_2003.bed\
	-b db_Meltonlab/GSM4171644_in_vitro-matured_SCbeta_ATAC_12h.peaks.bed\
	> r2d1_data_gwas/distance/SCbeta_12h_ATAC.tsv
```

```bash
bedtools closest -d\
	-a r2d1_data_gwas/gwas_hg19_biomart_2003.bed\
	-b db_Meltonlab/GSM4171645_in_vitro-matured_SCbeta_ATAC_72h.peaks.bed\
	> r2d1_data_gwas/distance/SCbeta_72h_ATAC.tsv
```

## For Meltonlab each cell type

Full code was wrote in `db_gwas/meltonlab_dist.sh`:

```bash
bedtools sort -i ./meltonlab/GSM4171636_PP2_ATAC_rep1.peaks.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > ./meltonlab_dist/PP2_ATAC_rep1.tsv
bedtools sort -i ./meltonlab/GSM4171637_EN_ATAC_rep1.peaks.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > ./meltonlab_dist/EN_ATAC_rep1.tsv
...
bedtools sort -i ./meltonlab/GSE140500_SCbeta_TE.byH3K27ac.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > ./meltonlab_dist/SCbeta_TE_enh.tsv
```

## Preparing Yizhou's β-cell ATAC-seq data

This result data mapped on **hg19** version. See details in `db_Yizhou/README.md` file. Using `bedtools intersect`:

```bash
bedtools intersect -u \
	-a YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled.pf_peaks.bed \
	-b YY005-beta-cell-1_R1.nodup.tn5.pr2_pooled.pf_peaks.bed \
	> YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed
```

`YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed` file has pseudo_pooled pr1 data overlapped with pr2 data. And `YY005-beta-cell-1_R1.nodup.tn5_pooled.pf_peaks.bed` file has union peaks from the four samples.

Here, we used the background corrected (peudo) and union (pooled) pr1 peaks overlapped with pr2 peaks.

## Overlapping T1D SNPs with Yizhou's β-cell ATAC-seq data

Using `bedtools closest`:

```bash
bedtools closest -d\
	-a r2d1_data_gwas/gwas_hg19_biomart_2003.bed\
	-b db_Yizhou/YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed\
	> r2d1_data_gwas/distance/SCbeta_YZ_pseudo_pool_ATAC.tsv
```

# 5. Overlapping the annotations

## General annotations overlapping

To identify the Roadmap and Encode annotations and the CDS regions overlapped T1D SNPs, result files of the general annotations from the `bedtools closest` function were used:

* `r2d1_data_gwas/distance/roadmap_enh.tsv`
* `r2d1_data_gwas/distance/roadmap_total.tsv`
* `r2d1_data_gwas/distance/encode_tfbs.tsv`
* `r2d1_data_gwas/distance/ensGene_cds_hg19.tsv`

Usage: `Rscript postgwas-exe.r --dbfilt dist --base <base file> --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
	--dbfilt dist ^
	--base r2d1_data_gwas/distance/roadmap_enh.tsv r2d1_data_gwas/distance/roadmap_total.tsv r2d1_data_gwas/distance/encode_tfbs.tsv r2d1_data_gwas/distance/ensGene_cds_hg19.tsv ^
	--out r2d1_data_gwas
```

> ** Run function: db_filter.r/distance_filt...
> File roadmap_enh... nrow= 20292.. done
>   Annotations occupied by SNPs  = [1] 3884
>   SNPs in annotations           = [1] 688
>   Write file: r2d1_data_gwas/snp_roadmap_enh_688.bed
> Job process: 0.6 sec
>
> File roadmap_total... nrow= 254205.. done
>   Annotations occupied by SNPs  = [1] 25819
>   SNPs in annotations           = [1] 2003
>   Write file: r2d1_data_gwas/snp_roadmap_total_2003.bed
> Job process: 4.9 sec
>
> File encode_tfbs... nrow= 4237.. done
>   Annotations occupied by SNPs  = [1] 2113
>   SNPs in annotations           = [1] 429
>   Write file: r2d1_data_gwas/snp_encode_tfbs_429.bed
> Job process: 0.1 sec
>
> File ensGene_cds_hg19... nrow= 5791.. done
>   Annotations occupied by SNPs  = [1] 58
>   SNPs in annotations           = [1] 61
>   Write file: r2d1_data_gwas/snp_ensGene_cds_hg19_61.bed
> Job process: 0.1 sec
>
> Job done: 2020-02-28 18:31:08 for 5.8 sec

### Enhancer in Pancreas Islets

```CMD
Rscript postgwas-exe.r ^
	--dbfilt dist ^
	--base r2d1_data_gwas/distance/roadmap_087_enh.tsv ^
	--out r2d1_data_gwas
```

> ** Run function: db_filter.r/distance_filt...
> File roadmap_087_enh... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 53
>   SNPs in annotations           = [1] 79
>   Write file: r2d1_data_gwas/snp_roadmap_087_enh_79.bed
> Job process: 0.2 sec
>
> Job done: 2020-03-01 18:23:21 for 0.2 sec

## For Roadmap each cell type

```CMD
Rscript postgwas-exe.r ^
	--dbfilt dist ^
	--base r2d1_data_gwas/roadmap_dist ^
	--out r2d1_data_gwas/roadmap_dist
```

> ** Run function: db_filter.r/distance_filt...
> File roadmap_001_enh... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 57
>   SNPs in annotations           = [1] 77
>   Write file: r2d1_data_gwas/roadmap_dist/snp_roadmap_001_enh_77.bed
> Job process: 0.3 sec
>
> ...
>
> File roadmap_129_enh... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 64
>   SNPs in annotations           = [1] 82
>   Write file: r2d1_data_gwas/roadmap_dist/snp_roadmap_129_enh_82.bed
> Job process: 0 sec
>
> Job done: 2020-03-04 12:42:23 for 6.7 secs

## ATAC-seq data overlapping

To identify ATAC-seq signal overlapped T1D SNPs, result files of the ATAC-seq data from the `bedtools closest` function were used:

* `r2d1_data_gwas/distance/Beta_ATAC.tsv`
* `r2d1_data_gwas/distance/SCbeta_ATAC.tsv`
* `r2d1_data_gwas/distance/SCbeta_12h_ATAC.tsv`
* `r2d1_data_gwas/distance/SCbeta_72h_ATAC.tsv`
* `r2d1_data_gwas/distance/SCbeta_in_vivo_ATAC_rep1.tsv`
* `r2d1_data_gwas/distance/SCbeta_in_vivo_ATAC_rep2.tsv`
* `r2d1_data_gwas/distance/SCbeta_YZ_pseudo_pool_ATAC.tsv`

Usage: `Rscript postgwas-exe.r --dbfilt dist --base <base file> --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
	--dbfilt dist ^
	--base r2d1_data_gwas/distance/Beta_ATAC.tsv r2d1_data_gwas/distance/SCbeta_ATAC.tsv r2d1_data_gwas/distance/SCbeta_12h_ATAC.tsv r2d1_data_gwas/distance/SCbeta_72h_ATAC.tsv r2d1_data_gwas/distance/SCbeta_in_vivo_ATAC_rep1.tsv r2d1_data_gwas/distance/SCbeta_in_vivo_ATAC_rep2.tsv r2d1_data_gwas/distance/SCbeta_YZ_pseudo_pool_ATAC.tsv ^
	--out r2d1_data_gwas
```

> ** Run function: db_filter.r/distance_filt...
> File Beta_ATAC... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 8
>   SNPs in annotations           = [1] 12
>   Write file: r2d1_data_gwas/snp_Beta_ATAC_12.bed
> Job process: 0.4 sec
>
> File SCbeta_ATAC... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 39
>   SNPs in annotations           = [1] 56
>   Write file: r2d1_data_gwas/snp_SCbeta_ATAC_56.bed
> Job process: 0.1 sec
>
> File SCbeta_12h_ATAC... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 30
>   SNPs in annotations           = [1] 40
>   Write file: r2d1_data_gwas/snp_SCbeta_12h_ATAC_40.bed
> Job process: 0.1 sec
>
> File SCbeta_72h_ATAC... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 28
>   SNPs in annotations           = [1] 39
>   Write file: r2d1_data_gwas/snp_SCbeta_72h_ATAC_39.bed
> Job process: 0.1 sec
>
> File SCbeta_in_vivo_ATAC_rep1... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 19
>   SNPs in annotations           = [1] 24
>   Write file: r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep1_24.bed
> Job process: 0.1 sec
>
> File SCbeta_in_vivo_ATAC_rep2... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 19
>   SNPs in annotations           = [1] 24
>   Write file: r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep2_24.bed
> Job process: 0.1 sec
>
> File SCbeta_YZ_pseudo_pool_ATAC... nrow= 3405.. done
>   Annotations occupied by SNPs  = [1] 38
>   SNPs in annotations           = [1] 58
>   Write file: r2d1_data_gwas/snp_SCbeta_YZ_pseudo_pool_ATAC_58.bed
> Job process: 0.1 sec
>
> Job done: 2020-02-28 18:31:44 for 0.8 sec

## For Meltonlab each cell type

```CMD
Rscript postgwas-exe.r ^
	--dbfilt dist ^
	--base r2d1_data_gwas/meltonlab_dist ^
	--out r2d1_data_gwas/meltonlab_dist
```

> ** Run function: db_filter.r/distance_filt...
> File Alpha_ATAC_rep1... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 8
>   SNPs in annotations           = [1] 12
>   Write file: r2d1_data_gwas/snp_Alpha_ATAC_rep1_12.bed
> Job process: 0.3 sec
>
> ...
>
>   20/21 being processed.
> File SCbeta_SE_enh... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 1
>   SNPs in annotations           = [1] 11
>   Write file: r2d1_data_gwas/snp_SCbeta_SE_enh_11.bed
> Job process: 0 sec
>
> File SCbeta_TE_enh... nrow= 2003.. done
>   Annotations occupied by SNPs  = [1] 6
>   SNPs in annotations           = [1] 10
>   Write file: r2d1_data_gwas/snp_SCbeta_TE_enh_10.bed
> Job process: 0 sec
>
> Job done: 2020-03-04 12:36:43 for 1.4 sec

## Regulome overlapping

Although Regulome data is old fashioned, it is still used for annotations. The T1D SNPs overlapped with Regulome high score (≥2b) annotation were identified through this code:

Usage: `Rscript postgwas-exe.r --dbfilt regulome --base <base files> --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
	--dbfilt regulome ^
	--base r2d1_data_gwas/gwas_hg19_biomart_2003.bed ^
	--lncrna db_gwas/regulome ^
	--out r2d1_data_gwas
```

> ** Run function: db_filter.r/regulome_filt...
> GWAS SNPs N   = [1] 2003
> 2 Regulome data load...
>   Read: db_gwas/regulome/dbSNP132.Category1.txt.gz.rds; dim = [1] 39432     5
>   Read: db_gwas/regulome/dbSNP132.Category2.txt.gz.rds; dim = [1] 407796      5
>
>   Regulome score >=2b, SNPs             = [1] 430528
>   Functional motifs (1a~2b - 1f_only)   = [1] 34705
>
>   Regulome >=2b, GWAS SNPs              = [1] 104
>   GWAS SNPs occupied in
>     functional motifs (1a~2b - 1f_only) = [1] 54
>
> Write file: r2d1_data_gwas/regulome_104.tsv
> Write file: r2d1_data_gwas/snp_regulome2b_104.bed
>
> Job done: 2020-03-01 00:54:10 for 7.9 sec

## GTEx overlapping



```CMD
Rscript postgwas-exe.r ^
	--dbfilt gtex_ovl ^
	--base r2d1_data_gwas/gwas_hg19_biomart_2003.bed ^
	--gtex db_gwas/gtex_signif_5e-08.rds ^
	--out r2d1_data_gwas
```

> ** Run function: db_filter.r/gtex_overlap...
> Input GWAS SNPs N = 2003
>   gtex_signif_5e-08.rds, dim    = [1] 30613850        9
>   Overlapped eQTL-gene pairs    = [1] 68822
>   eQTLs N               = [1] 1349
>   Associated eGenes     = [1] 248
>
> Write file: r2d1_data_gwas/gtex_signif_1349.tsv
>
>   GTEx eQTL BED, dim    = [1] 1349    4
>   eQTL SNP N            = [1] 1349
>
> Write file: r2d1_data_gwas/snp_gtex_1349.bed
> Job done: 2020-03-01 14:32:12 for 2.2 min

```CMD
Rscript postgwas-exe.r ^
	--dbfilt gtex_ovl ^
	--base r2d1_data_gwas/gwas_hg19_biomart_2003.bed ^
	--gtex db_gwas/gtex_signif_5e-08.rds ^
	--out r2d1_data_gwas ^
	--tissue Pancreas
```

> ** Run function: db_filter.r/gtex_overlap...
> Input GWAS SNPs N = 2003
>   gtex_signif_5e-08.rds, dim    = [1] 30613850        9
>   Overlapped eQTL-gene pairs    = [1] 68822
>
> [Option] Pancreas, dim = [1] 1482    9
>   eQTLs N               = [1] 493
>   Associated eGenes     = [1] 42
>
> Write file: r2d1_data_gwas/gtex_signif_Pancreas_493.tsv
>
>   GTEx eQTL BED, dim    = [1] 493   4
>   eQTL SNP N            = [1] 493
>
> Write file: r2d1_data_gwas/snp_gtex_Pancreas_493.bed
> Job done: 2020-03-01 14:34:51 for 2.2 min

## lncRNASNP overlapping

Finding input GWAS SNPs occupied at lncRNA open reading frame regions:

```CMD
Rscript postgwas-exe.r ^
	--dbfilt lnc_ovl ^
	--base r2d1_data_gwas/gwas_hg19_biomart_2003.bed ^
	--lncrna db_gwas/lncrna
	--out r2d1_data_gwas
```

> ** Run function: db_filter.r/lncrna_overlap...
> Input GWAS SNPs N = 2003
> 3 lncRNASNP2 data load...
>   Read: db_gwas/lncrna/lncRNASNP2_snplist.txt.rds;              dim = [1] 10205295        3
>   Read: db_gwas/lncrna/lncrnas.txt.rds;                         dim = [1] 141271      4
>   Read: db_gwas/lncrna/lncrna-diseases_experiment.txt.rds;      dim = [1] 753   3
>
> Summary =
>   lncRNA SNPs
> 1     52   88
>
>   Write file: r2d1_data_gwas/snp_lncrnasnp_88.bed
>   Write file: r2d1_data_gwas/lncrnasnp_88.tsv
> Job done: 2020-03-01 15:56:41 for 28.9 sec

# 6. Summary

## Roadmap cell type union list

Roadmap cell type IDs:

* Blood & T-cells: 33-48, 62
* HSC & B-cells: 29-32, 35-36, 46, 50-51
* Pancreas: 87, 98

```CMD

```



## Roadmap / ENCODE

Identifying "Enhancer Tfbs" SNPs by venn analysis of:

* `r2d1_data_gwas/snp_roadmap_enh_688.bed`
* `r2d1_data_gwas/snp_encode_tfbs_429.bed`

```CMD
Rscript postgwas-exe.r ^
	--dbvenn venn ^
	--base r2d1_data_gwas/snp_roadmap_enh_688.bed r2d1_data_gwas/snp_encode_tfbs_429.bed ^
	--out r2d1_data_gwas ^
	--fig r2d1_fig_gwas
```

> ** Run function: db_venn.r/venn...
>   Read: snp_roadmap_enh_688
>   Read: snp_encode_tfbs_429
>
> Figure draw:            r2d1_fig_gwas/venn_2_snp_roadmap_enh_688, snp_encode_tfbs_429.png
> Write TSV file:         r2d1_data_gwas/venn_2_snp_roadmap_enh_688, snp_encode_tfbs_429.tsv
> Write core BED file:    r2d1_data_gwas/snp_core_269.bed
>
> Job done: 2020-03-01 17:35:49 for 0.4 sec

Change the result file name: `snp_core_269.bed` to `snp_roadmap_encode_269.bed`.

## Roadmap Pancreas Islets / ENCODE

Identifying "Pancreas Islets Enhancer Tfbs" SNPs by venn analysis of:

* `r2d1_data_gwas/snp_roadmap_087_enh_79.bed`
* `r2d1_data_gwas/snp_encode_tfbs_429.bed`

```CMD
Rscript postgwas-exe.r ^
	--dbvenn venn ^
	--base r2d1_data_gwas/snp_roadmap_087_enh_79.bed r2d1_data_gwas/snp_encode_tfbs_429.bed ^
	--out r2d1_data_gwas ^
	--fig r2d1_fig_gwas
```

> ** Run function: db_venn.r/venn...
>   Read: snp_roadmap_087_enh_79
>   Read: snp_encode_tfbs_429
>
> Figure draw:            r2d1_fig_gwas/venn_2_snp_roadmap_087_enh_79, snp_encode_tfbs_429.png
> Write TSV file:         r2d1_data_gwas/venn_2_snp_roadmap_087_enh_79, snp_encode_tfbs_429.tsv
> Write core BED file:    r2d1_data_gwas/snp_core_50.bed
>
> Job done: 2020-03-01 18:28:50 for 0.4 sec

Change the result file name: `snp_core_50.bed` to `snp_roadmap_087_encode_50.bed`.

## Roadmap / ENCODE / Regulome

Identifying "high-probability" SNPs by venn analysis of:

* `r2d1_data_gwas/snp_roadmap_enh_688.bed`
* `r2d1_data_gwas/snp_encode_tfbs_429.bed`
* `r2d1_data_gwas/snp_regulome2b_104.bed`

```CMD
Rscript postgwas-exe.r ^
	--dbvenn venn ^
	--base r2d1_data_gwas/snp_roadmap_enh_688.bed r2d1_data_gwas/snp_encode_tfbs_429.bed r2d1_data_gwas/snp_regulome2b_104.bed ^
	--out r2d1_data_gwas ^
	--fig r2d1_fig_gwas
```

> ** Run function: db_venn.r/venn...
>   Read: snp_roadmap_enh_688
>   Read: snp_encode_tfbs_429
>   Read: snp_regulome2b_104
>
> ** Euler fitting... done.
>
> Figure draw:            r2d1_fig_gwas/euler_3_snp_roadmap_enh_688, snp_encode_tfbs_429, snp_regulome2b_104.png
> Write TSV file:         r2d1_data_gwas/venn_3_snp_roadmap_enh_688, snp_encode_tfbs_429, snp_regulome2b_104.tsv
> Write core BED file:    r2d1_data_gwas/snp_core_57.bed
>
> Job done: 2020-03-01 17:40:45 for 0.5 sec

Change the result file name: `snp_core_57.bed` to `snp_core_high_57.bed`

## GTEx / Roadmap_ENCODE / high_probability

Identifying eQTLs associated genes by venn analysis of:

* `r2d1_data_gwas/snp_gtex_1349.bed`
* `r2d1_data_gwas/snp_roadmap_encode_269.bed`
* `r2d1_data_gwas/snp_core_57.bed`

```CMD
Rscript postgwas-exe.r ^
	--dbvenn venn ^
	--base r2d1_data_gwas/snp_gtex_1349.bed r2d1_data_gwas/snp_roadmap_encode_269.bed r2d1_data_gwas/snp_core_57.bed ^
	--out r2d1_data_gwas ^
	--fig r2d1_fig_gwas
```

> ** Run function: db_venn.r/venn...
>   Read: snp_gtex_1349
>   Read: snp_roadmap_encode_269
>   Read: snp_core_57
>
> ** Euler fitting... done.
>
> Figure draw:            r2d1_fig_gwas/euler_3_snp_gtex_1349, snp_roadmap_encode_269, snp_core_57.png
> Write TSV file:         r2d1_data_gwas/venn_3_snp_gtex_1349, snp_roadmap_encode_269, snp_core_57.tsv
> Write core BED file:    r2d1_data_gwas/snp_core_47.bed
>
> Job done: 2020-03-01 17:46:14 for 0.5 sec

Change the result file name: `snp_core_47.bed` to `snp_core_high_eqtl_47.bed`

## Summary table

Roadmap results

* `r2d1_data_gwas/snp_roadmap_enh_688.bed`
* `r2d1_data_gwas/snp_roadmap_087_enh_79.bed`: Pancreas Islets
* `r2d1_data_gwas/snp_roadmap_total_2003.bed`

ENCODE result

* `r2d1_data_gwas/snp_encode_tfbs_429.bed`

RegulomeDB result

* `r2d1_data_gwas/snp_regulome2b_104.bed`

GTEx results

* `r2d1_data_gwas/snp_gtex_1349.bed`
* `r2d1_data_gwas/snp_gtex_Pancreas_493.bed`

lncRNASNP2 result

* `r2d1_data_gwas/snp_lncrnasnp_88.bed`

CDS region

* `r2d1_data_gwas/snp_ensGene_cds_hg19_61.bed`

ATAC-seq results

* `r2d1_data_gwas/snp_Beta_ATAC_12.bed`
* `r2d1_data_gwas/snp_SCbeta_12h_ATAC_40.bed`
* `r2d1_data_gwas/snp_SCbeta_72h_ATAC_39.bed`
* `r2d1_data_gwas/snp_SCbeta_ATAC_56.bed`
* `r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep1_24.bed`
* `r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep2_24.bed`
* `r2d1_data_gwas/snp_SCbeta_YZ_pseudo_pool_ATAC_58.bed`

```CMD
Rscript postgwas-exe.r ^
	--dbvenn venn ^
	--base r2d1_data_gwas/snp_roadmap_enh_688.bed r2d1_data_gwas/snp_roadmap_087_enh_79.bed r2d1_data_gwas/snp_roadmap_total_2003.bed r2d1_data_gwas/snp_encode_tfbs_429.bed r2d1_data_gwas/snp_regulome2b_104.bed r2d1_data_gwas/snp_gtex_1349.bed r2d1_data_gwas/snp_gtex_Pancreas_493.bed r2d1_data_gwas/snp_lncrnasnp_88.bed r2d1_data_gwas/snp_ensGene_cds_hg19_61.bed r2d1_data_gwas/snp_Beta_ATAC_12.bed r2d1_data_gwas/snp_SCbeta_12h_ATAC_40.bed r2d1_data_gwas/snp_SCbeta_72h_ATAC_39.bed r2d1_data_gwas/snp_SCbeta_ATAC_56.bed r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep1_24.bed r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep2_24.bed r2d1_data_gwas/snp_SCbeta_YZ_pseudo_pool_ATAC_58.bed ^
	--out r2d1_data_gwas ^
	--fig r2d1_fig_gwas
```

> ** Run function: db_venn.r/venn...
>   Read: snp_roadmap_enh_688
>   Read: snp_roadmap_087_enh_79
>   Read: snp_roadmap_total_2003
>   Read: snp_encode_tfbs_429
>   Read: snp_regulome2b_104
>   Read: snp_gtex_1349
>   Read: snp_gtex_Pancreas_493
>   Read: snp_lncrnasnp_88
>   Read: snp_ensGene_cds_hg19_61
>   Read: snp_Beta_ATAC_12
>   Read: snp_SCbeta_12h_ATAC_40
>   Read: snp_SCbeta_72h_ATAC_39
>   Read: snp_SCbeta_ATAC_56
>   Read: snp_SCbeta_in_vivo_ATAC_rep1_24
>   Read: snp_SCbeta_in_vivo_ATAC_rep2_24
>   Read: snp_SCbeta_YZ_pseudo_pool_ATAC_58
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
> Write TSV file:         r2d1_data_gwas/venn_16_snp_roadmap_enh_688-snp_SCbeta_YZ_pseudo_pool_ATAC_58.tsv
>
> [Message] If you need a core rsid BED file,
>         please input two or three files.
> Warning message:
> package 'eulerr' was built under R version 3.6.2

# #. DeapSEA results

