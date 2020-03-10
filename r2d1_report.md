# Post-GWAS analysis : Type 1 Diabetes

This is a log file for analyzing "Type I Diabetes Mellitus" (EFO0001359) from GWAS Catalog.

As of 2020-01-20, the GWAS Catalog contains 4,390 publications and 171,674 associations.
GWAS Catalog data is currently mapped to Genome Assembly GRCh38.p13 and dbSNP Build 153.

And LDlink results are returned with **hg19** coordinates.

## Brief result description

![](Fig1.png)



![](Fig2.png)



![](Fig3.png)



## Additional analysis

* De novo motif: TrawlerWeb [http://trawler.erc.monash.edu.au](http://trawler.erc.monash.edu.au/)
  * Dang et al., 2018, BMC Genomics, [10.1186/s12864-018-4630-0](https://doi.org/10.1186/s12864-018-4630-0)
  * BED format motif list entry (hg19/hg38)
* Core transcriptional Regulatory Circuitries: dbCoRC http://dbcorc.cam-su.org/
  * Huang et al., 2018, Nucleic Acids Res, [10.1093/nar/gkx796](https://dx.doi.org/10.1093%2Fnar%2Fgkx796), pmid [28977473](https://www.ncbi.nlm.nih.gov/pubmed/28977473)
  * Gene symbol list entry (gh19)
* LOLA
* Hi-C-coupled multimarker analysis tool: H-MAGMA
  * Sey et al., 2020, Nature Neuroscience, [10.1038/s41593-020-0603-0](https://doi.org/10.1038/s41593-020-0603-0)

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

> ** Run function ldlink_filter...
> Read download files... 152
> Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
>   line 1 did not have 11 elements
>   [ERROR] rs75793288
>   Read LDlink results           = [1] 303054     12
> Filtering by "r2 > 0.6 and Dprime = 1":
>   Filtered data dimension       = [1] 2147    3
>   Excluded no rsid elements     = [1] 25
> Basic summary of LDlink results:
>   **SNP Tier 1                    = 152**
>   **SNP Tier 2                    = 1851**
>   **SNP candidates                = 2003**
>   SNP source annotation table   = [1] 2003    2
>
> Add annotations:
>   LD block annotation... [1] 139
> Search biomart for SNP coordinates:
>   Query SNPs            = [1] 2003
>   Hg19 result table     = [1] 1992    4
>   Hg38 result table     = [1] 2003    4
>   Cytoband annotation... 2004.. done
>   Merged table          = [1] 2004   11
> Write file: r2d1_data_gwas/gwas_biomart.tsv
> Job done: 2020-03-05 01:41:43 for 20.5 sec

## Generating BED files for coordinates from hg19 and hg38

To prepare candidate SNPs, you have to fill the "NA" coordinate information from the **hg19** and **hg38** columns of the "--ldlink fl" function result file (e.g. `gwas_biomart.tsv`) by searching from the Ensembl biomart. As a note, I found hg19 SNP `rs71080526` have two hg38 SNPs, `rs67544786` and `rs1553647310`. Among those two SNPs, I selected `rs67544786` as same for the Ensemble homepage searching result. And I found the 11 SNPs have no coordinate from the Ensembl, I filled their coordinates from LDlink coordinate information (e.g. ld_coord column). Among the 11 hg19 SNP, `rs78898539` have no hg38 information unfortunately. In total, I found coordinate information of the total 2,003 hg19 SNPs and 2,002 hg38 SNPs. This result was saved as `gwas_snp_2003.bed`.

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
bedtools sort -i db/roadmap_enh.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > data/roadmap_enh.tsv # Error: std::bad_alloc <- memory shortage error
```

```bash
#bedtools sort -i db/roadmap_total.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > data/roadmap_total.tsv
```

```bash
#./bin/bedtools sort -i 2020_t1d/db/roadmap_087_enh.bed | ./bin/bedtools closest -d -a 2020_t1d/gwas_hg19_biomart_2003.bed -b stdin > 2020_t1d/data/roadmap_087_enh.tsv
```

```bash
bedtools sort -i db/wgEncodeRegTfbsClusteredV3.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > data/encode_tfbs.tsv
```

```bash
bedtools sort -i db/ensembl_gene_hg19.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > data/ensGene_hg19.tsv
```

```bash
bedtools sort -i db/ucsc_tbBrowser_ensGene_CDS_hg19.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > data/ensGene_cds_hg19.tsv
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

## Preparing β-cell ATAC-seq data

### Distance of Meltonlab's β-cell ATAC-seq data

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
bedtools closest -d
	-a gwas_hg19_biomart_2003.bed
	-b db/GSM4171638_SCbeta_ATAC_rep1.peaks.bed
	> data/SCbeta_ATAC.tsv
```

```bash
bedtools closest -d
	-a gwas_hg19_biomart_2003.bed
	-b db/GSM4171640_Beta_ATAC_rep1.peaks.bed
	> data/Beta_ATAC.tsv
```

```bash
bedtools closest -d
	-a gwas_hg19_biomart_2003.bed
	-b db/GSM4171642_in_vivo-matured_SCbeta_ATAC_rep1.peaks.bed
	> data/SCbeta_in_vivo_ATAC_rep1.tsv
```

```bash
bedtools closest -d
	-a gwas_hg19_biomart_2003.bed
	-b db/GSM4171643_in_vivo-matured_SCbeta_ATAC_rep2.peaks.bed
	> data/SCbeta_in_vivo_ATAC_rep2.tsv
```

```bash
bedtools closest -d
	-a gwas_hg19_biomart_2003.bed
	-b db/GSM4171644_in_vitro-matured_SCbeta_ATAC_12h.peaks.bed
	> data/SCbeta_12h_ATAC.tsv
```

```bash
bedtools closest -d
	-a gwas_hg19_biomart_2003.bed
	-b db/GSM4171645_in_vitro-matured_SCbeta_ATAC_72h.peaks.bed
	> data/SCbeta_72h_ATAC.tsv
```

### Distance of Meltonlab each cell type data

Full code was wrote in `db_gwas/meltonlab_dist.sh`:

```bash
bedtools sort -i ./meltonlab/GSM4171636_PP2_ATAC_rep1.peaks.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > ./meltonlab_dist/PP2_ATAC_rep1.tsv
bedtools sort -i ./meltonlab/GSM4171637_EN_ATAC_rep1.peaks.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > ./meltonlab_dist/EN_ATAC_rep1.tsv
...
bedtools sort -i ./meltonlab/GSE140500_SCbeta_TE.byH3K27ac.bed | bedtools closest -d -a gwas_hg19_biomart_2003.bed -b stdin > ./meltonlab_dist/SCbeta_TE_enh.tsv
```

### Preparing Yizhou's β-cell ATAC-seq data

This result data mapped on **hg19** version. See details in `db_Yizhou/README.md` file. Using `bedtools intersect`:

```bash
bedtools intersect -u
	-a YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled.pf_peaks.bed
	-b YY005-beta-cell-1_R1.nodup.tn5.pr2_pooled.pf_peaks.bed
	> YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed
```

`YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed` file has pseudo_pooled pr1 data overlapped with pr2 data. And `YY005-beta-cell-1_R1.nodup.tn5_pooled.pf_peaks.bed` file has union peaks from the four samples.

Here, we used the background corrected (peudo) and union (pooled) pr1 peaks overlapped with pr2 peaks.

### Distance of Yizhou's β-cell ATAC-seq data

Using `bedtools closest`:

```bash
bedtools closest -d
	-a gwas_hg19_biomart_2003.bed
	-b db/YY005-beta-cell-1_R1.nodup.tn5.pr1_pooled_pseudo_over.bed
	> data/SCbeta_YZ_pseudo_pool_ATAC.tsv
```

## Preparing Enhancer data of blood cells



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

**Enhancer in Pancreas Islets**

```CMD
#Rscript postgwas-exe.r ^
#	--dbfilt dist ^
#	--base r2d1_data_gwas/distance/roadmap_087_enh.tsv ^
#	--out r2d1_data_gwas
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

Stem cell-derived/primary β-cell ATAC-seq and ChIP-seq (H3K4me1 and H3K21ac) enhancer data.

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

## For Pritchardlab each cell type

25 blood immune cell types ATAC-seq data:

```CMD
Rscript postgwas-exe.r
	--dbfilt dist
	--base r2d1_data_gwas/pritchardlab_dist
	--out r2d1_data_gwas/pritchardlab_dist
```

> ...
>
> File 1011.Naive_Teffs.S... nrow= 2004.. done
>   Annotations occupied by SNPs  = [1] 257
>   SNPs in annotations           = [1] 413
>   Write file: r2d1_data_gwas/pritchardlab_dist/snp_1011.Naive_Teffs.S_413.bed
> Job process: 0.1 sec
>
> Job done: 2020-03-08 23:54:36 for 12.6 sec

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

GTEx version 8 includes 17,382 samples, 54 tissues and 948 donors. 

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

### Pancreas

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

### Whole blood

```CMD
Rscript postgwas-exe.r ^
	--dbfilt gtex_ovl ^
	--base r2d1_data_gwas/gwas_hg19_biomart_2003.bed ^
	--gtex db_gwas/gtex_signif_5e-08.rds ^
	--out r2d1_data_gwas ^
	--tissue Whole_Blood
```

> ** Run function: db_filter.r/gtex_overlap...
> Input GWAS SNPs N       = 2003
>   gtex_signif_5e-08.rds, dim    = [1] 30613850        9
>   Overlapped eQTL-gene pairs    = [1] 68822
>
> [Option] Whole_Blood, dim = [1] 2675    9
>   eQTLs N               = [1] 670
>   Associated eGenes     = [1] 72
>
> Write file: r2d1_data_gwas/gtex_signif_Whole_Blood_670.tsv
>
>   GTEx eQTL BED, dim    = [1] 670   4
>   eQTL SNP N            = [1] 670
>
> Write file: r2d1_data_gwas/snp_gtex_Whole_Blood_670.bed
> Job done: 2020-03-08 23:58:09 for 3.1 min

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

# 6. Table summary

## Roadmap union list

Roadmap cell type IDs:

- Move the Blood & T-cell files (IDs: 33-45,47-48, 62) from `roadmap_dist` to `roadmap_dist/blood_t_cells`
- Move the HSC & B-cell files (IDs: 29-32, 35-36, 46, 50-51) from `roadmap_dist` to `roadmap_dist/hsc_b_cells`
- Move the Pancreas files (IDs: 87, 98) from `roadmap_dist` to `roadmap_dist/pancreas`

### Blood & T-cell

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/roadmap_dist/bed_blood_t_cells ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_roadmap_033_enh_275
>   Read 2: snp_roadmap_034_enh_236
>   Read 3: snp_roadmap_035_enh_184
>   Read 4: snp_roadmap_036_enh_173
>   Read 5: snp_roadmap_037_enh_254
>   Read 6: snp_roadmap_038_enh_233
>   Read 7: snp_roadmap_039_enh_229
>   Read 8: snp_roadmap_040_enh_249
>   Read 9: snp_roadmap_041_enh_215
>   Read 10: snp_roadmap_042_enh_232
>   Read 11: snp_roadmap_043_enh_241
>   Read 12: snp_roadmap_044_enh_250
>   Read 13: snp_roadmap_045_enh_239
>   Read 14: snp_roadmap_047_enh_243
>   Read 15: snp_roadmap_048_enh_242
>   Read 16: snp_roadmap_062_enh_247
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Return function: db_venn.r/summ...
>   Returned union list dim       = [1] 420  20
>   Write a BED file: r2d1_data_gwas/snp_union_bed_blood_t_cells_420.bed
> Warning message:
> package 'eulerr' was built under R version 3.6.2

### HSC & B-cell

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/roadmap_dist/bed_hsc_b_cells ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_roadmap_029_enh_163
>   Read 2: snp_roadmap_030_enh_214
>   Read 3: snp_roadmap_031_enh_245
>   Read 4: snp_roadmap_032_enh_246
>   Read 5: snp_roadmap_046_enh_227
>   Read 6: snp_roadmap_050_enh_219
>   Read 7: snp_roadmap_051_enh_153
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Return function: db_venn.r/summ...
>   Returned union list dim       = [1] 363  11
>   Write a BED file: r2d1_data_gwas/snp_union_bed_hsc_b_cells_363.bed
> Warning message:
> package 'eulerr' was built under R version 3.6.2

### Pancreas

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/roadmap_dist/bed_pancreas ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_roadmap_087_enh_79
>   Read 2: snp_roadmap_098_enh_106
>
> [Message] Can't plot Euler plot.
>
> ** Return function: db_venn.r/summ...
>   Returned union list dim       = [1] 134   6
>   Write a BED file: r2d1_data_gwas/snp_union_bed_pancreas_134.bed
> Warning message:
> package 'eulerr' was built under R version 3.6.2

## Pritchardlab union list

25 Human blood cell types:

* 3 B-cells in rest/stimulated conditions:
  * Bulk B
  * Naive B
  * Memory B
* 3 NK-cells in rest/stimulated conditions:
  * Immature NK
  * Mature NK
  * Memory NK
* 10 CD4+ T-cells
  * Effector CD4+ T
  * Follicular T helper
  * Memory T effector
  * Memory T regulatory
  * Naive T regulatory
  * Naive T effector
  * Regulatory T
  * Th1 precursors
  * Th2 precursors
  * Th17 precursors
* 4 CD8+ T-cells
  * CD8+ T
  * Central memory CD8+ T
  * Effector memory CD8+ T
  * Naive CD8 T
* 5 Other
  * Gamma delta T
  * Monocytes
  * Myeloid
  * pDC
  * Plasmablasts

### B-cells

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/B-cells_rest ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 12 Files in the folder: r2d1_data_gwas/pritchardlab_dist/B-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Bulk_B.U_407
>   Read 2: snp_1001.Mem_B.U_407
>   Read 3: snp_1001.Naive_B.U_399
>   Read 4: snp_1002.Bulk_B.U_411
>   Read 5: snp_1002.Mem_B.U_406
>   Read 6: snp_1002.Naive_B.U_411
>   Read 7: snp_1003.Bulk_B.U_410
>   Read 8: snp_1003.Mem_B.U_407
>   Read 9: snp_1003.Naive_B.U_408
>   Read 10: snp_1004.Bulk_B.U_421
>   Read 11: snp_1004.Mem_B.U_426
>   Read 12: snp_1004.Naive_B.U_428
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 431  16
>   Write a BED file: r2d1_data_gwas/snp_union_B-cells_rest_431.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:37:38 for 3.1 sec

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/B-cells_stim ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 10 Files in the folder: r2d1_data_gwas/pritchardlab_dist/B-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Bulk_B.S_409
>   Read 2: snp_1001.Mem_B.S_408
>   Read 3: snp_1001.Naive_B.S_407
>   Read 4: snp_1002.Bulk_B.S_412
>   Read 5: snp_1002.Mem_B.S_411
>   Read 6: snp_1002.Naive_B.S_408
>   Read 7: snp_1003.Bulk_B.S_416
>   Read 8: snp_1003.Mem_B.S_414
>   Read 9: snp_1003.Naive_B.S_413
>   Read 10: snp_1010.Mem_B.S_429
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  14
>   Write a BED file: r2d1_data_gwas/snp_union_B-cells_stim_432.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:38:01 for 3.1 sec

### NK-cells

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/NK-cells_rest ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 15 Files in the folder: r2d1_data_gwas/pritchardlab_dist/NK-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Immature_NK.U_372
>   Read 2: snp_1001.Mature_NK.U_407
>   Read 3: snp_1001.Memory_NK.U_362
>   Read 4: snp_1002.Immature_NK.U_382
>   Read 5: snp_1002.Memory_NK.U_412
>   Read 6: snp_1003.Immature_NK.U_307
>   Read 7: snp_1003.Mature_NK.U_403
>   Read 8: snp_1003.Memory_NK.U_358
>   Read 9: snp_1004.Immature_NK.U_414
>   Read 10: snp_1004.Mature_NK.U_429
>   Read 11: snp_1004.Memory_NK.U_403
>   Read 12: snp_1008.Immature_NK.U_415
>   Read 13: snp_1008.Mature_NK.U_412
>   Read 14: snp_1008.Memory_NK.U_418
>   Read 15: snp_1010.Memory_NK.U_429
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  19
>   Write a BED file: r2d1_data_gwas/snp_union_NK-cells_rest_432.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:40:37 for 3.1 sec

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/NK-cells_stim ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 6 Files in the folder: r2d1_data_gwas/pritchardlab_dist/NK-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Mature_NK.S_408
>   Read 2: snp_1002.Mature_NK.S_409
>   Read 3: snp_1003.Mature_NK.S_415
>   Read 4: snp_1004.Mature_NK.S_429
>   Read 5: snp_1008.Mature_NK.S_414
>   Read 6: snp_1010.Mature_NK.S_430
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  10
>   Write a BED file: r2d1_data_gwas/snp_union_NK-cells_stim_432.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:40:54 for 3 sec

### CD4+ T-cells

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_rest ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 38 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Effector_CD4pos_T.U_407
>   Read 2: snp_1001.Follicular_T_Helper.U_408
>   Read 3: snp_1001.Memory_Teffs.U_387
>   Read 4: snp_1001.Memory_Tregs.U_385
>   Read 5: snp_1001.Naive_Teffs.U_395
>   Read 6: snp_1001.Regulatory_T.U_377
>   Read 7: snp_1001.Th1_precursors.U_397
>   Read 8: snp_1001.Th17_precursors.U_407
>   Read 9: snp_1001.Th2_precursors.U_405
>   Read 10: snp_1002.Effector_CD4pos_T.U_408
>   Read 11: snp_1002.Follicular_T_Helper.U_413
>   Read 12: snp_1002.Memory_Teffs.U_411
>   Read 13: snp_1002.Memory_Tregs.U_412
>   Read 14: snp_1002.Naive_Teffs.U_412
>   Read 15: snp_1002.Regulatory_T.U_413
>   Read 16: snp_1002.Th1_precursors.U_403
>   Read 17: snp_1002.Th17_precursors.U_401
>   Read 18: snp_1002.Th2_precursors.U_413
>   Read 19: snp_1003.Effector_CD4pos_T.U_398
>   Read 20: snp_1003.Follicular_T_Helper.U_369
>   Read 21: snp_1003.Memory_Teffs.U_402
>   Read 22: snp_1003.Memory_Tregs.U_416
>   Read 23: snp_1003.Naive_Teffs.U_406
>   Read 24: snp_1003.Regulatory_T.U_411
>   Read 25: snp_1003.Th1_precursors.U_359
>   Read 26: snp_1003.Th2_precursors.U_415
>   Read 27: snp_1004.Effector_CD4pos_T.U_427
>   Read 28: snp_1004.Follicular_T_Helper.U_420
>   Read 29: snp_1004.Memory_Teffs.U_417
>   Read 30: snp_1004.Memory_Tregs.U_429
>   Read 31: snp_1004.Naive_Teffs.U_425
>   Read 32: snp_1004.Naive_Tregs.U_419
>   Read 33: snp_1004.Regulatory_T.U_412
>   Read 34: snp_1004.Th1_precursors.U_424
>   Read 35: snp_1004.Th17_precursors.U_415
>   Read 36: snp_1004.Th2_precursors.U_408
>   Read 37: snp_1008.Naive_Tregs.U_426
>   Read 38: snp_1010.Follicular_T_Helper.U_430
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  42
>   Write a BED file: r2d1_data_gwas/snp_union_T-CD4-cells_rest_432.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:41:10 for 3.2 sec

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_stim ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 41 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Effector_CD4pos_T.S_410
>   Read 2: snp_1001.Follicular_T_Helper.S_407
>   Read 3: snp_1001.Memory_Teffs.S_406
>   Read 4: snp_1001.Memory_Tregs.S_409
>   Read 5: snp_1001.Naive_Teffs.S_406
>   Read 6: snp_1001.Regulatory_T.S_399
>   Read 7: snp_1001.Th1_precursors.S_402
>   Read 8: snp_1001.Th17_precursors.S_410
>   Read 9: snp_1001.Th2_precursors.S_405
>   Read 10: snp_1002.Effector_CD4pos_T.S_405
>   Read 11: snp_1002.Follicular_T_Helper.S_407
>   Read 12: snp_1002.Gamma_delta_T.S_412
>   Read 13: snp_1002.Memory_Teffs.S_393
>   Read 14: snp_1002.Memory_Tregs.S_409
>   Read 15: snp_1002.Naive_Teffs.S_407
>   Read 16: snp_1002.Regulatory_T.S_409
>   Read 17: snp_1002.Th1_precursors.S_407
>   Read 18: snp_1002.Th17_precursors.S_405
>   Read 19: snp_1002.Th2_precursors.S_407
>   Read 20: snp_1003.Effector_CD4pos_T.S_378
>   Read 21: snp_1003.Follicular_T_Helper.S_400
>   Read 22: snp_1003.Gamma_delta_T.S_403
>   Read 23: snp_1003.Memory_Teffs.S_406
>   Read 24: snp_1003.Memory_Tregs.S_416
>   Read 25: snp_1003.Naive_Teffs.S_400
>   Read 26: snp_1003.Regulatory_T.S_407
>   Read 27: snp_1003.Th1_precursors.S_412
>   Read 28: snp_1003.Th17_precursors.S_393
>   Read 29: snp_1003.Th2_precursors.S_378
>   Read 30: snp_1004.Follicular_T_Helper.S_422
>   Read 31: snp_1004.Gamma_delta_T.S_427
>   Read 32: snp_1004.Memory_Teffs.S_422
>   Read 33: snp_1004.Memory_Tregs.S_422
>   Read 34: snp_1004.Naive_Teffs.S_425
>   Read 35: snp_1004.Naive_Tregs.S_416
>   Read 36: snp_1004.Regulatory_T.S_422
>   Read 37: snp_1004.Th1_precursors.S_423
>   Read 38: snp_1004.Th17_precursors.S_426
>   Read 39: snp_1004.Th2_precursors.S_417
>   Read 40: snp_1010.Naive_Tregs.S_430
>   Read 41: snp_1011.Naive_Teffs.S_413
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  45
>   Write a BED file: r2d1_data_gwas/snp_union_T-CD4-cells_stim_432.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:41:40 for 3.2 sec

### CD8+ T-cells

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_rest ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 16 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.CD8pos_T.U_411
>   Read 2: snp_1001.Central_memory_CD8pos_T.U_406
>   Read 3: snp_1001.Effector_memory_CD8pos_T.U_411
>   Read 4: snp_1001.Naive_CD8_T.U_406
>   Read 5: snp_1002.CD8pos_T.U_413
>   Read 6: snp_1002.Central_memory_CD8pos_T.U_413
>   Read 7: snp_1002.Effector_memory_CD8pos_T.U_407
>   Read 8: snp_1002.Naive_CD8_T.U_411
>   Read 9: snp_1003.CD8pos_T.U_413
>   Read 10: snp_1003.Central_memory_CD8pos_T.U_386
>   Read 11: snp_1003.Effector_memory_CD8pos_T.U_390
>   Read 12: snp_1003.Naive_CD8_T.U_409
>   Read 13: snp_1004.CD8pos_T.U_425
>   Read 14: snp_1004.Central_memory_CD8pos_T.U_423
>   Read 15: snp_1004.Effector_memory_CD8pos_T.U_412
>   Read 16: snp_1004.Naive_CD8_T.U_422
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  20
>   Write a BED file: r2d1_data_gwas/snp_union_T-CD8-cells_rest_432.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:42:06 for 2.8 sec

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_stim ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
> 15 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.CD8pos_T.S_406
>   Read 2: snp_1001.Central_memory_CD8pos_T.S_407
>   Read 3: snp_1001.Effector_memory_CD8pos_T.S_410
>   Read 4: snp_1001.Naive_CD8_T.S_406
>   Read 5: snp_1002.CD8pos_T.S_410
>   Read 6: snp_1002.Central_memory_CD8pos_T.S_401
>   Read 7: snp_1002.Effector_memory_CD8pos_T.S_408
>   Read 8: snp_1002.Naive_CD8_T.S_412
>   Read 9: snp_1003.CD8pos_T.S_407
>   Read 10: snp_1003.Central_memory_CD8pos_T.S_409
>   Read 11: snp_1003.Effector_memory_CD8pos_T.S_381
>   Read 12: snp_1003.Naive_CD8_T.S_402
>   Read 13: snp_1004.Central_memory_CD8pos_T.S_410
>   Read 14: snp_1004.Effector_memory_CD8pos_T.S_420
>   Read 15: snp_1004.Naive_CD8_T.S_429
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  19
>   Write a BED file: r2d1_data_gwas/snp_union_T-CD8-cells_stim_432.bed
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 00:42:29 for 3.2 sec

## Meltonlab & YZ summary

### Stem cell-derived/primary β-cell ATAC-seq union list

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/meltonlab_dist/bed_beta_ATAC ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Beta_ATAC_12
>   Read 2: snp_Beta_ATAC_rep1_5
>   Read 3: snp_SCbeta_12h_ATAC_40
>   Read 4: snp_SCbeta_72h_ATAC_39
>   Read 5: snp_SCbeta_ATAC_56
>   Read 6: snp_SCbeta_ATAC_rep1_56
>   Read 7: snp_SCbeta_in_vivo_ATAC_rep1_24
>   Read 8: snp_SCbeta_in_vivo_ATAC_rep2_24
>   Read 9: snp_SCbeta_YZ_pseudo_pool_ATAC_58
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Return function: db_venn.r/summ...
>   Returned union list dim       = [1] 96 13
>   Write a BED file: r2d1_data_gwas/snp_union_bed_beta_ATAC_96.bed
> Warning message:
> package 'eulerr' was built under R version 3.6.2

### Stem cell-derived/primary β-cell enhancer union list

```CMD
Rscript postgwas-exe.r ^
  --dbvenn summ ^
  --base r2d1_data_gwas/meltonlab_dist/bed_beta_enh ^
  --out r2d1_data_gwas ^
  --uni TRUE
```

> ** Run function: db_venn.r/summ...
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Beta_SE_enh_22
>   Read 2: snp_Beta_TE_enh_35
>   Read 3: snp_SCbeta_SE_enh_11
>   Read 4: snp_SCbeta_TE_enh_10
>
> [Message] Can't plot Euler plot.
>
> ** Return function: db_venn.r/summ...
>   Returned union list dim       = [1] 73  8
>   Write a BED file: r2d1_data_gwas/snp_union_bed_beta_enh_73.bed
> Warning message:
> package 'eulerr' was built under R version 3.6.2

## Generating summary tables

### Base file list in `summary` folder

**CDS region**

- `r2d1_data_gwas/snp_ensGene_cds_hg19_61.bed`

**Roadmap results**

* `r2d1_data_gwas/summary/snp_roadmap_enh_688.bed`
* `r2d1_data_gwas/summary/snp_union_bed_blood_t_cells_420.bed`
* `r2d1_data_gwas/summary/snp_union_bed_hsc_b_cells_363.bed`
* `r2d1_data_gwas/summary/snp_union_bed_pancreas_134.bed`
* `r2d1_data_gwas/summary/snp_roadmap_087_enh_79.bed`: Pancreas Islets
* ~~`r2d1_data_gwas/snp_roadmap_total_2003.bed`~~

**ENCODE result**

* `r2d1_data_gwas/summary/snp_encode_tfbs_429.bed`

**RegulomeDB result**

* `r2d1_data_gwas/summary/snp_regulome2b_104.bed`

**GTEx results**

* `r2d1_data_gwas/snp_gtex_1349.bed`
* `r2d1_data_gwas/snp_gtex_Pancreas_493.bed`

**lncRNASNP2 result**

* `r2d1_data_gwas/snp_lncrnasnp_88.bed`

**ATAC-seq results**

* `r2d1_data_gwas/snp_Beta_ATAC_12.bed`
* `r2d1_data_gwas/snp_SCbeta_12h_ATAC_40.bed`
* `r2d1_data_gwas/snp_SCbeta_72h_ATAC_39.bed`
* `r2d1_data_gwas/snp_SCbeta_ATAC_56.bed`
* `r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep1_24.bed`
* `r2d1_data_gwas/snp_SCbeta_in_vivo_ATAC_rep2_24.bed`
* `r2d1_data_gwas/snp_SCbeta_YZ_pseudo_pool_ATAC_58.bed`

### Annotation file list

* GWAS: `r2d1_data_gwas/gwas_biomart_fill.tsv`
* ENCODE: `r2d1_data_gwas/distance/encode_tfbs.tsv`
* Nearest genes: `r2d1_data_gwas/distance/ensGene_hg19.tsv`
* Gene CDS: `r2d1_data_gwas/distance/ensGene_cds_hg19.tsv`
* GTEx: `r2d1_data_gwas/gtex_signif_1349.tsv`
* lncRNAs: `r2d1_data_gwas/lncrnasnp_88.tsv`

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/summary
	--out r2d1_data_gwas
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
	--ann_encd r2d1_data_gwas/distance/encode_tfbs.tsv
	--ann_near r2d1_data_gwas/distance/ensGene_hg19.tsv
	--ann_cds r2d1_data_gwas/distance/ensGene_cds_hg19.tsv
	--ann_gtex r2d1_data_gwas/gtex_signif_1349.tsv
	--ann_lnc r2d1_data_gwas/lncrnasnp_88.tsv
```

> ** Run function: db_venn.r/summ...
> 23 Files in the folder: r2d1_data_gwas/summary
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_01_ensGene_cds_hg19_61
>   Read 2: snp_02_roadmap_enh_688
>   Read 3: snp_03_roadmap_bed_pancreas_134
>   Read 4: snp_04_roadmap_bed_blood_t_cells_394
>   Read 5: snp_05_roadmap_bed_hsc_b_cells_372
>   Read 6: snp_06_roadmap_087_enh_79
>   Read 7: snp_07_encode_tfbs_429
>   Read 8: snp_08_regulome2b_104
>   Read 9: snp_09-1_Suh_YZ_SCbeta_ATAC_58
>   Read 10: snp_09_melton_bed_beta_ATAC_96
>   Read 11: snp_10_melton_bed_beta_enh_73
>   Read 12: snp_11_pritchard_T-CD4-cells_rest_432
>   Read 13: snp_12_pritchard_T-CD4-cells_stim_432
>   Read 14: snp_13_pritchard_T-CD8-cells_rest_432
>   Read 15: snp_14_pritchard_T-CD8-cells_stim_432
>   Read 16: snp_15_pritchard_NK-cells_rest_432
>   Read 17: snp_16_pritchard_NK-cells_stim_432
>   Read 18: snp_17_pritchard_B-cells_rest_431
>   Read 19: snp_18_pritchard_B-cells_stim_432
>   Read 20: snp_19_gtex_1349
>   Read 21: snp_20_gtex_Pancreas_493
>   Read 22: snp_21_gtex_Whole_Blood_670
>   Read 23: snp_22_lncrnasnp_88
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 1742   27
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   31
>   Write a CSV file: r2d1_data_gwas/summary_gwas.csv
>
>   ENCODE dim    = [1] 4237   13
>   Merge dim     = [1] 2662   26
>   Write a CSV file: r2d1_data_gwas/summary_encode.csv
>
>   Nearest gene dim      = [1] 2315    9
>   Search biomaRt... 220.. Cache found
> 213.. [1] 2315    5
>   CDS dim               = [1] 5791   11
>   Merge dim             = [1] 2424   31
>   Write a CSV file: r2d1_data_gwas/summary_nearest.csv
>
>   GTEx dim      = [1] 68822     9
>   Merge dim     = [1] 68822    30
>   Write a CSV file: r2d1_data_gwas/summary_gtex.csv
>
>   lncRNA dim    = [1] 132   7
>   Merge dim     = [1] 132  30
>   Write a CSV file: r2d1_data_gwas/summary_lncRNA.csv
>
> Job done: 2020-03-09 00:55:14 for 10.9 sec

### Meltonlab GWAS summary

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/meltonlab_dist/bed_beta_ATAC
	--out r2d1_data_gwas/meltonlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 9 Files in the folder: r2d1_data_gwas/meltonlab_dist/bed_beta_ATAC
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Beta_ATAC_12
>   Read 2: snp_Beta_ATAC_rep1_5
>   Read 3: snp_SCbeta_12h_ATAC_40
>   Read 4: snp_SCbeta_72h_ATAC_39
>   Read 5: snp_SCbeta_ATAC_56
>   Read 6: snp_SCbeta_ATAC_rep1_56
>   Read 7: snp_SCbeta_in_vivo_ATAC_rep1_24
>   Read 8: snp_SCbeta_in_vivo_ATAC_rep2_24
>   Read 9: snp_SCbeta_YZ_pseudo_pool_ATAC_58
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 96 13
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   17
>   Write a CSV file: r2d1_data_gwas/meltonlab_dist/summary_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:12:18 for 3 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/meltonlab_dist/bed_beta_enh
	--out r2d1_data_gwas/meltonlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 4 Files in the folder: r2d1_data_gwas/meltonlab_dist/bed_beta_enh
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_Beta_SE_enh_22
>   Read 2: snp_Beta_TE_enh_35
>   Read 3: snp_SCbeta_SE_enh_11
>   Read 4: snp_SCbeta_TE_enh_10
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 73  8
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   12
>   Write a CSV file: r2d1_data_gwas/meltonlab_dist/summary_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:12:36 for 3 sec

### Pritchardlab GWAS summary

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/B-cells_rest
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 12 Files in the folder: r2d1_data_gwas/pritchardlab_dist/B-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Bulk_B.U_407
>   Read 2: snp_1001.Mem_B.U_407
>   Read 3: snp_1001.Naive_B.U_399
>   Read 4: snp_1002.Bulk_B.U_411
>   Read 5: snp_1002.Mem_B.U_406
>   Read 6: snp_1002.Naive_B.U_411
>   Read 7: snp_1003.Bulk_B.U_410
>   Read 8: snp_1003.Mem_B.U_407
>   Read 9: snp_1003.Naive_B.U_408
>   Read 10: snp_1004.Bulk_B.U_421
>   Read 11: snp_1004.Mem_B.U_426
>   Read 12: snp_1004.Naive_B.U_428
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 431  16
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   20
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/B-cells_rest_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:25:34 for 2.8 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/B-cells_stim
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 10 Files in the folder: r2d1_data_gwas/pritchardlab_dist/B-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Bulk_B.S_409
>   Read 2: snp_1001.Mem_B.S_408
>   Read 3: snp_1001.Naive_B.S_407
>   Read 4: snp_1002.Bulk_B.S_412
>   Read 5: snp_1002.Mem_B.S_411
>   Read 6: snp_1002.Naive_B.S_408
>   Read 7: snp_1003.Bulk_B.S_416
>   Read 8: snp_1003.Mem_B.S_414
>   Read 9: snp_1003.Naive_B.S_413
>   Read 10: snp_1010.Mem_B.S_429
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  14
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   18
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/B-cells_stim_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:25:54 for 3.1 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/NK-cells_rest
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 15 Files in the folder: r2d1_data_gwas/pritchardlab_dist/NK-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Immature_NK.U_372
>   Read 2: snp_1001.Mature_NK.U_407
>   Read 3: snp_1001.Memory_NK.U_362
>   Read 4: snp_1002.Immature_NK.U_382
>   Read 5: snp_1002.Memory_NK.U_412
>   Read 6: snp_1003.Immature_NK.U_307
>   Read 7: snp_1003.Mature_NK.U_403
>   Read 8: snp_1003.Memory_NK.U_358
>   Read 9: snp_1004.Immature_NK.U_414
>   Read 10: snp_1004.Mature_NK.U_429
>   Read 11: snp_1004.Memory_NK.U_403
>   Read 12: snp_1008.Immature_NK.U_415
>   Read 13: snp_1008.Mature_NK.U_412
>   Read 14: snp_1008.Memory_NK.U_418
>   Read 15: snp_1010.Memory_NK.U_429
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  19
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   23
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/NK-cells_rest_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:26:11 for 3.1 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/NK-cells_stim
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 6 Files in the folder: r2d1_data_gwas/pritchardlab_dist/NK-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Mature_NK.S_408
>   Read 2: snp_1002.Mature_NK.S_409
>   Read 3: snp_1003.Mature_NK.S_415
>   Read 4: snp_1004.Mature_NK.S_429
>   Read 5: snp_1008.Mature_NK.S_414
>   Read 6: snp_1010.Mature_NK.S_430
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  10
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   14
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/NK-cells_stim_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:26:26 for 3 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_rest
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 38 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Effector_CD4pos_T.U_407
>   Read 2: snp_1001.Follicular_T_Helper.U_408
>   Read 3: snp_1001.Memory_Teffs.U_387
>   Read 4: snp_1001.Memory_Tregs.U_385
>   Read 5: snp_1001.Naive_Teffs.U_395
>   Read 6: snp_1001.Regulatory_T.U_377
>   Read 7: snp_1001.Th1_precursors.U_397
>   Read 8: snp_1001.Th17_precursors.U_407
>   Read 9: snp_1001.Th2_precursors.U_405
>   Read 10: snp_1002.Effector_CD4pos_T.U_408
>   Read 11: snp_1002.Follicular_T_Helper.U_413
>   Read 12: snp_1002.Memory_Teffs.U_411
>   Read 13: snp_1002.Memory_Tregs.U_412
>   Read 14: snp_1002.Naive_Teffs.U_412
>   Read 15: snp_1002.Regulatory_T.U_413
>   Read 16: snp_1002.Th1_precursors.U_403
>   Read 17: snp_1002.Th17_precursors.U_401
>   Read 18: snp_1002.Th2_precursors.U_413
>   Read 19: snp_1003.Effector_CD4pos_T.U_398
>   Read 20: snp_1003.Follicular_T_Helper.U_369
>   Read 21: snp_1003.Memory_Teffs.U_402
>   Read 22: snp_1003.Memory_Tregs.U_416
>   Read 23: snp_1003.Naive_Teffs.U_406
>   Read 24: snp_1003.Regulatory_T.U_411
>   Read 25: snp_1003.Th1_precursors.U_359
>   Read 26: snp_1003.Th2_precursors.U_415
>   Read 27: snp_1004.Effector_CD4pos_T.U_427
>   Read 28: snp_1004.Follicular_T_Helper.U_420
>   Read 29: snp_1004.Memory_Teffs.U_417
>   Read 30: snp_1004.Memory_Tregs.U_429
>   Read 31: snp_1004.Naive_Teffs.U_425
>   Read 32: snp_1004.Naive_Tregs.U_419
>   Read 33: snp_1004.Regulatory_T.U_412
>   Read 34: snp_1004.Th1_precursors.U_424
>   Read 35: snp_1004.Th17_precursors.U_415
>   Read 36: snp_1004.Th2_precursors.U_408
>   Read 37: snp_1008.Naive_Tregs.U_426
>   Read 38: snp_1010.Follicular_T_Helper.U_430
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  42
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   46
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_rest_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:26:41 for 3.4 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_stim
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 41 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.Effector_CD4pos_T.S_410
>   Read 2: snp_1001.Follicular_T_Helper.S_407
>   Read 3: snp_1001.Memory_Teffs.S_406
>   Read 4: snp_1001.Memory_Tregs.S_409
>   Read 5: snp_1001.Naive_Teffs.S_406
>   Read 6: snp_1001.Regulatory_T.S_399
>   Read 7: snp_1001.Th1_precursors.S_402
>   Read 8: snp_1001.Th17_precursors.S_410
>   Read 9: snp_1001.Th2_precursors.S_405
>   Read 10: snp_1002.Effector_CD4pos_T.S_405
>   Read 11: snp_1002.Follicular_T_Helper.S_407
>   Read 12: snp_1002.Gamma_delta_T.S_412
>   Read 13: snp_1002.Memory_Teffs.S_393
>   Read 14: snp_1002.Memory_Tregs.S_409
>   Read 15: snp_1002.Naive_Teffs.S_407
>   Read 16: snp_1002.Regulatory_T.S_409
>   Read 17: snp_1002.Th1_precursors.S_407
>   Read 18: snp_1002.Th17_precursors.S_405
>   Read 19: snp_1002.Th2_precursors.S_407
>   Read 20: snp_1003.Effector_CD4pos_T.S_378
>   Read 21: snp_1003.Follicular_T_Helper.S_400
>   Read 22: snp_1003.Gamma_delta_T.S_403
>   Read 23: snp_1003.Memory_Teffs.S_406
>   Read 24: snp_1003.Memory_Tregs.S_416
>   Read 25: snp_1003.Naive_Teffs.S_400
>   Read 26: snp_1003.Regulatory_T.S_407
>   Read 27: snp_1003.Th1_precursors.S_412
>   Read 28: snp_1003.Th17_precursors.S_393
>   Read 29: snp_1003.Th2_precursors.S_378
>   Read 30: snp_1004.Follicular_T_Helper.S_422
>   Read 31: snp_1004.Gamma_delta_T.S_427
>   Read 32: snp_1004.Memory_Teffs.S_422
>   Read 33: snp_1004.Memory_Tregs.S_422
>   Read 34: snp_1004.Naive_Teffs.S_425
>   Read 35: snp_1004.Naive_Tregs.S_416
>   Read 36: snp_1004.Regulatory_T.S_422
>   Read 37: snp_1004.Th1_precursors.S_423
>   Read 38: snp_1004.Th17_precursors.S_426
>   Read 39: snp_1004.Th2_precursors.S_417
>   Read 40: snp_1010.Naive_Tregs.S_430
>   Read 41: snp_1011.Naive_Teffs.S_413
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  45
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   49
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/T-CD4-cells_stim_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:27:03 for 3.5 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_rest
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 16 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_rest
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.CD8pos_T.U_411
>   Read 2: snp_1001.Central_memory_CD8pos_T.U_406
>   Read 3: snp_1001.Effector_memory_CD8pos_T.U_411
>   Read 4: snp_1001.Naive_CD8_T.U_406
>   Read 5: snp_1002.CD8pos_T.U_413
>   Read 6: snp_1002.Central_memory_CD8pos_T.U_413
>   Read 7: snp_1002.Effector_memory_CD8pos_T.U_407
>   Read 8: snp_1002.Naive_CD8_T.U_411
>   Read 9: snp_1003.CD8pos_T.U_413
>   Read 10: snp_1003.Central_memory_CD8pos_T.U_386
>   Read 11: snp_1003.Effector_memory_CD8pos_T.U_390
>   Read 12: snp_1003.Naive_CD8_T.U_409
>   Read 13: snp_1004.CD8pos_T.U_425
>   Read 14: snp_1004.Central_memory_CD8pos_T.U_423
>   Read 15: snp_1004.Effector_memory_CD8pos_T.U_412
>   Read 16: snp_1004.Naive_CD8_T.U_422
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  20
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   24
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_rest_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:27:29 for 3.4 sec

```CMD
Rscript postgwas-exe.r
	--dbvenn summ
	--base r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_stim
	--out r2d1_data_gwas/pritchardlab_dist
	--uni_save FALSE
	--ann_gwas r2d1_data_gwas/gwas_biomart_fill.tsv
```

> ** Run function: db_venn.r/summ...
> 15 Files in the folder: r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_stim
>
> ** Run function: db_venn.r/venn_bed...
>   Read 1: snp_1001.CD8pos_T.S_406
>   Read 2: snp_1001.Central_memory_CD8pos_T.S_407
>   Read 3: snp_1001.Effector_memory_CD8pos_T.S_410
>   Read 4: snp_1001.Naive_CD8_T.S_406
>   Read 5: snp_1002.CD8pos_T.S_410
>   Read 6: snp_1002.Central_memory_CD8pos_T.S_401
>   Read 7: snp_1002.Effector_memory_CD8pos_T.S_408
>   Read 8: snp_1002.Naive_CD8_T.S_412
>   Read 9: snp_1003.CD8pos_T.S_407
>   Read 10: snp_1003.Central_memory_CD8pos_T.S_409
>   Read 11: snp_1003.Effector_memory_CD8pos_T.S_381
>   Read 12: snp_1003.Naive_CD8_T.S_402
>   Read 13: snp_1004.Central_memory_CD8pos_T.S_410
>   Read 14: snp_1004.Effector_memory_CD8pos_T.S_420
>   Read 15: snp_1004.Naive_CD8_T.S_429
>
> [Message] Can't plot Venn diagram for more than 5 sets.
>
> [Message] Can't plot Euler plot.
>
> ** Back to function: db_venn.r/summ...
>   Returned union list dim       = [1] 432  19
>
>   [PASS] uni_save       = FALSE
>
>   GWAS dim      = [1] 2003   11
>   Merge dim     = [1] 2003   23
>   Write a CSV file: r2d1_data_gwas/pritchardlab_dist/T-CD8-cells_stim_gwas.csv
>   [PASS] Nearest gene summary.
>
> Job done: 2020-03-09 01:27:47 for 3.1 sec

# 7. Venn summary

## Roadmap Enhancer & ENOCDE Tfbs

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

# #. DeapSEA results



# #. TrawlerWeb results



# #. CRC model results

