# Type 1 Diabetes Analysis

This is a log file for analyzing "Type I Diabetes Mellitus" (EFO0001359) from GWAS Catalog.

As of 2020-01-20, the GWAS Catalog contains 4,390 publications and 171,674 associations.
GWAS Catalog data is currently mapped to Genome Assembly GRCh38.p13 and dbSNP Build 153.

And LDlink results are returned with **hg19** coordinates.



## 1. From the seed SNPs (EFO0001359) to candidate risk SNPs

### Generating pivot tables

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



### Filtering the GWAS SNPs

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



### Downloading LDlink data

To expand the seed SNPs to the neighbor SNPs within the LD block, LS associated SNPs from the seed SNPs are downloaded from the LDlink DB (Machiela et al, 2015, Bioinformatics, pmid 26139635).

Usage: `Rscript postgwas-exe.r --ldlink dn --base <input file> --out <out folder> --popul <options>`

```CMD
# Download LD SNPs from 152 query SNPs takes about ~1.2 hrs
Rscript postgwas-exe.r ^
    --ldlink dn ^
    --base db_db/gwas_5e-08_152.tsv ^
    --out db_gwas/ldlink ^
    --popul CEU TSI FIN GBR IBS
```

> ** Run function ldlink_dn... 152.. done
>   Files are moved to target folder: db_gwas/ldlink
> Job done: 2020-02-21 15:05:45 for 1.2 hr
> Warning message:
> package 'LDlinkR' was built under R version 3.6.2

### Filtering the LDlink data

Filtering the LDlink data by the criteria, r<sup>2</sup> > 0.6 and D' = 1.

Usage: `Rscript postgwas-exe.r --ldlink fl --base <input file> <input folder> --out <out folder> --r2d <option>`

```CMD
Rscript postgwas-exe.r ^
    --ldlink fl ^
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

### Generating BED files for coordinates from hg19 and hg38

To prepare candidate SNPs, you have to fill the "NA" coordinate information from the **hg19** and **hg38** columns of the "--ldlink fl" function result file (e.g. `gwas_biomart.tsv`) by searching from the Ensembl biomart. As a note, I found hg19 SNP `rs71080526` have two hg38 SNPs, `rs67544786` and `rs1553647310`. Among those two SNPs, I selected `rs67544786` as same for the Ensemble homepage searching result. And I found the 11 SNPs have no coordinate from the Ensembl, I filled their coordinates from LDlink results (e.g. ld_coord column). Among the 11 hg19 SNP, `rs78898539` have no information of hg38. Therefore, I found coordinate information of the total 2,003 hg19 SNPs and 2,002 hg38 SNPs. This result was saved as `gwas_snp_2003.bed`.

File name change: `db_gwas/gwas_biomart.bed` -> `db_gwas/gwas_biomart_fill.bed`

For further analysis, you have to generate **BED** format files for both **hg19** and **hg38**.

Usage: `Rscript postgwas-exe.r --ldlink bed --base <input file> <input folder> --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
    --ldlink bed ^
    --base r2d1_data_gwas/gwas_biomart.tsv ^
    --out r2d1_data_gwas
```

> ** Run function ldlink_bed... [1] 2003    8
> Write file:     r2d1_data_gwas/gwas_hg19_biomart_2003.bed
> Write file:     r2d1_data_gwas/gwas_hg38_biomart_2002.bed
> Job done: 2020-02-25 23:56:09 for 1.9 sec



## 2. Overlapping with ATAC-seq data

See details in `db_ATAC-seq/Source info.md`

