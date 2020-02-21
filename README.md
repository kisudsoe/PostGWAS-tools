# Type 1 Diabetes Analysis

This is a log file for analyzing "Type I Diabetes Mellitus" (EFO0001359) from GWAS Catalog.

As of 2020-01-20, the GWAS Catalog contains 4,390 publications and 171,674 associations.
GWAS Catalog data is currently mapped to Genome Assembly GRCh38.p13 and dbSNP Build 153.

In this analysis, **hg38** will applied.



## 1. From the seed SNPs (EFO0001359) to candidate risk SNPs

### Generating pivot tables

From the downloaded GWAS Catalog data, basic statistics of the trait are need to know. This function will generate three summary pivot tables.

Usage: `Rscript postgwas-exe.r --gwas <options> --base <input file> --out <out folder>`

```CMD
Rscript postgwas-exe.r ^
	--gwas trait gene study ^
    --base db_gwas/EFO0001359_2020-02-19.tsv ^
    --out data_gwas
```

> ** Run function trait_pivot:
> Write TRAITS pivot:     data_gwas/EFO0001359_2020-02-19_snps.tsv
>
> ** Run function gene_pivot:
> Write gene pivot:       data_gwas/EFO0001359_2020-02-19_genes.tsv
>
> ** Run function study_pivot:
> Write study pivot:      data_gwas/EFO0001359_2020-02-19_studies.tsv



### Filtering the GWAS SNPs

Using the downloaded GWAS Catalog file, SNPs are needed to filter by the P-value criteria.

Usage: `Rscript postgwas-exe.r --gwas <option> --base <input file> --out <out folder> --p.criteria <number>`

```CMD
Rscript postgwas-exe.r ^
	--gwas filter ^
    --base db_gwas/EFO0001359_2020-02-19.tsv ^
    --out data_gwas ^
    --p.criteria 5e-8
```

> ** Run function gwas_filt:
>   gwas dim              = [1] 330  38
>   gwas (5e-8)           = [1] 240  38
> Write gwas filter:      data_gwas/gwas_5e-8_198.tsv



### Download LDlink data

To expand the seed SNPs to the neighbor SNPs within the LD block, LS associated SNPs from the seed SNPs are downloaded from the LDlink DB (Machiela et al, 2015, Bioinformatics, pmid 26139635).

Usage: `Rscript postgwas-exe.r --ldlink dn --base <input file> --out <out folder> --popul <options>`

```CMD
# Download LD SNPs from 152 query SNPs takes about ~1.2 hrs
Rscript postgwas-exe.r ^
    --ldlink dn ^
    --base db_gwas/gwas_5e-08_152.tsv ^
    --out db_gwas/ldlink ^
    --popul CEU TSI FIN GBR IBS
```

> ** Run function ldlink_dn... 152.. done
>   Files are moved to target folder: db_gwas/ldlink
> Job done: 2020-02-21 15:05:45 for 1.2 hr
> Warning message:
> package 'LDlinkR' was built under R version 3.6.2

### Gathering and filtering LDlink data



Usage: `Rscript postgwas-exe.r --ldlink fl --base <input file> <input folder> --out <out folder> --r2d <option>`

```CMD
Rscript postgwas-exe.r ^
	--ldlink fl ^
    --base db_gwas/gwas_5e-08_152.tsv db_gwas/ldlink ^
    --out db_gwas/ldlink ^
    --r2d 1
```

> Read download files... 152
> Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
>   line 1 did not have 11 elements
> In addition: Warning message:
> package 'LDlinkR' was built under R version 3.6.2
>   rs75793288
>   Read LDlink results   = [1] 303054     12
> Filtering by "r2 > 0.6 and Dprime = 1":
>   Filtered data dimension       = [1] 2135    2
>   Excluded no rsid elements     = [1] 13
> Basic summary of LDlink results:
>   SNP Tier 1    = 152
>   SNP Tier 2    = 1851
>   SNP candidates        = 2003
> Search biomart for hg38 coordinates:
>   Query SNPs    = [1] 2003
>   Result table  = [1] 1954    4
>   Merged table  = [1] 2003    4
> Write file:     db_gwas/ldlink/gwas_biomart.bed
> Job done: 2020-02-21 16:21:07 for 9.4 sec



To prepare candidate SNPs, you have to fill coordinate information (**hg38**) of the 48 SNPs by search in Ensembl biomart SNP synonym. As a note, I found hg19 SNP `rs71080526` have two hg38 SNPs, `rs67544786` and `rs1553647310`. And another hg19 SNP `rs78898539` have no hg38 version. Therefore, I found coordinate information of the total 2,003 hg38 SNPs, and saved as `gwas_snp_2003.bed`.

File name change: `gwas_biomart.bed` -> `gwas_snp_2003.bed`



## 2. Overlapping with ATAC-seq data

See details in `db_ATAC-seq/Source info.md`

