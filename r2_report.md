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
Rscript postgwas-exe.r
  --ldlink filter
  --base db_gwas/gwas_5e-08_152.tsv db_gwas/ldlink
  --out r2_data
  --r2d 2
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
>   SNP Tier 1                    = 152
>   SNP Tier 2                    = 5738
>   SNP candidates                = 5890
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