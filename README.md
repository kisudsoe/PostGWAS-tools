# Post GWAS tools

```CMD
Version: 2020-03-06

Usage: Rscript postgwas-exe.r <Function calls>
    <--base file(s)> <--out folder> <options> <--debug>
    
Function calls:
    --gwas    This is a function call for GWAS Catalog data.
    --ldlink  This is a function call for LDlink data.
    --dbdown  This is a function call for downloading databases.
    --filter  This is a function call for filtering data.
    --dbvenn  This is a function call for venn analysis.

Global arguments:
    --base    Base input file is mendatory.
    --out     Out folder path is mendatory.
    --debug   TRUE/FALSE: Rich description for debugging. Default is FALSE.

Running functions with "--help" argument prints [Function] usage information.
```



## --gwas

```CMD
gwas_catalog, v2020-01-06
This is a function for GWAS Catalog data.

Usage: Rscript postgwas-exe.r --gwas <functions> --base <base file> --out <out folder> <...>

Functions:
    tarit   Generating pivot table for traits
    gene    Generating pivot table for genes
    study   Generating summary table for studies
    filt    Filtering SNPs by P-values

Global arguments:
    --base  <EFO0001359.tsv>
            One base TSV file is mendatory.
    --out   <Default:data folder>
            Out files target path is mendatory. Default is "data" folder.

Required arguments:
    --p.criteria 5e-8
            An argument for "--pivot filt". Default is 5e-8.
```



## --ldlink

```CMD
gwas_ldlink, v2020-01-21
This is a function for LDlink data.

Usage: Rscript postgwas-exe.r --ldlink <Function> --base <base file> --out <out folder> <...>
    --ldlink <Functions: dn/fl>

Functions:
    down     This is a function for LDlink data download.
    filter   This is a function for LDlink data filter.
    bed      This is a function for generating two BED files (hg19 and hg38).

Global arguments:
    --base   <EFO0001359.tsv>
             One base TSV file is mendatory.
    --out    <default: data>
             Out folder path is mendatory. Default is "data" folder.

Required arguments:
    --popul  <CEU TSI FIN GBR IBS ...>
             An argument for the "--ldlink dn". One or more population option have to be included.
    --r2d    <1/2/3/4>
             An argument for the "--ldlink fl". Choose one number among these options:
                1) r2>0.6 and Dprime=1  <- The most stringent criteria.
                2) r2>0.6               <- Usual choice to define LD association.
                3) Dprime=1
                4) r2>0.6 or Dprime=1
```



## --dbdown

```
db_download, v2020-02-27
This is a function call for downloading databases
    Roadmap, ENCODE, RegulomeDB, GTEx v8, and lncRNASNP2

Usage: Rscript postgwas-exe.r --dbdown <function> --out <out folder>

Functions:
    roadmap   Downloading Roadmap data (hg19).
    encode    Downloading ENCODE data (hg19).
    regulome  Downloading RegulomeDB data (≥2b).
    gtex      Downloading GTEx v8 data.
    lncrna    Downloading lncRNASNP2 data.
    genes     Downloading Ensembl Biomart Gene coordinates (hg19/hg38).

Global arguments:
    --out     <out folder>
              Download folder path is mendatory. Default is "db" folder.

Required arguments:
    --hg      <hg19/hg38>
              A required argument for the "genes" function. Choose one human genome version.
```



## --dbfilt

```CMD
db_filter, v2020-03-10
This is a function call for filtering data.

Usage: Rscript postgwas-exe.r --dbfilt <function> --base <base file(s)> --out <out folder> <...>

Functions:
    roadmap   Filtering Roadmap data by enhancer tags.
    gtex      Filtering GTEx data by eQTL p-value.
    gtex_ovl  Overlapping the GTEx data with the input GWAS SNPs.
    dist      Filtering distance data from Bedtools closest function.
    regulome  Filtering and overlapping by Regulome score ≥2b.
    lnc_ovl   Overlapping the lncRNASNP2 data with the input GWAS SNPs.

Global arguments:
    --base    <base file/folder>
              Base file/folder path is mendatory.
    --out     <out folder>
              Out folder path is mendatory. Default is "db" folder.

Required arguments:
    --ctype   <cell type id>
              An optional argument for the "roadmap" function.
              See cell-type number information at
              https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types.
    --enh     <default: TRUE>
              An optional argument for the "roadmap" function to filter enhancer regions.
    --sep     <default: FALSE>
              An optional argument for the "roadmap" function to generate cell-type seperated results.
    --meta    <roadmap meta file path>
              An optional argument for the "roadmap", "dist" function.
              For "roadmap" function, this argument needs "--sep TRUE" argument.
              Output file will be organized by the cell types.
    --pval    <p-value threshold>
              A required argument for the "gtex" function to filter significant eQTLs.
    --gtex    <Filtered GTEx RDS file path>
              A required argument for the "gtex_ovl" function to overlap GWAS SNPs with GTEx data.
    --tissue  <GTEx tissue name>
              An optional argument for the "gtex_ovl" function to filter a specific tissue.
    --regulm  <Regulome data folder>
              A required argument for the "regulome" function to load the data.
    --lncrna  <lncRNASNP2 data folder>
```



## --dbvenn

```CMD
db_venn, v2020-03-06
This is a function call for venn analysis of filtered DB data.

Usage: Rscript postgwas-exe.r --dbvenn <function> --base <base files> --out <out folder> --fig <figure out folder>

Functions:
	venn        Venn analysis of rsids.
	summ        Generating summary table with annotations.

Global arguments:
	--base      <base files>
			    At least 2 Base BED file paths are mendatory.
			    If the base file number is over 4, then venn diagram is not generated.
	--out       <out folder>
			    Out folder path is mendatory. Default is "data" folder.

Required arguments:
	--fig       <figure out folder>
			    An optional argument for the "venn" function to save figure file.
			    If no figure out path is designated, no venn figure will generated.
	--uni_list  <default: FALSE>
				An optional argument for the "venn" function to return the union SNP list.
	--dir_only  <default: FALSE>
				An optional argument for the "summ" function to get file paths only in the subfoler.
	--uni_save  <default: TRUE>
				An optional argument for the "summ" function to save to union SNP list as a BED file.
	--ann_gwas  <GWAS annotation TSV file>
			    An optional argument for the "summ" function. Add GWAS annotations to the summary table 1.
	--ann_encd  <ENCODE annotation dist file>
			    An optional argument for the "summ" function. Add ENCODE annotations to the summary table 2.
	--ann_near  <Nearest gene annotation dist file>
			    An optional argument for the "summ" function. Add nearest gene annotations to the summary table 3.
	--ann_cds   <Gene CDS annotation dist file>
			    An optional argument for the "summ" function. Add gene CDS annotations to the summary table 3.
	--ann_gtex  <GTEx eQTL annotation TSV file>
			    An optional argument for the "summ" function. Add GTEx eQTL annotations to the summary table 4.
	--ann_lnc   <lncRNA annotation TSV file>
			    An optional argument for the "summ" function. Add lncRNA annotations to the summary table 5.
```

