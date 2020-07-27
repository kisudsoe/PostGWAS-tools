# Post GWAS tools

```bash
Version: 2020-06-26

Usage:
    Rscript postgwas-exe.r --gwas <Function> --base <base file(s)> --out <out folder> [options:--p.criteria]
    Rscript postgwas-exe.r --ldlink <Function> --base <base file(s)> --out <out folder> [options:--popul --r2 --dprime]
    Rscript postgwas-exe.r --dbdown <Function> --out <out folder> [option:--hg]
    Rscript postgwas-exe.r --dbfilt <Function> --base <base file/folder> --out <out folder> [option:--hg]
    ...

    
Function calls:
    --gwas    A function for GWAS Catalog data.
    --ldlink  A function for LDlink data.
    --dbdown  A function for downloading databases.
    --dbfilt  A function for filtering data.
    --dbvenn  A function for venn analysis.
    --dbgene  A function for gene analysis.
    --dbcomp  A function for PCA analysis for datasets.

Global arguments:
    --base    <Base input file path>
              This is mendatory.
    --out     <Out folder path>
              This is mendatory.
    --debug   <default: FALSE>
              TRUE/FALSE: Rich description for debugging.
              This is optional.

Running functions with "--help" argument prints [Function] usage information.
```



# Quick start

This function is not developed yet.



# Step-specific tools

## --gwas: src/gwas_catalog.r

```bash
gwas_catalog, v2020-06-26
This is a function for GWAS Catalog data.

Usage:
    Rscript postgwas-exe.r --gwas trait --base <base TSV file> --out <out folder>
    Rscript postgwas-exe.r --gwas gene --base <base TSV file> --out <out folder>
    Rscript postgwas-exe.r --gwas study --base <base TSV file> --out <out folder>
    Rscript postgwas-exe.r --filter trait --base <base TSV file> --out <out folder> --p.criteria 5e-8


Functions:
    trait   Generating pivot table for traits.
    gene    Generating pivot table for genes.
    study   Generating summary table for studies.
    filter  Filtering SNPs by P-values.

Global arguments:
    --base  <EFO0001359.tsv>
            An input base TSV file downloaded from GWAS Catalog is mendatory.
    --out   <Default:data folder>
            Out files target directory path that is mendatory. Default is "data" folder.

Required arguments:
    --p.criteria <5e-8>
            An argument for "--pivot filt". Default is 5e-8.
```



## --ldlink: src/gwas_ldlink.r

```bash
gwas_ldlink, v2020-06-26
This is a function for LDlink data.

Usage:
    Rscript postgwas-exe.r --ldlink down --base <base file> --out <out folder> --popul <CEU TSI FIN GBR IBS ...>
    Rscript postgwas-exe.r --ldlink filter --base <base file> <ldlink dir path> --out <out folder> --r2 0.6 --dprime 1
    Rscript postgwas-exe.r --ldlink filter --base <base file> --out <out folder> --r2 0.5
    Rscript postgwas-exe.r --ldlink bed --base <base file> --out <out folder>


Functions:
    down     This is a function for LDlink data download.
    filter   This is a function for LDlink data filter.
    bed      This is a function for generating two BED files (hg19 and hg38).

Global arguments:
    --base   <EFO0001359.tsv>
             One base TSV file is mendatory.
             A TSV file downloaded from GWAS catalog for "down" and "filter" functions.
             A TSV file processed from "filter" function for "bed" function.
    --out    <default: data>
             Out folder path is mendatory. Default is "data" folder.

Required arguments:
    --popul  <CEU TSI FIN GBR IBS ...>
             An argument for the "--ldlink ddown". One or more population option have to be included.
    --r2     An argument for the "--ldlink filter". Set a criteria for r2 over.
    --dprime An argument for the "--ldlink filter". Set a criteria for dprime over.
```



## --dbdown: src/db_download.r

```bash
db_download, v2020-07-22
This is a function call for downloading databases
    Roadmap, ENCODE, RegulomeDB, GTEx v8, and lncRNASNP2

Usage:
    Rscript postgwas-exe.r --dbdown roadmap --out <out folder>
    Rscript postgwas-exe.r --dbdown encode --out <out folder>
    Rscript postgwas-exe.r --dbdown regulome --out <out folder>
    Rscript postgwas-exe.r --dbdown gtex --out <out folder>
    Rscript postgwas-exe.r --dbdown lncrna --out <out folder>
    Rscript postgwas-exe.r --dbdown gene --out <out folder> --hg hg19
    Rscript postgwas-exe.r --dbdown gene --out <out folder> --hg hg38
    Rscript postgwas-exe.r --dbdown genebed --base <Rsid list TSV file path> --out <out folder> --hg hg19
    Rscript postgwas-exe.r --dbdown genebed --base <Rsid list TSV file path> --out <out folder> --hg hg38


Functions:
    roadmap   Downloading Roadmap data (hg19).
    encode    Downloading ENCODE data (hg19).
    regulome  Downloading RegulomeDB data (≥2b, hg19).
    gtex      Downloading GTEx v8 data (hg38).
    lncrna    Downloading lncRNASNP2 data (hg38).
    gene      Downloading Ensembl Biomart Gene coordinates (hg19/hg38).
    genebed   Downloading seed SNP coordinates from biomaRt

Global arguments:
    --out     <out folder>
              Download folder path is mendatory. Default is "db" folder.

Function-specific arguments:
    --base    <Rsid list file path>
    --hg      <hg19/hg38>
              A required argument for the "gene" function. Choose one human genome version.
```



## --bedtools: src/bedtools.r

```bash
bedtools, v2020-07-22
This is a function call for generating bedtools command.

Usage:
    Rscript postgwas-exe.r --bedtools bash --base <base file> --out <out folder>


Function:
    bash    Generating bash command scripts to run bedtools

Global arguments:
    --base  <base file path>
            Mendatory. For bash function.
    --out   <out folder path>
            Mendatory. For bash function.
```



## --dbfilt: src/db_filter.r

```bash
db_filter, v2020-07-23
This is a function call for filtering data.

Usage:
    Rscript postgwas-exe.r --dbfilt ucsc --base <CDS> <Gene> <Promoter> --out <out folder>
    Rscript postgwas-exe.r --dbfilt roadmap --base <base file(s)> --out <out folder> --enh TRUE
    Rscript postgwas-exe.r --dbfilt roadmap --base <base file(s)> --out <out folder> --enh FALSE
    Rscript postgwas-exe.r --dbfilt roadmap --base <base file(s)> --out <out folder> --enh TRUE --sep TRUE


Functions:
    ucsc        Compile the ucsc downloaded promoter/gene/cds region annotations.
    roadmap     Filtering Roadmap data by enhancer tags.
    gtex        Filtering GTEx data by eQTL p-value.
    gtex_ovl    Overlapping the GTEx data with the input GWAS SNPs.
    hic_bed     --base <HiCCUPS files> --out <out folder>
                Converting the HiCCUPS data to the BED format.
    dist        Filtering distance data from Bedtools closest function.
    regulome    Filtering and overlapping by Regulome score ≥2b.
    lnc_ovl     Overlapping the lncRNASNP2 data with the input GWAS SNPs.

Global arguments:
    --base      <base file/folder>
                Base file/folder path is mendatory.
                For ucsc function, you have to input three UCSC downloaded BED files by this order:
                  [1] cds region file, [2] whole gene region file, [3] proximal promoter region file
    --out       <out folder>
                Out folder path is mendatory. Default is "db" folder.

Required arguments:
    --ctype     <cell type id>
                An optional argument for the "roadmap" function.
                See cell-type number information at
                https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types.
    --enh       <default: TRUE>
                An optional argument for the "roadmap" function to filter enhancer regions.
    --sep       <default: FALSE>
                An optional argument for the "roadmap" function to generate cell-type seperated results.
    --meta      <roadmap meta file path>
                An optional argument for the "roadmap", "dist" function.
                For "roadmap" function, this argument needs "--sep TRUE" argument.
                Output file will be organized by the cell types.
    --infotype  <default: FALSE>
                An optional argument for the "dist" function.
                If input as "ucsc", three output files will be generated.
                  [1] cds region, [2] whole gene region, [3] proximal promoter region
                If input as "tags", This option allows to generate seperated output files by the tags.
    --pval      <p-value threshold>
                A required argument for the "gtex" function to filter significant eQTLs.
    --gtex      <Filtered GTEx RDS file path>
                A required argument for the "gtex_ovl" function to overlap GWAS SNPs with GTEx data.
    --tissue    <GTEx tissue name>
                An optional argument for the "gtex_ovl" function to filter a specific tissue.
    --regulm    <Regulome data folder>
                A required argument for the "regulome" function to load the data.
    --lncrna    <lncRNASNP2 data folder>
```



## --dbvenn: src/db_venn.r

```bash
db_venn, v2020-07-22
This is a function call for venn analysis of filtered DB data.

Usage: Rscript postgwas-exe.r --dbvenn <function> --base <base files> --out <out folder> --fig <figure out folder>

Functions:
	venn        Venn analysis of rsids.
	summ        Generating summary table with annotations.

Global arguments:
	--base      <base files>
			    At least 2 Base BED file paths are mendatory.
			    For the venn function,if the base file number is over 4, then venn diagram is not generated.
	--out       <out folder>
			    Out folder path is mendatory. Default is "data" folder.

Required arguments:
	--fig       <figure out folder>
			    An optional argument for the "venn" function to save figure file.
			    If no figure out path is designated, no venn figure will generated.
	--uni_list  <default: FALSE>
				An optional argument for the "venn" function.
				To return the union SNP list.
	--sub_dir   <default: FALSE>
				An optional argument for the "summ" function.
				To get file paths grouped by the subfoler.
	--uni_save  <default: TRUE>
				An optional argument for the "summ" function.
				To save to union SNP list as a BED file.
	--ann_gwas  <GWAS annotation TSV file>
			    An optional argument for the "summ" function.
				Add GWAS annotations to the summary table 1.
	--ann_encd  <ENCODE annotation dist file>
			    An optional argument for the "summ" function.
				Add ENCODE annotations to the summary table 2.
	--ann_near  <Nearest gene annotation dist file>
			    An optional argument for the "summ" function.
				Add nearest gene annotations to the summary table 3.
	--ann_cds   <Gene CDS annotation dist file>
			    An optional argument for the "summ" function.
				Add gene CDS annotations to the summary table 3.
	--ann_gtex  <GTEx eQTL annotation TSV file>
			    An optional argument for the "summ" function.
				Add GTEx eQTL annotations to the summary table 4.
	--ann_lnc   <lncRNA annotation TSV file>
			    An optional argument for the "summ" function.
				Add lncRNA annotations to the summary table 5.
```



## --dbgene: src/db_gene.r

```bash
Version: v2020-04-29
This is a function call for gene analysis to compile the Hi-C and eQTL data.

Usage: Rscript postgwas-exe.r --dbgene <function> --base <base files> --out <out folder> [options]

Functions:
    hic_pair    Extract Hi-C linked SNP-Gene pairs.
    gtex_pair   Extract eQTL linked SNP-Gene pairs.
    summary     Summarizing GWAS SNPs-Gene pairs.
    david_go    Summarizing DAVID GO analysis result.
    pivot_gene  Pivotting the gene summary pair table to gene level summary.

Global arguments:
    --base      <base files>
                For "hic_pair" function, two input files are:
                  [1] gwas_dist file, [2] gene_dist file
    --out       <out folder>
                Out folder path is mendatory. Default is "data" folder.
                If the "bed=TRUE" option put in "hic_pair", out folder could be two.

Required arguments:
    --bed       <default:FALSE>
                An optional argument for the "hic_pair" function.
                To save as BED format file.
    --nearest   <nearest gene summary file path>
                An optional argument for the "summary" function.
                Add nearest gene summary table to the hic gene summary.
    --criteria  <default:0.05>
                An optional argument for the "david_go" function.
    --stat      <default:fdr>
                An optional argument for the "david_go" function.
                Either --fdr or --pval have to choose.
    --dataset   An optional argument for the "david_go" function.
                Add filtering dataset name.
```