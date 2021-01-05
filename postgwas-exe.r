#!/usr/bin/env Rscript

# Help Messages ----------
help_message = '
Version: 2020-07-22

Usage:
    Rscript postgwas-exe.r --gene   <Function>
    Rscript postgwas-exe.r --gwas   <Function> --base <base file(s)> --out <out folder> [options:--p.criteria]
    Rscript postgwas-exe.r --ldlink <Function> --base <base file(s)> --out <out folder> [options:--popul --r2 --dprime]
    Rsciprt postgwas-exe.r --utils  <Function> --base <base file>    --out <out folder>
    Rscript postgwas-exe.r --dbdown <Function> --out  <out folder> [option:--hg]
    Rscript postgwas-exe.r --dbfilt <Function> --base <base file/folder> --out <out folder> [option:--hg]
    Rscript postgwas-exe.r --dbvenn <Function> --base <base file/folder> --out <out folder> [optioons:--sub_dir,--uni_save,--ann_gwas,--ann_encd,--ann_near,--ann_eqtl]
    ...

    
Function calls:
    --gwas    A function for GWAS Catalog data.
    --ldlink  A function for LDlink data.
    --utils   A function for utils such as gene  generating bedtools bash file.
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
'


# Command Input ----------
commandLine = commandArgs(trailingOnly=T)
comm = paste(unlist(commandLine),collapse=' ')
listoptions = unlist(strsplit(comm,'--'))[-1]
args = lapply(listoptions,function(x) {
    arg = unlist(strsplit(x,' '))[-1]
    if(length(arg)==0) { return(TRUE)
    } else return(arg)
})
args_names = sapply(listoptions,function(x) {
    option = unlist(strsplit(x,' '))[1]
})
args_names = unlist(args_names)
names(args) = args_names


# Run Function ----------
call = args_names[1]
if(call=='help') {
    cat(help_message); quit()
} else if(call=='gwas') {
    source('/src/gwas_catalog.r')
    gwas_catalog(args)
} else if(call=='ldlink') {
    source('src/gwas_ldlink.r')
    gwas_ldlink(args)
} else if(call=='utils') {
    source('/src/utils.r')
    utils(args)
} else if(call=='dbdown') {
    source('/src/db_download.r')
    db_download(args)
} else if(call=='dbfilt') {
    source('/src/db_filter.r')
    db_filter(args)
} else if(call=='dbvenn') {
    source('/src/db_venn.r')
    db_venn(args)
} else if(call=='dbgene') {
    source('/src/db_gene.r')
    db_gene(args)
} else if(call=='dbcomp') {
    source('/src/db_compare.r')
    db_compare(args)
} else {
    error_m = paste0('\n** [ERROR] There is no such function: "',commandLine[1],'". Try again.\n')
    cat(error_m)
}
