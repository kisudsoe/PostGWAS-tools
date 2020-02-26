#!/usr/bin/env Rscript

# Help Messages ----------
help_message = c('
Version: 2020-02-20

Usage: Rscript postgwas-exe.r [Functions]
    <--base file(s)> <--out folder> <options> <--debug>
    
Functions:
    --gwas    This is a function for GWAS Catalog data.
    --ldlink  This is a function for LDlink data.
    --bw2bed  This is a function to convert BigWig to BED format.

Global arguments:
    --base    Base input file is mendatory.
    --out     Out folder path is mendatory.
    --debug   TRUE/FALSE: Rich description for debugging. Default is FALSE.
    
Running functions with --help or -h arguments prints usage information for the [Function].
')

# Command Input ----------
commandLine = commandArgs(trailingOnly=T)
help   = (sum(c('--help')   %in% commandLine[1])>0)
gwas   = (sum(c('--gwas')   %in% commandLine[1])>0)
ldlink = (sum(c('--ldlink') %in% commandLine[1])>0)

comm = paste(unlist(commandLine),collapse=' ')
listoptions = unlist(strsplit(comm,'--'))[-1]
args = lapply(listoptions,function(x) {
    arg = unlist(strsplit(x,' '))[-1]
    if(length(arg)==0) { return(TRUE)
    } else return(arg)
})
args.names = sapply(listoptions,function(x) {
    option = unlist(strsplit(x,' '))[1]
})
names(args) = unlist(args.names)

# Run Function ----------
if(help) { cat(help_message)
} else {
    if(gwas) {
        source('src/gwas_catalog.r')
        gwas_catalog(args)
    } else if (ldlink) {
        source('src/gwas_ldlink.r')
        gwas_ldlink(args)
    } else {
        paste0('There is no such function: "',commandLine[1],'". Try again.')
    }
}