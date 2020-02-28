#!/usr/bin/env Rscript

# Help Messages ----------
help_message = c('
Version: 2020-02-20

Usage: Rscript postgwas-exe.r <Function calls>
    <--base file(s)> <--out folder> <options> <--debug>
    
Function calls:
    --gwas    This is a function call for GWAS Catalog data.
    --ldlink  This is a function call for LDlink data.
    --dbdown  This is a function call for downloading databases.
    --filter  This is a function call for filtering data.

Global arguments:
    --base    Base input file is mendatory.
    --out     Out folder path is mendatory.
    --debug   TRUE/FALSE: Rich description for debugging. Default is FALSE.
    
Running functions with --help or -h arguments prints usage information for the [Function].
')

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
if(call=='help') { cat(help_message)
} else if(call=='gwas') {
    source('src/gwas_catalog.r')
    gwas_catalog(args)
} else if(call=='ldlink') {
    source('src/gwas_ldlink.r')
    gwas_ldlink(args)
} else if(call=='dbdown') {
    source('src/db_download.r')
    db_download(args)
} else if(call=='dbfilt') {
    source('src/db_filter.r')
    db_filter(args)
} else {
    paste0('There is no such function: "',commandLine[1],'". Try again.')
}