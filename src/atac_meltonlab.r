help_messages = function() {
    help_message = '
atac_meltonlab, v2020-02-26
This is a function for ATAC-seq data.

Usgae: Rscript postgwas-exe.r --atac <functions> --base <base file>
    --atac  <Functions: bw2bed>

Functions:
    bw2bed  Converting BigWig to BED format.

Global arguments:
    --base  <folder path>

Required arguments:
    none
' %>% cat
    quit()
}

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
bw2bed_multi = function(
    folder_paths = NULL,  # Folder path(s).
    debug        = F
) {
    # Function specific library
    suppressMessages(library(rtracklayer))

    paste0('\n** Run function bw2bed...\n') %>% cat
    for(i in 1:length(folder_paths)) {
        source('src/pdtime.r'); t0=Sys.time()
        bw2bed(folder_paths[i],debug)
        paste0(pdtime(t0,2),'\n') %>% cat
    }
}

bw2bed = function(
    folder_path = NULL, # Folder path (single).
    debug       = F
) {
    paste0('Folder ',folder_path,'... ') %>% cat
    fs = list.files(folder_path)
    length(fs) %>% print

    # Generate names list
    names = lapply(fs,function(x) {
        x1 = strsplit(x,"_") %>% unlist
        which_x = grep("UM|merged|normalized",x1)
        x1 = x1[!x1 %in% x1[which_x]]
        x2 = paste0(x1,collapse="_")
        return(x2)
    }) %>% unlist

    # Convert format from bw to bed
    for(i in 1:length(names)) {
        paste0('  Generate BED table of ',names[i],'.. ') %>% cat
        f_bw = import.bw(paste0(folder_path,'/',fs[i]))
        n = length(f_bw@seqnames)
        name = paste0(names[i],'_signal',c(1:n))

        # Extract data for BED format
        paste0('with dim = ') %>% cat
        f_bed1 = data.frame(
            chr   = f_bw@seqnames,
            f_bw@ranges,
            score = f_bw@elementMetadata
        )
        f_bed2 = cbind(f_bed1,name)
        f_bed = f_bed2[,c(1:3,6,5)] # Reorder columns
        f_name = paste0(folder_path,'/',names[i],'.bed')
        dim(f_bed) %>% print
        write.table(f_bed,f_name,row.names=F,col.names=F,quote=F,sep='\t')
        paste0('    Write bed file: ',f_name,'\n') %>% cat
    }
}

atac_meltonlab = function(
    args = args
) {
    if(length(args$help)>0)    help = args$help
    if(help)                   help_messages()

    if(length(args$base)>0)    b_path = args$base
    if(length(args$debug)>0) { debug = args$debug
    } else                     debug = FALSE

    source('src/pdtime.r'); t0=Sys.time()
    if(args$atac == 'bw2bed') {
        bw2bed_multi(b_path,debug)
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}