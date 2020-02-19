## Version ##
# 1.0 written by Seungsoo Kim for Post-GWAS analysis pipeline

## Usage ##
# gwassnpsummary: a function for generating pivot tables of GWAS Catalog data
# Usage: gwassnpsummary(gwas, out, debug=F)
#   gwas = 'db_gwas/EFO0001359_2020-02-19.tsv' # GWAS Catalog file path
#   out  = 'data_gwas'                         # out files target path

library(dplyr)
library(tidyr)
gwassnpsummary = function(
    gwas  = NULL,   # GAWS Catalog file address
    out   = 'data', # Output folder path
    debug = F
) {
    # Read GWAS Catalog file
    gdata = read.delim(gwas)
    file_nm = tools::file_path_sans_ext(gwas %>% basename)
    if(debug) {
        paste0('Read file, ',basename(gwas),'\t= ') %>% cat
        dim(gdata) %>% print
    }
    
    # Generate unique TRAITS
    TRAITS = paste0(gdata$DISEASE.TRAIT,"_",gdata$PUBMEDID)
    if(debug) {
        paste0('  TRAITS length\t= ') %>% cat
        print(length(TRAITS))
    }
    
    # Shrink data
    gdata2 = data.frame(
        SNPS        = gdata$SNPS,
        MAPPED_GENE = gdata$MAPPED_GENE,
        TRAITS      = TRAITS,
        P.VALUE     = gdata$P.VALUE
    ) %>% unique
    if(debug) {
        paste0('  gdata2 dim\t= ') %>% cat
        dim(gdata2) %>% print
    }
    
    # Reshape data 1: Spread by TRAITS
    pivot1 = gdata2 %>%
        group_by(SNPS,MAPPED_GENE,TRAITS) %>%
        summarize(Min_P = min(P.VALUE)) %>%
        spread('TRAITS','Min_P')
    if(debug) {
        paste0('  pivot1 dim\t= ') %>% cat
        dim(pivot1) %>% print
    }
    
    # Reshape data 2: Min_P
    pivot2 = gdata2 %>%
        group_by(SNPS,MAPPED_GENE) %>%
        summarize(Min_P = min(P.VALUE))
    if(debug) {
        paste0('  pivot2 dim\t= ') %>% cat
        dim(pivot2) %>% print
    }
    
    # Count P.VALUE number in a row
    pivot1_val = pivot1[,c(-1,-2)]
    ncols = ncol(pivot1_val)
    P_num = apply(pivot1_val,1,function(row) {
        ncols - sum(is.na(row))
    })
    if(debug) {
        paste0('  P_num length\t= ') %>% cat
        length(P_num) %>% print
    }
    
    # Save 1: Pivot table
    pivot = data.frame(
        pivot2,
        P_num,
        pivot1_val
    )
    f_name = paste0(out,'/',file_nm,'.tsv')
    write.table(pivot,f_name,sep='\t',quote=F,row.names=F)
    paste0('Write pivot:\t',f_name,'\n') %>% cat
    
    # Reshape data 3: SNP-Gene-Min_P
    pivot3_li = apply(pivot2,1,function(row) {
        geness  = strsplit(row[2]," - ") %>% unlist %>% as.character
        genes   = strsplit(geness,", ") %>% unlist %>% unique
        genes_n = length(genes)
        out = data.frame(
            SNPS  = rep(row[1],genes_n),
            GENES = genes,
            Min_P = rep(row[3],genes_n)
        )
        return(out)
    })
    pivot3 = data.table::rbindlist(pivot3_li)
    if(debug) {
        paste0('  pivot3 dim\t= ') %>% cat
        dim(pivot3) %>% print
    }
    
    # Save 2: Pivot3 table
    f_name2 = paste0(out,'/',file_nm,'_genes.tsv')
    write.table(pivot3,f_name2,sep='\t',quote=F,row.names=F)
    paste0('Write pivot3:\t',f_name2,'\n') %>% cat
}