help_message = '
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
'

## Load libraries ##
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

## Functions Start ##
gwas_filt = function(
    gdata      = NULL,   # GWAS Catalog file path
    out        = 'data', # Output folder path
    p_criteria = 5e-8,   # P-value criteria
    debug      = F
) {
    # Read GWAS Catalog file
    snp = gdata
    paste0('  gwas dim\t\t= ') %>% cat
    dim(snp) %>% print
    
    # Filtering by p-value criteria
    p_criteria = p_criteria %>% as.character %>% as.numeric
    paste0('  gwas (',p_criteria,')\t\t= ') %>% cat
    snp_ = subset(snp,`P.VALUE`<p_criteria)
    dim(snp_) %>% print
    
    # Generate return table
    coord  = paste0('chr',snp_$`CHR_ID`,':',snp_$`CHR_POS`); head(coord)
    snp_df = data.frame(
        rsid       = snp_$`SNPS`,
        coord_hg38 = coord,pval=snp_$`P.VALUE`,
        cytoband   = snp_$`REGION`) %>% unique
    snp_n  = unique(snp_df$rsid) %>% length
    
    f_name = paste0(out,'/gwas_',p_criteria,'_',snp_n,'.tsv')
    write.table(snp_df,f_name,sep='\t',quote=F,row.names=F)
    paste0('Write gwas filter:\t',f_name,'\n') %>% cat
}

study_pivot = function(
    gdata   = NULL,   # GWAS Catalog file path
    out     = 'data', # Output folder path
    file_nm = NULL,   # Output file name
    debug   = F
) {
    # Generate AUTHORS
    AUTHORS = paste0(gdata$FIRST.AUTHOR,' (',gdata$PUBMEDID,')')
    if(debug) {
        paste0('  AUTHORS length\t= ') %>% cat
        print(length(AUTHORS))
    }
    
    # Generate SAMPLES
    SAMPLES = paste0(
        'INIT: ', gdata$INITIAL.SAMPLE.SIZE,
        " //REP: ", gdata$REPLICATION.SAMPLE.SIZE)
    
    # Shrink data: SAMPLE info
    studies = data.frame(
        AUTHORS = AUTHORS,
        DATE    = gdata$DATE,
        SAMPLES = SAMPLES,
        TRAITS  = gdata$DISEASE.TRAIT,
        SNPS    = gdata$SNPS,
        Pvalue  = gdata$P.VALUE
    ) %>% unique
    if(debug) {
        paste0('  studies dim\t\t= ') %>% cat
        dim(studies) %>% print
    }
    
    # Reshape data 1: No filter (P <1e-5)
    pivot_1e5 = studies %>%
        group_by(AUTHORS,DATE,SAMPLES,TRAITS) %>%
        tally(name='SNP_N') #%>% # A wrapper for summarise that will either call n() or sum(n)
        #summarize(SNP_n = count(SNPS)) #%>%
        #spread('TRAITS','SNP_n')
    if(debug) {
        paste0('  pivot (<1e-5) dim\t= ') %>% cat
        dim(pivot_1e5) %>% print
    }
    
    # Reshape data 2: Filtered by P <5e-8
    studies2  = subset(studies,Pvalue < 5e-8)
    pivot_5e8 = studies2 %>%
        group_by(AUTHORS,DATE,SAMPLES,TRAITS) %>%
        tally(name='SNP_N') #%>% # A wrapper for summarise that will either call n() or sum(n)
        #summarize(SNP_n = count(SNPS)) #%>%
        #spread('TRAITS','SNP_n')
    if(debug) {
        paste0('  pivot (<5e-8) dim\t= ') %>% cat
        dim(pivot_5e8) %>% print
    }
    
    # Merge and Save: Pivot table
    pivot = left_join(pivot_1e5, pivot_5e8, by=c(
        'AUTHORS' = 'AUTHORS',
        'DATE'    = 'DATE',
        'SAMPLES' = 'SAMPLES',
        'TRAITS'  = 'TRAITS'
    ))
    colnames(pivot)[5:6] = c('SNP_1e-5','SNP_5e-8')
    if(debug) {
        paste0('  pivot merge dim\t= ') %>% cat
        dim(pivot_5e8) %>% print
    }
    f_name = paste0(out,'/',file_nm,'_studies.tsv')
    write.table(pivot,f_name,sep='\t',quote=F)
    paste0('Write study pivot:\t',f_name,'\n') %>% cat
}

gene_pivot = function(
    gdata   = NULL,   # GWAS Catalog data
    out     = 'data', # Output folder path
    file_nm = NULL,   # Output file name
    debug   = F
) {
    # Shrink data
    gdata2 = data.frame(
        SNPS        = gdata$SNPS,
        MAPPED_GENE = gdata$MAPPED_GENE,
        P.VALUE     = gdata$P.VALUE
    ) %>% unique
    if(debug) {
        paste0('  gdata2 dim\t= ') %>% cat
        dim(gdata2) %>% print
    }
    
    # Reshape data 2: Min_P
    pivot1 = gdata2 %>%
        group_by(SNPS,MAPPED_GENE) %>%
        summarize(Min_P = min(P.VALUE))
    if(debug) {
        paste0('  pivot1 dim\t= ') %>% cat
        dim(pivot1) %>% print
    }
    
    # Reshape data 3: SNP-Gene-Min_P
    pivot_li = apply(pivot1,1,function(row) {
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
    pivot = data.table::rbindlist(pivot_li)
    if(debug) {
        paste0('  pivot dim\t= ') %>% cat
        dim(pivot) %>% print
    }
    
    # Save: Pivot table
    f_name = paste0(out,'/',file_nm,'_genes.tsv')
    write.table(pivot,f_name,sep='\t',quote=F,row.names=F)
    paste0('Write gene pivot:\t',f_name,'\n') %>% cat
}

trait_pivot = function(
    gdata   = NULL,   # GAWS Catalog data
    out     = 'data', # Output folder path
    file_nm = NULL,   # Output file name
    debug   = F
) {
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
    
    # Save: Pivot table
    pivot = data.frame(
        pivot2,
        P_num,
        pivot1_val
    )
    f_name = paste0(out,'/',file_nm,'_snps.tsv')
    write.table(pivot,f_name,sep='\t',quote=F,row.names=F)
    paste0('Write TRAITS pivot:\t',f_name,'\n') %>% cat
}
## Functions End ##

## __INITE__ function
gwas_summ = function(
    gwas       = NULL,   # GAWS Catalog file path
    out        = 'data', # Output folder path
    p_criteria = 5e-8,   # P-value criteria
    pivot      = c('trait','gene','study','filter'), # Generate pivot table options
    debug      = F
) {
    # Read GWAS Catalog file
    gdata = read.delim(gwas,stringsAsFactors=F)
    file_nm = tools::file_path_sans_ext(gwas %>% basename)
    if(debug) {
        paste0('Read file, ',basename(gwas),'\t= ') %>% cat
        dim(gdata) %>% print
    }

    # Generate TRAITS pivot table
    if('trait' %in% pivot) {
        paste0('\n** Run function trait_pivot:\n') %>% cat
        trait_pivot(gdata, out, file_nm, debug)
    }
    
    # Generate SNP-Gene-Min_P table
    if('gene'  %in% pivot) {
        paste0('\n** Run function gene_pivot:\n') %>% cat
        gene_pivot(gdata, out, file_nm, debug)
    }
    
    # Generate Study summary table
    if('study' %in% pivot) {
        paste0('\n** Run function study_pivot:\n') %>% cat
        study_pivot(gdata, out, file_nm, debug)
    }
    
    # Generate filtered SNP table
    if('filter' %in% pivot) {
        paste0('\n** Run function gwas_filt:\n') %>% cat
        gwas_filt(gdata, out, p_criteria, debug)
    }
}

gwas_catalog = function(
    args = args
) {
    if(length(args$help)>0) {       help    = args$help
    } else                          help    = FALSE
    if(help)                        cat(help_message)
    
    if(length(args$base)>0)         gwas  = args$base
    if(length(args$out)>0)          out   = args$out
    if(length(args$debug)>0) {      debug = args$debug
    } else                          debug = FALSE
        
    if(length(args$gwas)>0)         pivot = args$gwas
    if(length(args$p.criteria)>0) { p_criteria = args$p.criteria
    } else                          p_criteria = 5e-8
    
    gwas_summ(gwas,out,p_criteria,pivot,debug)
}