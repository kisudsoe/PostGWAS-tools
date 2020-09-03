help_message = '
gwas_catalog, v2020-08-11
This is a function for GWAS Catalog data.

Usage:
    Rscript postgwas-exe.r --gwas trait --base <base TSV file> --out <out folder>
    Rscript postgwas-exe.r --gwas gene --base <base TSV file> --out <out folder>
    Rscript postgwas-exe.r --gwas study --base <base TSV file> --out <out folder>
    Rscript postgwas-exe.r --gwas trait gene study --base <base TSV file> --out <out folder>
    Rscript postgwas-exe.r --gwas filter --base <base TSV file> --out <out folder> --p.criteria 5e-8


Functions:
    trait       Generating pivot table for traits.
    gene        Generating pivot table for genes.
    study       Generating summary table for studies.
    filter      Filtering SNPs by P-values.

Global arguments:
    --base      <EFO0001359.tsv>
                An input base TSV file downloaded from GWAS Catalog is mendatory.
    --out       <Default:data folder>
                Out files target directory path that is mendatory. Default is "data" folder.

Required arguments:
    --p.criteria <5e-8>
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
    paste0('  gwas dim = ') %>% cat
    dim(snp) %>% print
    
    # Filtering by p-value criteria
    p_criteria = p_criteria %>% as.character %>% as.numeric
    paste0('  gwas (',p_criteria,') = ') %>% cat
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
    paste0('Write gwas filter: ',f_name,'\n') %>% cat
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
        paste0('  AUTHORS length = ') %>% cat
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
        paste0('  studies dim = ') %>% cat
        dim(studies) %>% print
    }
    
    # Reshape data 1: No filter (P <1e-5)
    pivot_1e5 = studies %>%
        group_by(AUTHORS,DATE,SAMPLES,TRAITS) %>%
        tally(name='SNP_N') #%>% # A wrapper for summarise that will either call n() or sum(n)
        #summarize(SNP_n = count(SNPS)) #%>%
        #spread('TRAITS','SNP_n')
    if(debug) {
        paste0('  pivot (<1e-5) dim = ') %>% cat
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
        paste0('  pivot (<5e-8) dim = ') %>% cat
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
        paste0('  pivot merge dim = ') %>% cat
        dim(pivot_5e8) %>% print
    }
    f_name = paste0(out,'/',file_nm,'_studies.tsv')
    write.table(pivot,f_name,sep='\t',quote=F)
    paste0('Write study pivot: ',f_name,'\n') %>% cat
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
        paste0('  gdata2 dim = ') %>% cat
        dim(gdata2) %>% print
    }
    
    # Reshape data 2: Min_P
    gdata2$MAPPED_GENE[gdata2$MAPPED_GENE==""] = NA # debug 2020-08-10
    pivot1 = gdata2 %>%
        group_by(SNPS,MAPPED_GENE) %>%
        summarize(Min_P = min(P.VALUE))
    if(debug) {
        paste0('  pivot1 dim = ') %>% cat
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
        paste0('  pivot dim = ') %>% cat
        dim(pivot) %>% print
    }
    
    # Save: Pivot table
    f_name = paste0(out,'/',file_nm,'_genes.tsv')
    write.table(pivot,f_name,sep='\t',quote=F,row.names=F)
    paste0('Write gene pivot: ',f_name,'\n') %>% cat
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
        paste0('  TRAITS length = ') %>% cat
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
        paste0('  gdata2 dim = ') %>% cat
        dim(gdata2) %>% print
    }
    
    # Reshape data 1: Spread by TRAITS
    pivot1 = gdata2 %>%
        group_by(SNPS,MAPPED_GENE,TRAITS) %>%
        summarize(Min_P = min(P.VALUE)) %>%
        spread('TRAITS','Min_P')
    if(debug) {
        paste0('  pivot1 dim = ') %>% cat
        dim(pivot1) %>% print
    }
    
    # Reshape data 2: Min_P
    pivot2 = gdata2 %>%
        group_by(SNPS,MAPPED_GENE) %>%
        summarize(Min_P = min(P.VALUE))
    if(debug) {
        paste0('  pivot2 dim = ') %>% cat
        dim(pivot2) %>% print
    }
    
    # Count P.VALUE number in a row
    pivot1_val = pivot1[,c(-1,-2)]
    ncols = ncol(pivot1_val)
    P_num = apply(pivot1_val,1,function(row) {
        ncols - sum(is.na(row))
    })
    if(debug) {
        paste0('  P_num length = ') %>% cat
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
    paste0('Write TRAITS pivot: ',f_name,'\n') %>% cat
}
## Functions End ##

## __INITE__ function
gwas_catalog = function(
    args = args
) {
    if(length(args$help)>0) {       help    = args$help
    } else                          help    = FALSE
    if(help) {                      cat(help_message); quit() }
    
    if(length(args$base)>0)         gwas  = args$base
    if(length(args$out)>0)          out   = args$out
    if(length(args$debug)>0) {      debug = args$debug
    } else                          debug = FALSE
        
    if(length(args$gwas)>0)         pivot = args$gwas
    if(length(args$p.criteria)>0) { p_criteria = args$p.criteria
    } else                          p_criteria = 5e-8
    
    # Read GWAS Catalog file
    gdata = read.delim(gwas,stringsAsFactors=F)
    file_nm = tools::file_path_sans_ext(gwas %>% basename)
    if(debug=="TRUE"|debug=="T") debug=TRUE
    if(debug) {
        paste0('Read file, ',basename(gwas),' = ') %>% cat
        dim(gdata) %>% print
    }

    # Generate out folder
    ifelse(!dir.exists(out), dir.create(out),'')

    source('src/pdtime.r'); t0=Sys.time()
    if('trait' %in% args$gwas) {
        # Generate TRAITS pivot table
        paste0('\n** Run function trait_pivot:\n') %>% cat
        trait_pivot(gdata, out, file_nm, debug)
    }
    if('gene' %in% args$gwas) {
        # Generate SNP-Gene-Min_P table
        paste0('\n** Run function gene_pivot:\n') %>% cat
        gene_pivot(gdata, out, file_nm, debug)
    }
    if('study' %in% args$gwas) {
        # Generate Study summary table
        paste0('\n** Run function study_pivot:\n') %>% cat
        study_pivot(gdata, out, file_nm, debug)
    } else if(args$gwas == 'filter') {
        # Generate filtered SNP table
        paste0('\n** Run function gwas_filt:\n') %>% cat
        gwas_filt(gdata, out, p_criteria, debug)
    }
    paste0(pdtime(t0,1),'\n\n') %>% cat
}