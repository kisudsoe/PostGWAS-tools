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
    snp         Generating pivot table for snps.
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
    return(pivot)
}

snp_pivot = function(
    gdata   = NULL,   # GWAS Catalog data
    out     = 'data', # Output folder path
    file_nm = NULL,   # Output file name
    debug   = F
) {
    # Shrink data
    gdata2 = data.frame(
        SNPS         = gdata$SNPS,
        MAPPED_GENE  = gdata$MAPPED_GENE,
        P.VALUE      = gdata$P.VALUE,
        PMID         = gdata$PUBMEDID,
        FIRST.AUTHOR = gdata$FIRST.AUTHOR,
        DATE         = gdata$DATE
    ) %>% unique
    if(debug) {
        paste0('  gdata2 dim = ') %>% cat
        dim(gdata2) %>% print
    }
    
    # Reshape data 2: Min_P
    gdata2$MAPPED_GENE[gdata2$MAPPED_GENE==""] = NA # debug 2020-08-10
    pivot1 = gdata2 %>%
        group_by(SNPS,MAPPED_GENE,PMID,FIRST.AUTHOR,DATE) %>%
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
            SNPS          = row[1],
            GENES         = genes,
            PMID          = row[3],
            FIRST.AUTHOR  = row[4],
            DATE          = row[5],
            Min_P         = row[6]
        )
        return(out)
    })
    pivot = data.table::rbindlist(pivot_li)
    if(debug) {
        paste0('  pivot dim = ') %>% cat
        dim(pivot) %>% print
    }
    
    return(pivot)
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

    return(pivot)
}


generate_pivot = function(
    gwas = NULL,
    out = NULL,
    pivot = NULL,
    debug = FALSE
) {
    paste0('\n** Run function generate_pivot in gwas_catalog.r:\n') %>% cat
    f_gwas = list.files(gwas, full.names=T)
    f_n = length(f_gwas)
    if(f_n==0) { 
        f_gwas = gwas
        paste0('Read ', f_gwas,' = ') %>% cat
    } else {
        paste0(f_n,' files input -> [') %>% cat
    }
    f_n = length(f_gwas) # debug 6/10/2021
    gwas_li=list()
    for(i in 1:f_n) {
        paste0('.') %>% cat
        # Read GWAS Catalog file
        gwas_li[[i]] = read.delim(f_gwas[i],stringsAsFactors=F)
        file_nm = tools::file_path_sans_ext(f_gwas[i] %>% basename)
        if(debug=="TRUE"|debug=="T") debug=TRUE
        if(debug) {
            paste0('Read file, ',basename(f_gwas[i]),' = ') %>% cat
            dim(gdata) %>% print
        }
    }
    if(f_n>1) paste0('] = ') %>% cat
    gwas_merged  = data.table::rbindlist(gwas_li)
    dim(gwas_merged) %>% print

    # Generate TRAITS pivot table
    dir_nm = tools::file_path_sans_ext(gwas %>% basename)
    f_name = paste0(out,'/',dir_nm,c('_trait','_snp','_study'),'.tsv')
    if('trait' %in% pivot) {
        #paste0('\n** Run trait_pivot:\n') %>% cat
        trait_pivot = trait_pivot(gwas_merged, out, file_nm, debug)
        write.table(trait_pivot, f_name[1],sep='\t',quote=F)
        paste0('Write study pivot: ',f_name[1],'\n') %>% cat
    }
    # Generate SNP-Gene-Min_P table
    if('snp' %in% pivot) {
        #paste0('\n** Run snp_pivot:\n') %>% cat
        snps_pivot = snp_pivot(gwas_merged, out, file_nm, debug)
        write.table(snps_pivot,  f_name[2],sep='\t',quote=F,row.names=F)
        paste0('Write snp pivot: ',f_name[2],'\n') %>% cat
    }
    # Generate Study summary table
    if('study' %in% pivot) {
        #paste0('\n** Run study_pivot:\n') %>% cat
        study_pivot = study_pivot(gwas_merged, out, file_nm, debug)
        write.table(study_pivot,f_name[3],sep='\t',quote=F,row.names=F)
        paste0('Write traits pivot: ',f_name[3],'\n') %>% cat
    }
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
    
    # Generate out folder
    source('src/pdtime.r'); t0=Sys.time()
    ifelse(!dir.exists(out), dir.create(out),'')
    #if(args$gwas %in% c('trait','gene','study')) {
    if(all(args$gwas %in% c('trait','snp','study'))) {
        # Generate pivot tables
        generate_pivot(gwas, out, pivot, debug)
    } else if(args$gwas[1] == 'filter') {
        # Generate filtered SNP table
        paste0('\n** Run function gwas_filt:\n') %>% cat
        gdata = read.delim(gwas,stringsAsFactors=F)
        gwas_filt(gdata, out, p_criteria, debug)
    }
    paste0(pdtime(t0,1),'\n\n') %>% cat
}