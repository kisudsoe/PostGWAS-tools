help_messages = function() {
    help_message = '
gwas_ldlink, v2020-01-21
This is a function for LDlink data.

Usage: Rscript postgwas.r --ldlink [options] <-b base file> <-o out folder> <...>
    --ldlink [Options: dn/fl]
        dn   This is a function for LDlink data download.
        fl   This is a function for LDlink data filter.

Required arguments:
    --base   <EFO0001359.tsv>
             One base TSV file is mendatory.
    --out    <data>
             Out folder path is mendatory. Default is "data" folder.
    --popul  <CEU TSI FIN GBR IBS ...>
             A download dependent argument. One or more population option have to be included.
    --r2d    <1/2/3/4>
             A filter dependent argument. Choose either number option.
             1) r2>0.6 and Dprime=1  <- The most stringent criteria.
             2) r2>0.6               <- Usual choice to define LD association.
             3) Dprime=1
             4) r2>0.6 or Dprime=1
'
    help_message %>% cat
    quit()
}

## Load libraries ##
suppressMessages(library(LDlinkR))
suppressMessages(library(dplyr))
suppressMessages(library(biomaRt))

## Functions Start ##
ldlink_fl = function(
    snp_path = NULL,   # GWAS file path
    ld_path  = NULL,   # LDlink download folder path
    out      = 'data', # out folder path
    r2d      = NULL,   # LDlink filter criteria. R2>0.6 and/or D'=1
    debug    = F
) {
    # Read downloaded files
    paste0('Read download files... ') %>% cat
    snpdf     = read.delim(snp_path) %>% unique
    snpids    = snpdf$rsid %>% unique
    col_names = c('No','RS_Number','Coord','Alleles','MAF',
        'Distance','Dprime','R2','Correlated_Alleles','RegulomeDB','Function')
    ldlink    = paste0(ld_path,'/',snpids,'.txt')
    snptb     = data.frame(snpids=snpids, ldlink=ldlink)
    
    paste0(nrow(snptb),'\n') %>% cat
    ldlink_li = apply(snptb,1,function(row) {
        tb1 = try(
            read.table(as.character(row[2]),sep='\t',
                header=F,skip=1,stringsAsFactors=F,
                col.names=col_names) )
        if('try-error' %in% class(tb1)) {
            paste0('  ',row[1],'\n') %>% cat
            return(NULL)
        } else { # If no errors occurred,
            tb2 = data.frame(SNPid=rep(row[1],nrow(tb1)),tb1)
            return(tb2)
        }
    })
    ldlink_df = data.table::rbindlist(ldlink_li) %>% unique
    paste0('  Read LDlink results\t= ') %>% cat; dim(ldlink_df) %>% print

    # Filter the LDlink data
    if(r2d==1) {
        cat('Filtering by "r2 > 0.6 and Dprime = 1":\n')
        ldlink_1 = subset(ldlink_df,R2>0.6 & Dprime==1)
    } else if(r2d==2) {
        cat('Filtering by "r2 > 0.6":\n')
	    ldlink_1 = subset(ldlink_df,R2>0.6) # r2 > 0.6
    } else if(r2d==3) {
        cat('Filtering by "Dprime = 1":\n')
	    ldlink_1 = subset(ldlink_df,Dprime==1) # D' = 1
    } else if(r2d==4) {
        cat('Filtering by "r2 > 0.6 or Dprime = 1":\n')
	    ldlink_1 = subset(ldlink_df,R2>0.6 | Dprime==1) # r2 > 0.6 or D' = 1
    } else cat('Which filtering option is not supported.\n')
    ldlink_2 = data.frame(
        gwasSNPs = ldlink_1$SNPid,
        ldSNPs   = ldlink_1$RS_Number
    ) %>% unique
    paste0('  Filtered data dimension \t= ') %>% cat; dim(ldlink_2) %>% print
    ldlink_ = ldlink_2[!ldlink_2$`ldSNPs` %in% c("."),] # Exclude no rsid elements
    ex = nrow(ldlink_2[ldlink_2$`ldSNPs` %in% c("."),])
    paste0('  Excluded no rsid elements\t= ') %>% cat; print(ex)
    
    # Basic summary of LDlink results
    paste0('Basic summary of LDlink results:\n') %>% cat
    paste0('  SNP Tier 1\t\t\t= ',length(snpids),'\n') %>% cat
    snp_t2 = setdiff(ldlink_$ldSNPs,snpids) %>% unique
    paste0('  SNP Tier 2\t\t\t= ',length(snp_t2),'\n') %>% cat
    snp_cand = union(ldlink_$ldSNPs,snpids) %>% unique
    paste0('  SNP candidates\t= ',length(snp_cand),'\n') %>% cat

    # Search biomart to get coordinates
    paste0('Search biomart for hg38 coordinates:\n') %>% cat
    paste0('  Query SNPs\t= ') %>% cat; length(snp_cand) %>% print
    hg38_snp = useMart(biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
    snp_attr = c("refsnp_id","chr_name","chrom_start","chrom_end")
    snps_bio = getBM(
        attributes = snp_attr,
        filters    = "snp_filter",
        values     = snp_cand,
        mart       = hg38_snp) %>% unique
    colnames(snps_bio) = c('rsid','chr','start','end')
    snps_bio_     = subset(snps_bio,chr %in% c(1:22,'X','Y'))
    snps_bio_[,2] = paste0('chr',snps_bio_[,2])
    snps_bio_[,3] = as.numeric(as.character(snps_bio_[,3]))-1
    paste0('  Result table\t= ') %>% cat; dim(snps_bio_) %>% print

    # Merge the biomart result with the GWAS SNP list 
    snps_ = data.frame(rsid=snp_cand)
    snps_merge = merge(snps_,snps_bio_,by='rsid',all.x=TRUE)
    snp_bed = data.frame(
        chr   = snps_merge$chr,
        start = snps_merge$start,
        end   = snps_merge$end,
        rsid  = snps_merge$rsid
    )
    paste0('  Merged table\t= ') %>% cat; dim(snp_bed) %>% print

    # Save as BED file
    f_name = paste0(out,'/gwas_biomart.bed')
    write.table(snp_bed,f_name,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('Write file:\t',f_name,'\n') %>% cat
}

ldlink_dn = function(
    snp_path = NULL,   # gwassnp_summ: 'filter' result file path
    out      = 'data', # out folder path
    popul    = NULL,   # population filter option for LDlink
    debug    = F
) {
    # Download from LDlink
    paste0('\n** Run function ldlink_dn... ') %>% cat
    snps = read.delim(snp_path)
    rsid = snps$rsid %>% unique
    paste0(rsid%>%length,'.. ') %>% cat
    token = '669e9dc0b428' # Seungsoo Kim's token
    LDproxy_batch(snp=rsid, pop=popul, r2d='d', token=token, append=F)
    paste0('done\n') %>% cat
    
    # Rename downloaded file
    f_name  = paste0(rsid,'.txt')
    f_name1 = paste0(out,'/',f_name)
    file.rename(f_name,f_name1)
    paste0('  Files are moved to target folder:\t',out,'\n') %>% cat
}

gwas_ldlink = function(
    args = args
) {
    if(length(args$help)>0)     help    = args$help
    if(help) help_messages()
    
    if(length(args$base)>0)     b_path  = args$base
    if(length(args$out)>0)      out     = args$out
    if(length(args$debug)>0) {  debug   = args$debug
    } else                      debug   = FALSE
    
    if(length(args$popul)>0)    popul   = args$popul
    if(length(args$r2d)>0)      r2d     = args$r2d
    
    if(args$ldlink == 'dn') {
        source('src/pdtime.r'); t0=Sys.time()
        ldlink_dn(b_path,out,popul,debug)
        paste0(pdtime(t0,1),'\n') %>% cat
    } else if(args$ldlink == 'fl') {
        source('src/pdtime.r'); t0=Sys.time()
        ldlink_fl(b_path[1],b_path[2],out,r2d,debug)
        paste0(pdtime(t0,1),'\n') %>% cat
    }
}