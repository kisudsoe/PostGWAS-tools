help_message = '
db_filter, v2020-02-28
This is a function call for filtering data.

Usage: Rscript postgwas-exe.r --dbfilt <function> --base <base file(s)> --out <out folder> <...>

Functions:
    roadmap   Filtering Roadmap data by enhancer tags.
    gtex      Filtering GTEx data by eQTL p-value.
    gtex_ovl  Overlapping the GTEx data with the input GWAS SNPs.
    dist      Filtering distance data from Bedtools closest function.
    regulome  Filtering and overlapping by Regulome score ≥2b.
    lnc_ovl   Overlapping the lncRNASNP2 data with the input GWAS SNPs.

Global arguments:
    --base    <base file/folder>
              Base file/folder path is mendatory.
    --out     <out folder>
              Out folder path is mendatory. Default is "db" folder.

Required arguments:
    --ctype   <cell type id>
              An optional argument for the "roadmap" function.
              See cell-type number information at
              https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types.
    --enh     <default: TRUE>
              An optional argument for the "roadmap" function to filter enhancer regions.
    --sep     <default: FALSE>
              An optional argument for the "roadmap" function to generate cell-type seperated results.
    --pval    <p-value threshold>
              A required argument for the "gtex" function to filter significant eQTLs.
    --gtex    <Filtered GTEx RDS file path>
              A required argument for the "gtex_ovl" function to overlap GWAS SNPs with GTEx data.
    --tissue  <GTEx tissue name>
              An optional argument for the "gtex_ovl" function to filter a specific tissue.
    --regulm  <Regulome data folder>
              A required argument for the "regulome" function to load the data.
    --lncrna  <lncRNASNP2 data folder>
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
lncrna_overlap = function(
    snp_path = NULL,  # Input GWAS SNP file path
    lnc_path = NULL,  # lncRNASNP2 data downloaded folder path
    out      = 'data' # Out folder path
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/lncrna_overlap...\n') %>% cat
    snp = read.delim(snp_path,header=F)
    colnames(snp) = c('chr','start','end','rsid')
    dbsnp  = gsub('(.*)_.*','\\1',snp[,4])
    snp_df = cbind(snp,dbsnp)
    paste0('Input GWAS SNPs N = ',length(dbsnp),'\n') %>% cat

    # Read the lncRNASNP2 data files
    slnc_path = paste0(lnc_path,'/lncRNASNP2_snplist.txt.rds')
    ann_path  = paste0(lnc_path,'/lncrnas.txt.rds')
    dis_path  = paste0(lnc_path,'/lncrna-diseases_experiment.txt.rds')
    
    paste0('3 lncRNASNP2 data load...\n') %>% cat
    paste0('  Read: ',slnc_path) %>% cat
    snplnc = readRDS(slnc_path)
    paste0(';\t\tdim = ') %>% cat; dim(snplnc) %>% print

    paste0('  Read: ',ann_path) %>% cat
    ann = readRDS(ann_path)
    colnames(ann)[1] = 'lncRNA'
    paste0(';\t\t\t\tdim = ') %>% cat; dim(ann) %>% print

    paste0('  Read: ',dis_path) %>% cat
    dis = readRDS(dis_path)
    paste0(';\tdim = ') %>% cat; dim(dis) %>% print
    
    # Overlapping lncRNA data with input GWAS SNPs
    snp_lnc = merge(snp_df,snplnc,by='dbsnp')
    paste0('\nSummary =\n') %>% cat
    data.frame(
        lncRNA = unique(snp_lnc$lncRNA) %>% length,
        SNPs   = unique(snp_lnc$dbsnp) %>% length
    ) %>% print

    # Save as a BED file
    f_name1 = paste0(out,'/snp_lncrnasnp_',unique(snp_lnc$dbsnp)%>%length,'.bed')
    write.table(snp_lnc[,2:5],f_name1,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('\n  Write file: ',f_name1,'\n') %>% cat

    # Annotating GWAS SNPs with lncRNAs
    snp_lnc_ann = merge(snp_lnc[,c(1,6)],ann,by='lncRNA') %>% unique
    snp_lnc_ann_dis = merge(snp_lnc_ann,dis,by='lncRNA',all.x=T) %>% unique

    # Save as a TSV file
    f_name2 = paste0(out,'/lncrnasnp_',unique(snp_lnc_ann_dis$dbsnp)%>%length,'.tsv')
    write.table(snp_lnc_ann_dis,f_name2,row.names=F,quote=F,sep='\t')
    paste0('  Write file: ',f_name2,'\n') %>% cat
}

regulome_filt = function(
    snp_path = NULL,  # Input GWAS SNP file path
    reg_path = NULL,  # Regulome data downloaded folder path
    out      = 'data' # Out folder path
) {
    # Read a GWAS SNP BED file input
    paste0('\n** Run function: db_filter.r/regulome_filt...\n') %>% cat
    snp   = read.delim(snp_path,header=F)
    rsids = gsub('(.*)_.*','\\1',snp[,4])
    paste0('Input GWAS SNPs N\t= ') %>% cat; length(rsids) %>% print

    # Read the Regulome data files
    f_reg = list.files(reg_path)
    paths = paste0(reg_path,'/',f_reg)
    n = length(paths)
    paste0(n,' Regulome data load...\n') %>% cat
    reg_li = lapply(paths,function(path) {
        paste0('  Read: ',path,'; ') %>% cat
        f_ext = tools::file_ext(path)
        if(f_ext=='gz') out =  try(read.delim(gzfile(path),header=F))
        else if(f_ext=='rds') out = try(readRDS(path))
        else out = try(read.delim(path,header=F))
        colnames(out) = c('chr','id','rsid','description','level')
        paste0('dim = ') %>% cat; dim(out) %>% print
        return(out)
    })
    reg_df = data.table::rbindlist(reg_li)
    
    # Filter by score ≥2b
    reg_1f_only = subset(reg_df,level == '1f') %>% unique
    reg_2b = subset(reg_df,level %in% c("1a","1b","1c","1d","1e","1f","2a","2b")) %>% unique
    paste0('\n  Regulome score >=2b, SNPs\t\t= ') %>% cat
    nrow(reg_2b) %>% print
    paste0('  Functional motifs (1a~2b - 1f only)\t= ') %>% cat
    nrow(reg_2b)-nrow(reg_1f_only) %>% print

    # GWAS SNP counts
    snp_1f_only = subset(reg_1f_only,rsid %in% rsids) %>% unique
    snp_2b      = subset(reg_2b,rsid %in% rsids) %>% unique
    paste0('\n  Regulome >=2b, GWAS SNPs\t\t= ') %>% cat
    nrow(snp_2b) %>% print
    paste0('  GWAS SNPs occupied in
    functional motifs (1a~2b - 1f only)\t= ') %>% cat
    nrow(snp_2b) - nrow(snp_1f_only) %>% print

    # Save as TSV file
    f_name1 = paste0(out,'/regulome_',nrow(snp_2b),'.tsv')
    write.table(snp_2b,f_name1,row.names=F,quote=F,sep='\t')
    paste0('\nWrite file: ',f_name1,'\n') %>% cat

    # Save as BED file
    snp_rsid = data.frame(snp,rsid=rsids)
    snp_bed  = subset(snp_rsid,rsid %in% snp_2b$rsid) %>% unique
    f_name2  = paste0(out,'/snp_regulome2b_',nrow(snp_bed),'.bed')
    write.table(snp_bed[,1:4],f_name2,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('Write file: ',f_name2,'\n\n') %>% cat
}

distance_filt = function(
    f_path = NULL,   # Bedtools closest result paths
    out    = 'data', # Out folder path
    debug
) {
    # Preparing..
    f_name = tools::file_path_sans_ext(basename(f_path))
    paste0('File ',f_name,'... ') %>% cat
    rd = data.table::fread(f_path) %>% as.data.frame
    paste0('nrow= ',nrow(rd),'.. ') %>% cat

    # Extract data
    pos = apply(rd,1,function(row) paste0(row[5],':',row[6],'-',row[7])) %>% unlist
    rd_ = cbind(rd[,1:4],pos,rd[,c(8,ncol(rd))])
    colnames(rd_) = c('chr','start','end','rsid','ann_pos','tag','dist')
    dis = rd_$dist %>% as.character %>% as.numeric
    rd_df = cbind(rd_[,1:6],dis) %>% unique
    paste0('done\n') %>% cat

    # Filtering by distance 0 from the annotations
    rd_df2 = subset(rd_df,dis==0) %>% unique
    paste0('  Annotations occupied by SNPs\t= ') %>% cat; unique(rd_df2$ann_pos) %>% length %>% print
    paste0('  SNPs in annotations\t\t= ') %>% cat; unique(rd_df2$rsid) %>% length %>% print
    
    # Save as a TSV file
    #f_name1 = paste0(out,'/',f_name,'_filt.tsv')
    #write.table(rd_df2,f_name1,row.names=F,quote=F,sep='\t')
    #paste0('  Write file: ',f_name1,'\n') %>% cat

    # Save as a BED file
    snp_bed = rd_df2[,1:4] %>% unique
    snp_n   = unique(snp_bed$rsid) %>% length
    f_name2 = paste0(out,'/snp_',f_name,'_',snp_n,'.bed')
    write.table(snp_bed,f_name2,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('  Write file: ',f_name2,'\n') %>% cat
}

distance_filt_multi = function(
    f_paths = NULL,
    out     = 'data',
    debug
) {
    paste0('\n** Run function: db_filter.r/distance_filt...\n') %>% cat
    # If the base path is folder, get the file list
    if(length(f_paths)==1) {
        paths = list.files(f_paths,full.name=T)
        if(length(paths)==0) paths = f_paths
    } else paths = f_paths

    # Run function by each file path
    n = length(paths)
    o=lapply(c(1:n),function(i) {
        source('src/pdtime.r'); t0=Sys.time()
        f_path = paths[i]
        if(i%%10==0) paste0('  ',i,'/',n,' being processed.\n') %>% cat
        if(n<10) paste0('  ',i,'/',n,' ',path,'\n') %>% cat

        # Run function
        distance_filt(f_path,out,debug)
        paste0(pdtime(t0,2),'\n\n') %>% cat
    })
}

gtex_overlap = function(
    f_path = NULL,   # Inpust GWAS SNPs file path
    f_gtex    = NULL,   # Filtered GTEx RDS file path
    out       = 'data', # Out folder path
    tissue_nm = NULL,   # Optional tissue name
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/gtex_overlap...\n') %>% cat
    snp = read.delim(f_path,header=F)
    colnames(snp) = c('chr','start','end','ann')
    rsids = gsub('(.*)_.*','\\1',snp[,4])
    paste0('Input GWAS SNPs N\t= ',length(rsids),'\n') %>% cat

    # Load filtered GTEx RDS file
    f_ext = tools::file_ext(f_gtex)
    if(f_ext=='gz') {
        gtex = read.delim(gzfile(f_gtex),header=T)
    } else if(f_ext=='rds') {
        gtex = readRDS(f_gtex)
    } else {
        paste0('[ERROR] Unknown GTEx file format. gz/rds file is needed.')
        quit
    }
    colnames(gtex)[9] = 'rsid'
    paste0('  ',basename(f_gtex),', dim\t= ') %>% cat
    dim(gtex) %>% print

    # Overlapping eQTLs
    eqtls = subset(gtex, rsid %in% rsids)
    paste0('  Overlapped eQTL-gene pairs\t= ') %>% cat
    nrow(eqtls) %>% print

    # Filter by tissue (optional)
    if(!is.null(tissue_nm)) {
        paste0('\n[Option] ',tissue_nm) %>% cat
        eqtls = subset(eqtls,tissue == tissue_nm)
        f_name1 = paste0(out,'/gtex_signif_',tissue_nm,'_',unique(eqtls$rsid)%>%length,'.tsv')
        paste0(', dim = ') %>% cat; dim(eqtls) %>% print
    } else {
        f_name1 = paste0(out,'/gtex_signif_',unique(eqtls$rsid)%>%length,'.tsv')
    }
    paste0('  eQTLs N\t\t= ') %>% cat; unique(eqtls$rsid) %>% length %>% print
    paste0('  Associated eGenes\t= ') %>% cat; unique(eqtls$gene_id) %>% length %>% print

    # Save as a TSV file
    write.table(eqtls,f_name1,row.names=F,quote=F,sep='\t')
    paste0('\nWrite file: ',f_name1,'\n') %>% cat

    # Generate as BED file
    snp_rsid = data.frame(snp,rsid=rsids)
    snp_bed  = subset(snp_rsid,rsid %in% eqtls$rsid)[,1:4] %>% unique
    paste0('\n  GTEx eQTL BED, dim\t= ') %>% cat
    dim(snp_bed) %>% print
    paste0('  eQTL SNP N\t\t= ') %>% cat
    unique(snp_bed$ann) %>% length %>% print

    # Filter by tissue (optional)
    if(!is.null(tissue_nm)) {
        f_name2 = paste0(out,'/snp_gtex_',tissue_nm,'_',unique(snp_bed$ann)%>%length,'.bed')
    } else f_name2 = paste0(out,'/snp_gtex_',unique(snp_bed$ann)%>%length,'.bed')

    # Save as a BED file
    write.table(snp_bed,f_name2,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('\nWrite file: ',f_name2,'\n') %>% cat
}

gtex_filt = function(
    f_path = NULL, # RDS file path
    out    = 'db', # Out folder path
    pval   = NULL, # P-value criteria
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/gtex_filt...\n') %>% cat
    pval = pval %>% as.character %>% as.numeric
    paste0('  P-value threshold\t= ') %>% cat; print(pval)

    # Load GTEx BED files
    gtex_df = readRDS(f_path)
    paste0('  GTEx data, dim\t= ') %>% cat; dim(gtex_df) %>% print

    # Filtering GTEx data
    gtex_sig = subset(gtex_df,pval_nominal<pval)
    data.frame(gtex_sig$pval_nominal) %>% summary %>% print
    paste0('  GTEx <',pval,', dim\t= ') %>% cat; dim(gtex_sig) %>% print

    # Loading and merging annotation file
    f_name = paste0(out,'/gtex_signif_',pval,'.rds')
    saveRDS(gtex_sig,f_name)
    paste0('\nWrite file: ',f_name,'\n') %>% cat
}

roadmap_filt = function(
    f_path = NULL,  # Download folder path
    out    = 'db',  # Out folder path
    ctype  = NULL,  # Optional: Cell type ID
    enh    = TRUE,  # Optional: Filtering enhancer (Default: TRUE)
    sep    = FALSE, # Optional: Cell type separated results (Defualt: FALSE)
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/roadmap_filt...\n') %>% cat
    ifelse(!dir.exists(out), dir.create(out),
        paste0('  Directory already exists: ',out,'\n') %>% cat)

    if(length(ctype)>0) {
        cid    = formatC(ctype,width=3,flag='0') %>% as.character # Convert number to '###' format
        f_name = paste0(out,'/roadmap_',ctype,'_enh.bed')
    } else {
        cid    = formatC(c(1:129),width=3,flag='0') %>% as.character # 001~129
        if(enh) { f_name = paste0(out,'/roadmap_enh.bed')
        } else    f_name = paste0(out,'/roadmap_total.bed')
    }

    # Read Roadmap files
    paste0('  Reading files..\n') %>% cat
    n=length(cid)
    road_li = lapply(c(1:n), function(i) {
        if(i%%10==0) paste0('    ',i,'/',n,' being processed.\n') %>% cat
        #path = paste0(f_path,'E',cid[i],'_25_imputed12marks_hg38lift_dense.bed.rds') # hg38 by liftover
        path = paste0(f_path,'/E',cid[i],'_25_imputed12marks_dense.bed.rds') # hg19 original
        road = try(readRDS(path))
        if('try-error' %in% class(road)) { # if file is not found.
            paste0('  ',path,' - file not found.\n') %>% cat
            road_enh = NULL
        } else if(enh) { # If file exist.
            road_enh = subset(road,
                name %in% c("13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc"))
        } else road_enh = road
        if(n<10) {
            paste0('  ',path,' ') %>% cat
            dim(road_enh) %>% print
        }
        if(sep) {
            road_enh = road_enh[,c(1:3,6)]
            write.table(road_enh,paste0(out,'/roadmap_',cid[i],'_enh.bed'),
                row.names=F,col.names=F,quote=F,sep='\t')
        } else {
            return(road_enh)
        }
    })

    if(sep) {
        paste0('  Finished processing ',n,' files.\n\n') %>% cat
    } else {
        paste0('  Finished reading and filtering ',n,' files.\n\n') %>% cat
        road_enh = data.table::rbindlist(road_li)
        rm(road_li)

        # Save file
        road_enh = road_enh[,c(1:3,6)]
        write.table(road_enh,f_name,row.names=F,col.names=F,quote=F,sep='\t')
        paste0('Write file: ',f_name,'\n') %>% cat
    }
}

db_filter = function(
    args = NULL
) {
    if(length(args$help)>0) {   help     = args$help
    } else                      help     = FALSE
    if(help) {                  cat(help_message); quit() }

    # Global arguments
    if(length(args$base)>0)     b_path   = args$base
    if(length(args$out)>0)      out      = args$out
    if(length(args$debug)>0) {  debug    = args$debug
    } else                      debug    = FALSE

    # Reguired arguments
    if(length(args$ctype)>0) {  ctype    = args$ctype
    } else                      ctype    = NULL
    if(length(args$enh)>0) {    enh      = args$enh
    } else                      enh      = TRUE
    if(length(args$sep)>0) {    sep      = args$sep
    } else                      sep      = FALSE
    if(length(args$pval)>0)     pval     = args$pval
    if(length(args$gtex)>0)     gtex     = args$gtex
    if(length(args$tissue)>0) { tissue   = args$tissue
    } else                      tissue   = NULL
    if(length(args$regulm)>0)   reg_path = args$regulm
    if(length(args$lncrna)>0)   lnc_path = args$lncrna

    # Run function
    source('src/pdtime.r'); t0=Sys.time()
    if(args$dbfilt == 'roadmap') {
        roadmap_filt(b_path,out,ctype,enh,sep,debug)
    } else if(args$dbfilt == 'gtex') {
        gtex_filt(b_path,out,pval,debug)
    } else if(args$dbfilt == 'gtex_ovl') {
        gtex_overlap(b_path,gtex,out,tissue,debug)
    } else if(args$dbfilt == 'dist') {
        paste0('Input file N\t= ') %>% cat; length(b_path) %>% print
        distance_filt_multi(b_path,out,debug)
    } else if(args$dbfilt == 'regulome') {
        regulome_filt(b_path,reg_path,out)
    } else if(args$dbfilt == 'lnc_ovl') {
        lncrna_overlap(b_path,lnc_path,out)
    } else {
        paste0('[Error] There is no such function "',args$dbfilt,'" in db_filter: ',
            paste0(args$ldlink,collapse=', '),'\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}