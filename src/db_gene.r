help_message = '
Version: v2020-03-25
This is a function call for gene analysis to compile the Hi-C and eQTL data.

Usage: Rscript postgwas-exe.r --dbgene <function> --base <base files> --out <out folder>

Functions:
    hic_pair   Extract Hi-C linked SNP-Gene pairs.
    gtex_pair  Extract eQTL linked SNP-Gene pairs.
    summary    Summarizing GWAS SNPs-Gene pairs.

Global arguments:
    --base     <base files>
               For "hic_pair" function, two input files are:
                 [1] gwas_dist file, [2] gene_dist file
    --out      <out folder>
               Out folder path is mendatory. Default is "data" folder.

Required arguments:
    --bed      <default:FALSE>
               An optional argument for the "hic_pair" function.
               To save as BED format file.
    --nearest  An optional argument for the "summary" function.
               Add nearest gene summary table to the hic gene summary.
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
summary = function(
    f_paths = NULL,   # Input snp-gene pair files
    out     = 'data', # Out folder path
    nearest = NULL,   # Nearest gene summary table
    debug
) {
    # Load function-specific library
    suppressMessages(library(biomaRt))

    # Preparing...
    paste0('\n** Run function: db_filter.r/summary... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),''); 'ready\n' %>% cat

    # If the base path is folder, get the file list
    dir_name = basename(f_paths)
    paths1 = list.files(f_paths,full.name=T)
    n = length(paths1)
    paste0(n,' Files/folders input.\n') %>% cat

    paths = NULL
    for(i in 1:n) {
        paths2 = list.files(paths1[i],full.name=T)
        if(length(paths2)==0) {
            paths = c(paths,paths1[i])
            paste0('  ',i,' ',paths1[i],'\n') %>% cat
        } else {
            paths = c(paths,paths2)
            paste0('  ',i,' ',length(paths2),' files in the ',
                paths1[i] %>% basename,'\n') %>% cat
        }
    }
    n = length(paths)
    paste0('Total ',n,' files are input.\n') %>% cat

    # Read pair files
    enst_li = list(); snp_enst_li = list(); hic_ann_li = list()
    ensg_li = list(); snp_ensg_li = list(); gtex_ann_li = list()
    j = 0; k = 0
    for(i in 1:n) {
        tb = try(read.delim(paths[i],stringsAsFactors=F))
        if('try-error' %in% class(tb)) {
            paste0('  [ERROR] ',f_paths[i],'\n') %>% cat
            return()
        }
        tb_ncol = ncol(tb)
        # grep ENST
        enst = grep("ENST",tb[1,]$ensgid)
        ensg = grep("ENSG",tb[1,]$ensgid)
        if(length(enst)>0) {
            j = j+1
            enst_li[[j]] = tb$ensgid
            snp_enst_li[[j]] = tb[,1:2]
            hic_ann_li[[j]] = tb[,3:tb_ncol]
        } else if(length(ensg)>0) {
            k = k+1
            ensg_li[[k]] = tb$ensgid
            snp_ensg_li[[k]] = tb[,1:2]
            gtex_ann_li[[k]] = tb[,3:tb_ncol]
        } else paste0('[ERROR] Cannot find ENST/ENSG IDs.')
    }

    # Extract Enstids
    paste0('\n  Extracting Enstids... ') %>% cat
    enstids = enst_li %>% unlist %>% unique
    paste0(length(enstids),'.. ') %>% cat

    # Search biomaRt for enstids
    paste0('biomaRt... ') %>% cat
    ensembl   = useMart('ensembl',dataset='hsapiens_gene_ensembl')
    gene_attr1 = c('ensembl_transcript_id','ensembl_gene_id','hgnc_symbol')#,'description')
    gene_enst  = getBM(
        attributes = gene_attr1,
        filters = "ensembl_transcript_id",
        values = enstids,
        mart = ensembl
    ) %>% unique
    colnames(gene_enst)[1:2] = c('enstid','ensgid')
    dim(gene_enst) %>% print

    # Extract Ensgids
    paste0('  Extracting Ensgids... ') %>% cat
    ensgids   = ensg_li %>% unlist %>% unique
    paste0(length(ensgids),'.. ') %>% cat

    # Search biomaRt for ensgids
    paste0('biomaRt... ') %>% cat
    gene_attr2 = c('ensembl_gene_id','hgnc_symbol')#,'description')
    gene_ensg = getBM(
        attributes = gene_attr2,
        filters = "ensembl_gene_id",
        values = ensgids,
        mart = ensembl
    ) %>% unique
    colnames(gene_ensg)[1] = c('ensgid')
    dim(gene_ensg) %>% print

    # Merging enst-ensg pair
    paste0('\n  Merging enstid-ensgid pairs... ') %>% cat
    enstids_df = data.frame(enstid=enstids)
    enst_pair  = merge(enstids_df,gene_enst,by='enstid',all.x=T)
    
    ensgids_df = data.frame(ensgid=ensgids)
    ensg_pair_ = merge(ensgids_df,gene_ensg,by='ensgid',all.x=T)
    enstid_    = rep(NA,nrow(ensg_pair_))
    ensg_pair  = data.frame(
        enstid = enstid_,
        ensg_pair_
    )
    gene_ens = rbind(enst_pair,ensg_pair) %>% unique
    dim(gene_ens) %>% print
    paste0('  enstid na sum = ',sum(is.na(gene_ens$enstid)),'\n') %>% cat
    paste0('  ensgid na sum = ',sum(is.na(gene_ens$ensgid)),'\n') %>% cat

    # Parsing biomaRt gene name
    #paste0('\n  Parsing gene name... ') %>% cat
    #gene_name = lapply(gene_ens$description,function(x)
    #    strsplit(x,"\\ \\[")[[1]][1]) %>% unlist
    #gene_ens$description = gene_name
    #paste0(length(gene_ens$ensgid%>% unique),'.. ') %>% cat
    #dim(gene_ens) %>% print

    # Merge snp-enst and snp-ensg table
    paste0('\n  Merging SNPs to Enst/Ensg IDs... ') %>% cat
    snp_enst_df  = data.table::rbindlist(snp_enst_li) %>% unique
    colnames(snp_enst_df) = c('rsid_hic','enstid')
    snp_enst_mrg = merge(gene_ens,snp_enst_df,by='enstid',all.x=T)

    snp_ensg_df  = data.table::rbindlist(snp_ensg_li) %>% unique
    colnames(snp_ensg_df) = c('rsid_eqtl','ensgid')
    snp_mrg      = merge(snp_enst_mrg,snp_ensg_df,by='ensgid',all.x=T) %>% unique
    dim(snp_mrg) %>% print
    
    # out 1. Save as CSV file
    f_name1 = paste0(out,'/',dir_name,'_pairs.csv')
    write.csv(snp_mrg,f_name1,row.names=F)
    paste0('  Write file: ',f_name1,'\n') %>% cat
    return()
    
    # Merge hic and gtex annotations
    paste0('\n  Merging gene-hic/gtex annotations... ')
    enstids_df = data.frame(enstid=enstids)
    enst_pair  = merge(enstids_df,gene_ens,by='enstid',all.x=T)
    ensgids    = union(ensgids_,enst_pair$ensgid) %>% unique
    ensg_pair  = merge()

    # out 2. Save as TSV file
    return()
}

gtex_pair = function(
    f_paths = NULL,
    out     = 'data',
    debug
) {
    # Load function-specific library
    suppressMessages(library(tidyr))

    # Prepare...
    paste0('\n** Run function: db_gene.r/gtex_pair... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),''); 'ready\n' %>% cat

    # Read GTEx data and reshape
    paste0('  GTEx dim\t= ') %>% cat
    gtex = read.delim(f_paths,stringsAsFactors=F)
    dim(gtex) %>% print

    # Extract Ensgid
    paste0('  Extracting Ensgid.. ') %>% cat
    ensgid = lapply(gtex$gene_id,function(x) strsplit(x,"\\.")[[1]][1]) %>% unlist
    gtex$ensgid = ensgid
    paste0('done\n') %>% cat

    # Spread GTEx pair by tissue
    paste0('  Spreading GTEx data.. ') %>% cat
    gtex_reshape = gtex %>% spread(tissue,value=slope)
    n = ncol(gtex_reshape)
    gtex_pair = gtex_reshape[,c(7:n)]
    dim(gtex_pair) %>% print

    # Save as TSV file
    gene_n = unique(gtex_pair$ensgid) %>% length
    f_name = paste0(out,'/gtex_eqtl_',gene_n,'.tsv')
    write.table(gtex_pair,f_name,sep='\t',row.names=F,quote=F)
    paste0('  Write file: ',f_name,'\n') %>% cat
}

hic_pair = function(
    f_paths = NULL,   # Input [1] gwas_dist file, [2] gene_dist file
    out     = 'data', # Out folder path
    bed     = FALSE,   # Return as BED format
    debug
) {
    # Load function-specific library
    suppressMessages(library(biomaRt))
    if(bed=='TRUE') { bed = TRUE
    } else if(bed=='FALSE') bed = FALSE

    # Prepare...
    paste0('\n** Run function: db_gene.r/hic... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),''); 'ready\n' %>% cat

    # Read base files and filter by distance=0
    colnm = c('chr','start','end','name','hic_chr','hic_start','hic_end','hic_loop','dist')
    paste0('  ',f_paths[1],', length= ') %>% cat
    gwas = data.table::fread(f_paths[1],header=F)
    colnames(gwas) = colnm
    paste0(nrow(gwas)) %>% cat
    gwas0 = subset(gwas,dist==0) %>% unique
    gwas_0 = gwas0[,c(4,8)]
    paste0(',\toverlap= ',nrow(gwas_0),'\n') %>% cat

    paste0('  ',f_paths[2],', length= ') %>% cat
    gene = data.table::fread(f_paths[2],header=F)
    colnames(gene) = colnm
    paste0(nrow(gene)) %>% cat
    gene0 = subset(gene,dist==0) %>% unique
    gene_0 = gene0[,c(4,8)]
    paste0(',\toverlap= ',nrow(gene_0),'\n') %>% cat

    # Merge the gwas and gene data
    paste0('\n  Process gwas_loop.. ') %>% cat
    gwas_loop = gwas_0$hic_loop
    gwas_loop_id = lapply(gwas_loop,function(x) strsplit(x,"\\.")[[1]][1]) %>% unlist
    paste0(nrow(gwas_loop_id),'.. done\n') %>% cat

    paste0('  Process gene_loop.. ') %>% cat
    gene_loop = gene_0$hic_loop
    gene_loop_id = lapply(gene_loop,function(x) strsplit(x,"\\.")[[1]][1]) %>% unlist
    paste0(nrow(gene_loop_id),'.. done\n') %>% cat

    paste0('  Process merge.. ') %>% cat
    n = nrow(gwas_0); paste0(n,'.. ') %>% cat
    pair_li = lapply(c(1:n),function(i) {
        gwas_i = gwas0[i,]
        which_row = which(gene_loop_id==gwas_loop_id[i])
        loop = gene_loop[which_row]
        loop_g = gene_0[which_row,]
        if(nrow(loop_g)>0) {
            loop_g_ = subset(loop_g,hic_loop!=gwas_0$hic_loop[i])
            if(nrow(loop_g_)>0) {
                if(bed) { # if option bed = TRUE,
                    out = data.frame(
                        chr   = gwas_i$chr,
                        start = gwas_i$start,
                        end   = gwas_i$end,
                        rsid  = gwas_i$name
                    )
                } else { # if option bed = FALSE,
                    m = nrow(loop_g_)
                    genes = lapply(loop_g_$name,function(x)
                        strsplit(x,"\\_")[[1]][1]) %>% unlist
                    tags  = lapply(loop_g_$name,function(x) {
                        tag_v = strsplit(x,"\\_")[[1]]
                        l = length(tag_v)
                        paste0(tag_v[2:l],collapse="_") %>% return
                    }) %>% unlist
                    out = data.frame(
                        rsid      = rep(gwas_0$name,m),
                        ensgid    = genes,
                        gene_tag  = tags,
                        gwas_loop = rep(gwas_0$hic_loop,m),
                        gene_loop = loop_g_$hic_loop
                    )
                }
                return(out)
            } else return(NULL)
        } else return(NULL)
    })
    pair_df = data.table::rbindlist(pair_li) %>% unique
    dim(pair_df) %>% print

    # Save as TSV file
    file_nm = tools::file_path_sans_ext(basename(f_paths[1]))
    if(bed) {
        snp_n  = unique(pair_df$rsid) %>% length
        f_name = paste0(out,'/snp_hic_',file_nm,'_',snp_n,'.bed')
        write.table(pair_df,f_name,sep='\t',row.names=F,col.names=F,quote=F)
    } else {
        gene_n = unique(pair_df$ensgid) %>% length
        f_name = paste0(out,'/hic_',file_nm,'_',gene_n,'.tsv')
        write.table(pair_df,f_name,sep='\t',row.names=F,quote=F)
    }
    paste0('  Write file: ',f_name,'\n') %>% cat
}

## __INIT__ ##
db_gene = function(
    args = NULL
) {
    # Get help
    if(length(args$help)>0) { help = args$help
    } else                    help = FALSE
    if(help) {
        cat(help_message); quit()
    }

    # Global arguments
    if(length(args$base)>0)     b_path   = args$base
    if(length(args$out)>0)      out      = args$out
    if(length(args$debug)>0) {  debug    = args$debug
    } else                      debug    = FALSE

    # Required arguments
    if(length(args$bed)>0) {    bed      = args$bed
    } else                      bed      = FALSE


    # Run function
    source('src/pdtime.r'); t0=Sys.time()
    if(args$dbgene == 'hic_pair') {
        hic_pair(b_path,out,bed,debug)
    } else if(args$dbgene =='gtex_pair') {
        gtex_pair(b_path,out,debug)
    } else if(args$dbgene == 'summary') {
        summary(b_path,out,debug)
    } else {
        paste0('[Error] There is no such function "',args$dbgene,'" in gene: ',
            paste0(args$dbgene,collapse=', '),'\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}