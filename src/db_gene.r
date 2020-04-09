help_message = '
Version: v2020-03-25
This is a function call for gene analysis to compile the Hi-C and eQTL data.

Usage: Rscript postgwas-exe.r --dbgene <function> --base <base files> --out <out folder>

Functions:
    hic_pair   Extract Hi-C linked SNP-Gene pairs.
    gtex_pair  Extract eQTL linked SNP-Gene pairs.
    summary    Summarizing GWAS SNPs-Gene pairs.
    david_go   Summarizing DAVID GO analysis result.

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
    --nearest  <nearest gene summary file path>
               An optional argument for the "summary" function.
               Add nearest gene summary table to the hic gene summary.
    --fdr      <default:0.05>
               An optional argument for the "david_go" function.
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
david_go = function(
    go_path = NULL,   # Input DAVID GO analysis result file path
    out     = 'data', # Out folder path
    fdr     = 0.05,   # GO FDR criteria
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_gene.r/david_go... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),''); 'ready\n' %>% cat

    # Read go result file
    paste0('  Read DAVID result file\t') %>% cat
    go = read.delim(go_path)
    dim(go) %>% print
    
    paste0('  FDR criteria\t= ',fdr,'\n') %>% cat
    paste0('  Filtered GO terms\t= ') %>% cat
    go_ = subset(go,FDR<fdr)
    dim(go_) %>% print
    
    paste0('  Extract the gene list... ') %>% cat
    genes_li = lapply(go_$Genes,function(x) { strsplit(x%>%as.character,"\\, ")[[1]] })
    names(genes_li) = go_$Term
    genes_n = unlist(genes_li) %>% unique %>% length
    paste0(genes_n,'.. done\n') %>% cat
    
    paste0('  Venn analysis... ') %>% cat
    source('src/venn_analysis.r')
    venn_tf = venn_analysis(genes_li,wfig=F,wfile=F)
    dim(venn_tf) %>% print
    
    f_name = paste0(out,'/go_analysis.csv')
    write.csv(venn_tf,f_name,quote=F,row.names=F)
    paste0('  Write file: ',f_name,'\n') %>% cat
}

summary = function(
    f_paths = NULL,   # Input snp-gene pair files
    out     = 'data', # Out folder path
    nearest = NULL,   # Nearest gene summary table
    debug
) {
    # Load function-specific library
    suppressMessages(library(biomaRt))
    suppressMessages(library(tools))

    # Preparing...
    paste0('\n** Run function: db_gene.r/summary... ') %>% cat
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
    if(!is.null(nearest)) paste0('Nearest file: ',nearest,'\n') %>% cat

    # Read pair files
    ensg_li = list(); snp_ensg_li = list(); tb_li = list()
    for(i in 1:n) {
        tb = try(read.delim(paths[i],stringsAsFactors=F))
        if('try-error' %in% class(tb)) {
            paste0('  [ERROR] ',f_paths[i],'\n') %>% cat
            return()
        }
        ensg_li[[i]] = tb$ensgid
        snp_ensg = paste0(tb$ensgid,'-',tb$rsid)
        tb$ensgid_rsid = snp_ensg
        snp_ensg_li[[i]] = snp_ensg
        tb_li[[i]] = tb
    }

    # Read nearest gene file
    if(!is.null(nearest)) {
        tb = read.csv(nearest,stringsAsFactors=F)
        colnames(tb)[2] = 'ensgid'
        tb_ncol = ncol(tb)
        n = length(ensg_li)
        ensg_li[[n+1]] = tb$ensgid
        snp_ensg = paste0(tb$ensgid,'-',tb$rsid)
        tb$ensgid_rsid = snp_ensg
        snp_ensg_li[[n+1]] = snp_ensg
        tb_li[[n+1]] = tb
    }

    # Extract Ensgids
    paste0('  Extracting Ensgids... ') %>% cat
    ensgids = ensg_li %>% unlist %>% unique
    paste0(length(ensgids),'.. ') %>% cat

    # Search biomaRt for ensgids
    paste0('biomaRt... ') %>% cat
    ensembl   = useMart('ensembl',dataset='hsapiens_gene_ensembl')
    gene_attr2 = c('ensembl_gene_id','hgnc_symbol','description')
    gene_ensg = getBM(
        attributes = gene_attr2,
        filters = "ensembl_gene_id",
        values = ensgids,
        mart = ensembl
    ) %>% unique
    colnames(gene_ensg)[1] = c('ensgid')
    dim(gene_ensg) %>% print

    # Parsing biomaRt gene name
    paste0('\n  Parsing gene name... ') %>% cat
    gene_name = lapply(gene_ensg$description,function(x)
        strsplit(x,"\\ \\[")[[1]][1]) %>% unlist
    gene_ensg$description = gene_name
    paste0(length(gene_ensg$ensgid%>% unique),'.. ') %>% cat
    dim(gene_ensg) %>% print

    # Merge gene annotation to snp_ensg list
    paste0('  ENSGid-Rsid list... ') %>% cat
    ensgid_rsid = snp_ensg_li %>% unlist %>% unique
    paste0(length(ensgid_rsid),'.. ') %>% cat
    ensgid_split_li = lapply(ensgid_rsid,function(x) {
        row = strsplit(x,"\\-")[[1]]
        data.frame(
            ensgid = row[1],
            rsid = row[2]
        )
    })
    ensgid_split_df = data.table::rbindlist(ensgid_split_li)
    ensgid_rsid_df  = data.frame(ensgid_rsid,ensgid_split_df)
    dim(ensgid_rsid_df) %>% print

    # Merge 1: gene annotations
    paste0('  Merging biomaRt annotations.. ') %>% cat
    ensgid_rsid_df1 = merge(ensgid_rsid_df,gene_ensg,by='ensgid',all=T) %>% unique
    dim(ensgid_rsid_df1) %>% print

    # Merge 2: eqtl, hic, nearest
    paste0('  Merging eQTL, Hi-C, nearest genes... ') %>% cat
    col_names1 = basename(paths) %>% file_path_sans_ext
    if(!is.null(nearest)) {
        col_names = c(col_names1,'nearest_genes')
    } else col_names = col_names1
    n = length(snp_ensg_li)
    for(i in 1:n) {
        m = length(snp_ensg_li[[i]])
        df = data.frame(
            ensg_rsid = snp_ensg_li[[i]],
            t_f = rep('TRUE',m)
        )
        colnames(df) = c('ensgid_rsid',col_names[i])

        # Merge the df to the ensgid_rsid_df1
        ensgid_rsid_df1 = merge(ensgid_rsid_df1,df,by='ensgid_rsid',all=T)
    }
    ensgid_rsid_df2 = ensgid_rsid_df1 %>% unique
    dim(ensgid_rsid_df2) %>% print

    # Merge 3: gtex annotations
    paste0('  Merging GTEx eQTL SNP slopes... ') %>% cat
    ensgid_rsid_df3 = merge(ensgid_rsid_df2,tb_li[[1]],by='ensgid_rsid',all=T) %>% unique
    dim(ensgid_rsid_df3) %>% print

    # Save as CSV file
    f_name1 = paste0(out,'/',dir_name,'_pairs.csv')
    write.csv(ensgid_rsid_df3,f_name1,row.names=F)
    paste0('  Write file: ',f_name1,'\n') %>% cat
    return()

    ## Legacy codes: biomaRt enstid searching ##
    # Extract Enstids
    #paste0('\n  Extracting Enstids... ') %>% cat
    #enstids = enst_li %>% unlist %>% unique
    #paste0(length(enstids),'.. ') %>% cat

    # Search biomaRt for enstids
    #paste0('biomaRt... ') %>% cat
    #ensembl   = useMart('ensembl',dataset='hsapiens_gene_ensembl')
    #gene_attr1 = c('ensembl_transcript_id','ensembl_gene_id','hgnc_symbol','description')
    #gene_enst  = getBM(
    #    attributes = gene_attr1,
    #    filters = "ensembl_transcript_id",
    #    values = enstids,
    #    mart = ensembl
    #) %>% unique
    #colnames(gene_enst)[1:2] = c('enstid','ensgid')
    #dim(gene_enst) %>% print
    #head(gene_enst) %>% print
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
    paste0('done\n') %>% cat

    # Spread GTEx pair by tissue
    paste0('  Spreading GTEx data.. ') %>% cat
    ensg_rsid = paste0(ensgid,'-',gtex$rsid)
    gtex_ = data.frame(
        #ensg_rsid = ensg_rsid,
        rsid = gtex$rsid,
        ensgid = ensgid,
        tissue = gtex$tissue,
        slope = gtex$slope
    )
    gtex_reshape = spread(gtex_,key=tissue,value=slope,fill='-')
    n = ncol(gtex_reshape)
    dim(gtex_reshape) %>% print

    # Save as TSV file
    gene_n = unique(gtex_reshape$ensgid) %>% length
    f_name = paste0(out,'/gtex_eqtl_',gene_n,'.tsv')
    write.table(gtex_reshape,f_name,sep='\t',row.names=F,quote=F)
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
    if(length(args$help)>0) {    help    = args$help
    } else                       help    = FALSE
    if(help) { cat(help_message); quit() }

    # Global arguments
    if(length(args$base)>0)      b_path  = args$base
    if(length(args$out)>0)       out     = args$out
    if(length(args$debug)>0) {   debug   = args$debug
    } else                       debug   = FALSE

    # Required arguments
    if(length(args$bed)>0) {     bed     = args$bed
    } else                       bed     = FALSE
    if(length(args$nearest)>0) { nearest = args$nearest
    } else                       nearest = NULL
    if(length(args$fdr)>0) {     fdr     = args$fdr %>% as.numeric
    } else                       fdr     = 0.05

    # Run function
    source('src/pdtime.r'); t0=Sys.time()
    if(args$dbgene == 'hic_pair') {
        hic_pair(b_path,out,bed,debug)
    } else if(args$dbgene =='gtex_pair') {
        gtex_pair(b_path,out,debug)
    } else if(args$dbgene == 'summary') {
        summary(b_path,out,nearest,debug)
    } else if(args$dbgene == 'david_go') {
        david_go(b_path,out,fdr,debug)
    } else {
        paste0('[Error] There is no such function "',args$dbgene,'" in gene: ',
            paste0(args$dbgene,collapse=', '),'\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}