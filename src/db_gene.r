help_message = '
Version: v2020-04-29
This is a function call for gene analysis to compile the Hi-C and eQTL data.

Usage: Rscript postgwas-exe.r --dbgene <function> --base <base files> --out <out folder> [options]

Functions:
    hic_pair    Extract Hi-C linked SNP-Gene pairs.
    gtex_pair   Extract eQTL linked SNP-Gene pairs.
    summary     Summarizing GWAS SNPs-Gene pairs.
    david_go    Summarizing DAVID GO analysis result.
    pivot_gene  Pivotting the gene summary pair table to gene level summary.

Global arguments:
    --base      <base files>
                For "hic_pair" function, two input files are:
                  [1] gwas_dist file, [2] gene_dist file
    --out       <out folder>
                Out folder path is mendatory. Default is "data" folder.
                If the "bed=TRUE" option put in "hic_pair", out folder could be two.

Required arguments:
    --bed       <default:FALSE>
                An optional argument for the "hic_pair" function.
                To save as BED format file.
    --nearest   <nearest gene summary file path>
                An optional argument for the "summary" function.
                Add nearest gene summary table to the hic gene summary.
    --criteria  <default:0.05>
                An optional argument for the "david_go" function.
    --stat      <default:fdr>
                An optional argument for the "david_go" function.
                Either --fdr or --pval have to choose.
    --dataset   An optional argument for the "david_go" function.
                Add filtering dataset name.
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
pivot_gene = function(
    summ_path = NULL,   # Input gene summary pair TSV table
    out       = 'data', # Out folder path
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_gene.r/pivot_gene... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),''); 'ready\n' %>% cat

    # Read gene summary pair file
    paste0('  Read summary_gene_pair TSV file... ') %>% cat
    pair = read.delim(summ_path)
    dim(pair) %>% print

    # Extract Ensgids
    paste0('  Extract ensgids... ') %>% cat
    ensgid_pair = pair[,1]
    pair_val = pair[,-1]
    n = ncol(pair_val)
    ensgids_li = lapply(c(1:n),function(i) {
        which_i  = which(pair_val[,i]==TRUE)
        ensgid_pair[which_i] %>% unique
    })

    # Generate union table
    paste0('union table.. ') %>% cat
    ensgid_ = ensgid_pair %>% unique
    m = length(ensgids_li)
    df = NULL
    for(i in 1:m) {
        df = cbind(df,ensgids_li[[i]][match(ensgid_,ensgids_li[[i]])])
    }
    df1 = (df!='')           # Transcform values to TRUE, if ID exists.
    df1[is.na(df)] = FALSE   # Transform NA to FALSE value
    df1 = as.data.frame(df1) # Make 'df1' to data.frame
    colnames(df1) = colnames(pair_val)

    df2 = data.frame(
        ensgid = ensgid_,
        df1
    )
    dim(df2) %>% print
    
    # Save as CSV file
    f_name = paste0(out,'/summary_gene_pivot.csv')
    write.csv(df2,f_name,row.names=F)
    paste0('  Write file: ',f_name,'\n') %>% cat
}

david_go = function(
    go_path  = NULL,   # Input DAVID GO analysis result file path
    out      = 'data', # Out folder path
    criteria = 0.05,   # Optional: Significance criteria
    stat     = 'fdr',  # Optional: 'fdr','pval'
    dataset  = NULL,   # Optional: Filtering dataset name
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_gene.r/david_go... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),''); 'ready\n' %>% cat

    # Read go result file
    paste0('  Read DAVID result file\t') %>% cat
    go = read.delim(go_path)
    dim(go) %>% print

    if(!is.null(dataset)) {
        paste0('  Filtering dataset = ',dataset,'... ') %>% cat
        go = subset(go,Geneset==dataset)
        dim(go) %>% print
    }
    
    if(stat=='fdr') {
        paste0('\n  FDR criteria\t= ',criteria,'\n') %>% cat
        go_ = subset(go,FDR<criteria)
    } else if(stat=='pval') {
        paste0('\n  Pval criteria\t= ',criteria,'\n') %>% cat
        go_ = subset(go,PValue<criteria)
    } else {
        '\n** [ERROR] Need fdr/pval criteria.\n' %>% cat
        quit()
    }
    paste0('  Filtered GO terms\t= ') %>% cat
    dim(go_) %>% print
    
    paste0('\n  Extract the gene list... ') %>% cat
    genes_li = lapply(go_$Genes,function(x) { strsplit(x%>%as.character,"\\, ")[[1]] })
    names(genes_li) = go_$Term
    genes_n = unlist(genes_li) %>% unique %>% length
    paste0(genes_n,'.. done\n') %>% cat
    print(genes_li)
    
    paste0('  Venn analysis... ') %>% cat
    source('src/venn_analysis.r')
    venn_tf = venn_analysis(genes_li,wfig=F,wfile=F)
    dim(venn_tf) %>% print
    
    if(is.null(dataset)) {
        f_name = paste0(out,'/go_analysis.',stat,'-',criteria,'.csv')
    } else f_name = paste0(out,'/go_analysis.',dataset,'.',stat,'-',criteria,'.csv')
    write.csv(venn_tf,f_name,quote=F,row.names=F)
    paste0('  Write file: ',f_name,'\n') %>% cat
}

summary_gene = function(
    f_paths = NULL,   # Input snp-gene pair files
    out     = 'data', # Out folder path
    nearest = NULL,   # Nearest gene summary table
    debug
) {
    # Load function-specific library
    suppressMessages(library(biomaRt))
    suppressMessages(library(tools))

    # Preparing...
    paste0('\n** Run function: db_gene.r/summary_gene... ') %>% cat
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
    paste0('  Read gene-snp pair files... ') %>% cat
    ensg_li = list(); snp_ensg_li = list(); tb_li = list()
    for(i in 1:n) {
        #paste0('  ',i,' ',paths[i],'\n') %>% cat # for debug
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
    'done\n' %>% cat

    # Read nearest gene file
    if(!is.null(nearest)) {
        paste0('  Read nearest gene file... ') %>% cat
        tb = read.csv(nearest,stringsAsFactors=F)
        colnames(tb)[2] = 'ensgid'
        tb_ncol = ncol(tb)
        n = length(ensg_li)
        ensg_li[[n+1]] = tb$ensgid
        snp_ensg = paste0(tb$ensgid,'-',tb$rsid)
        tb$ensgid_rsid = snp_ensg
        snp_ensg_li[[n+1]] = snp_ensg
        tb_li[[n+1]] = tb
        dim(tb) %>% print
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
            rsid = row[2] )
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
    df = NULL
    for(i in 1:n) {
        df = cbind(df,snp_ensg_li[[i]][match(ensgid_rsid,snp_ensg_li[[i]])])
    }
    df1 = (df!='')           # Transcform values to TRUE, if ID exists.
    df1[is.na(df)] = FALSE   # Transform NA to FALSE value
    df1 = as.data.frame(df1) # Make 'df1' to data.frame from
    colnames(df1) = col_names

    df2 = data.frame(
        ensgid_rsid = ensgid_rsid,
        df1
    )
    ensgid_rsid_df1 = merge(ensgid_rsid_df1,df2,by='ensgid_rsid',all=T)
    ensgid_rsid_df2 = ensgid_rsid_df1 %>% unique
    dim(ensgid_rsid_df2) %>% print

    # Merge 3: gtex annotations
    paste0('  Merging GTEx eQTL SNP slopes... ') %>% cat
    ensgid_rsid_df3 = merge(ensgid_rsid_df2,tb_li[[1]],by='ensgid_rsid',all=T) %>% unique
    dim(ensgid_rsid_df3) %>% print

    # Save as CSV file
    f_name1 = paste0(out,'/',dir_name,'_pairs.csv')
    write.csv(ensgid_rsid_df3,f_name1,row.names=F)
    paste0('\n  Write file: ',f_name1,'\n') %>% cat
    return()

    ## Legacy codes: biomaRt enstid searching ##
    ## Don't forget to add enstid searching module later
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
    bed     = FALSE,  # Return as BED format
    debug
) {
    # Load function-specific library
    suppressMessages(library(biomaRt))
    if(bed=='TRUE') { bed = TRUE
    } else if(bed=='FALSE') bed = FALSE

    # Prepare...
    paste0('\n** Run function: db_gene.r/hic... ') %>% cat
    out_n = length(out)
    if(out_n==1) {
        ifelse(!dir.exists(out), dir.create(out),'')
    } else if(out_n==2) {
        ifelse(!dir.exists(out[1]), dir.create(out[1]),'')
        ifelse(!dir.exists(out[2]), dir.create(out[2]),'')
    } else {
        paste0('[ERROR] out folder should be one. Or with "bed=TRUE" option, it can be two.')
        quit()
    }
    'ready\n' %>% cat

    # Read base files and filter by distance=0
    colnm1 = c('chr','start','end','name','hic_chr','hic_start','hic_end','hic_loop')
    colnm2 = c('dist')
    paste0('  ',f_paths[1],', length= ') %>% cat
    gwas = data.table::fread(f_paths[1],header=F)
    k = ncol(gwas)
    colnames(gwas)[1:8] = colnm1
    colnames(gwas)[k] = colnm2
    paste0(nrow(gwas)) %>% cat
    gwas0 = subset(gwas,dist==0) %>% unique
    gwas_0 = gwas0[,c(4,8)]
    paste0(',\toverlap= ',nrow(gwas_0),'\n') %>% cat

    paste0('  ',f_paths[2],', length= ') %>% cat
    gene = data.table::fread(f_paths[2],header=F)
    k = ncol(gene)
    colnames(gene)[1:8] = colnm1
    colnames(gene)[k] = colnm2
    paste0(nrow(gene)) %>% cat
    gene0 = subset(gene,dist==0) %>% unique
    gene_0 = gene0[,c(4,8)]
    paste0(',\toverlap= ',nrow(gene_0),'\n') %>% cat

    # Merge the gwas and gene data
    paste0('\n  Process gwas_loop.. ') %>% cat
    gwas_loop = gwas_0$hic_loop
    gwas_loop_id = lapply(gwas_loop,function(x) strsplit(x,"\\.")[[1]][1]) %>% unlist
    paste0(length(gwas_loop_id),'.. done\n') %>% cat

    paste0('  Process gene_loop.. ') %>% cat
    gene_loop = gene_0$hic_loop
    gene_loop_id = lapply(gene_loop,function(x) strsplit(x,"\\.")[[1]][1]) %>% unlist
    paste0(length(gene_loop_id),'.. done\n') %>% cat

    paste0('  Process merge for TSV.. ') %>% cat
    n = nrow(gwas_0); paste0(n,'.. ') %>% cat
    tsv_li = lapply(c(1:n),function(i) {
        gwas_i = gwas0[i,]
        which_row = which(gene_loop_id==gwas_loop_id[i])
        loop = gene_loop[which_row]
        loop_g = gene_0[which_row,]
        if(nrow(loop_g)>0) {
            loop_g_ = subset(loop_g,hic_loop!=gwas_0$hic_loop[i])
            if(nrow(loop_g_)>0) {
                m = nrow(loop_g_)
                genes = lapply(loop_g_$name,function(x)
                    strsplit(x,"\\_")[[1]][1]) %>% unlist
                tags  = lapply(loop_g_$name,function(x) {
                    tag_v = strsplit(x,"\\_")[[1]]
                    l = length(tag_v)
                    paste0(tag_v[2:l],collapse="_") %>% return
                }) %>% unlist
                out = data.frame(
                    rsid      = rep(gwas_i$name,m),
                    ensgid    = genes,
                    gene_tag  = tags,
                    gwas_loop = rep(gwas_i$hic_loop,m),
                    gene_loop = loop_g_$hic_loop
                )
                return(out)
            } else return(NULL)
        } else return(NULL)
    })
    tsv_df = data.table::rbindlist(tsv_li) %>% unique
    dim(tsv_df) %>% print

    # Save as TSV file
    file_nm = tools::file_path_sans_ext(basename(f_paths[1]))
    gene_n = unique(tsv_df$ensgid) %>% length
    f_name = paste0(out[1],'/hic_',file_nm,'_',gene_n,'.tsv')
    write.table(tsv_df,f_name,sep='\t',row.names=F,quote=F)
    paste0('  Write file: ',f_name,'\n') %>% cat

    if(bed) {
        paste0('\n  Process extract for BED.. ') %>% cat
        paste0(n,'.. ') %>% cat
        bed_li = lapply(c(1:n),function(i) {
            gwas_i = gwas0[i,]
            which_row = which(gene_loop_id==gwas_loop_id[i])
            loop = gene_loop[which_row]
            loop_g = gene_0[which_row,]
            if(nrow(loop_g)>0) {
                loop_g_ = subset(loop_g,hic_loop!=gwas_0$hic_loop[i])
                if(nrow(loop_g_)>0) {
                    out = data.frame(
                        chr   = gwas_i$chr,
                        start = gwas_i$start,
                        end   = gwas_i$end,
                        rsid  = gwas_i$name
                    )
                    return(out)
                } else return(NULL)
            } else return(NULL)
        })
        bed_df = data.table::rbindlist(bed_li) %>% unique
        dim(bed_df) %>% print

        snp_n  = unique(bed_df$rsid) %>% length
        f_name = paste0(out[2],'/snp_hic_',file_nm,'_',snp_n,'.bed')
        write.table(bed_df,f_name,sep='\t',row.names=F,col.names=F,quote=F)
        paste0('  Write file: ',f_name,'\n') %>% cat
    }
}

## __INIT__ ##
db_gene = function(
    args = NULL
) {
    # Get help
    if(length(args$help)>0) {     help     = args$help
    } else                        help     = FALSE
    if(help) { cat(help_message); quit() }

    # Global arguments
    if(length(args$base)>0)       b_path   = args$base
    if(length(args$out)>0)        out      = args$out
    if(length(args$debug)>0) {    debug    = args$debug
    } else                        debug    = FALSE

    # Required arguments
    if(length(args$bed)>0) {      bed      = args$bed
    } else                        bed      = FALSE
    if(length(args$nearest)>0) {  nearest  = args$nearest
    } else                        nearest  = NULL
    if(length(args$criteria)>0) { criteria = args$criteria %>% as.numeric
    } else                        criteria = 0.05
    if(length(args$stat)>0) {     stat     = args$stat
    } else                        stat     = 'fdr'
    if(length(args$dataset)>0) {  dataset  = args$dataset
    } else                        dataset  = NULL

    # Run function
    source('src/pdtime.r'); t0=Sys.time()
    if(args$dbgene == 'hic_pair') {
        hic_pair(b_path,out,bed,debug)
    } else if(args$dbgene =='gtex_pair') {
        gtex_pair(b_path,out,debug)
    } else if(args$dbgene == 'summary') {
        summary_gene(b_path,out,nearest,debug)
    } else if(args$dbgene == 'david_go') {
        david_go(b_path,out,criteria,stat,dataset,debug)
    } else if(args$dbgene == 'pivot_gene') {
        pivot_gene(b_path,out,debug)
    } else {
        paste0('[Error] There is no such function "',args$dbgene,'" in gene: ',
            paste0(args$dbgene,collapse=', '),'\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}