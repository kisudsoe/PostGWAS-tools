help_message = '
utils, v2020-08-01
This is a function call for utils such as generating bedtools command.

Usage:
    Rscript postgwas-exe.r --utils genetraits --gwas <gwas file> --phewas <phewas file> --out <out folder>
    Rscript postgwas-exe.r --utils gene --base <base file> --trait <trait file> --out <out folder>
    Rscript postgwas-exe.r --utils bash --base <base file> --ann <anntation folder> --out <out folder>


Function:
    genetraits  Preparing GWAS and PheWAS traits as rsid-trait-gene relationships.
    gene        Finding traits and their risk SNPs from input gene list.
    bash        Generating bash command scripts to run bedtools.

Global arguments:
    --base      <base file path>
                Mendatory. For bash function.
    --out       <out folder path>
                Mendatory. For bash function.

Required arguments:
    --ann       <Functional annotation folder path>
                Mendatory. For bash function.
    --gwas      <GWAS file path>
                Mendatory. For genetraits function.
    --phewas    <PheWAS file path>
                Mendatory. For genetraits function.
    --age       <age-related traits meta file path>
                Mendatory. For genetraits function.
    --trait     <Trait file path>
                Mendatory. For gene function.
'


## Load global libraries ##
suppressMessages(library(dplyr))


## Functions Start ##
bash_script = function(
    b_path   = NULL, # Input SNP TSV file
    ann_path = NULL, # Folder path for storing functional annotation data
    out      = NULL, # Output folder path 
    debug    = FALSE
) {
    # Preparing...
    paste0('\n** Run function: utils.r/bash_script... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    out_genome = paste0(out,'/genome_dist')
    ifelse(!dir.exists(out_genome), dir.create(out_genome),'')
    out_roadmap = paste0(out,'/roadmap_dist')
    ifelse(!dir.exists(out_roadmap), dir.create(out_roadmap),'')
    'Ready\n' %>% cat

    bedtools_sort = 'bedtools sort -i '
    bedtools_closest_a = 'bedtools closest -d -a '
    bedtools_closest_b = ' -b '
    bedtools_out = ' > '

    # gene_dist files
    genome_dist_f = c(
        paste0(ann_path,'/ensembl_gene_hg19.bed'),
        paste0(ann_path,'/wgEncodeRegTfbsClusteredV3.bed'),
        paste0(ann_path,'/ucsc_annot.bed')
    )

    # roadmap files
    roadmap_f1 = list.files(paste0(ann_path,'/roadmap_enh'),full.names=F)
    roadmap_f1 = c(roadmap_f1,'roadmap_enh_merge.bed') # Add roadmap merge file
    roadmap_f = list.files(paste0(ann_path,'/roadmap_enh'),full.names=T)
    roadmap_f = c(roadmap_f,paste0(ann_path,'/roadmap_enh_merge.bed')) # Add roadmap merge file

    # Generate genome_dist bash script
    genome_dist_out = c(
        'nearest_gene.tsv',
        'encode_tfbs.tsv',
        'ucsc_annot.tsv'
    )

    out1 = paste0(out,'/genome_dist/',genome_dist_out)

    bash1 = paste0(
        bedtools_sort, genome_dist_f, ' | ',
        bedtools_closest_a, b_path,
        bedtools_closest_b, 'stdin',
        bedtools_out, out1
    )

    # Generate roadmap bash script
    out2 = paste0(out,'/roadmap_dist/',tools::file_path_sans_ext(roadmap_f1),'.tsv')

    bash2 = paste0(
        bedtools_sort, roadmap_f, ' | ',
        bedtools_closest_a, b_path,
        bedtools_closest_b, 'stdin',
        bedtools_out, out2
    )

    # Combine and save the bash scripts
    tag1 = 'printf "\n  1) Genome_dist.. "'
    tag2 = 'printf "done\n  2) Roadmap_dist.. \n"'
    tag3 = 'printf "done\n"'
    bash = c(tag1,bash1,tag2,bash2,tag3)
    out_base = tools::file_path_sans_ext(basename(out))
    f_name = paste0('dist_',out_base,'.sh')
    out_f = file(f_name,"wb") # Set file as Unix type
    write.table(bash,file=out_f,row.names=F,col.names=F,quote=F,sep='')
    paste0('Write bash file: ',f_name,'\n') %>% cat
}


gene = function(
    b_path = NULL,
    out    = NULL,
    f_trait    = NULL,
    debug  = FALSE
) {
    # Preparing...
    paste0('\n** Run function: utils.r/bash_script... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    'Done\n' %>% cat

    # Read files
    paste0('Read, ',b_path,' = ') %>% cat
    input = read.delim(b_path)
    dim(input) %>% print

    paste0('Read, ',f_trait,' = ') %>% cat
    trait = read.delim(f_trait)
    dim(trait) %>% print

    # Subset  by ensgid
    paste0('  Query ',nrow(input),' ensgids, found # = ') %>% cat
    trait_sub = subset(trait,Ensgid %in% input$Ensgid)
    dim(trait_sub) %>% print

    # Write TSV file
    n = trait_sub$Rsid %>% unique %>% length
    f_name = paste0(out,'/input_snp_',n,'.tsv')
    write.table(trait_sub,f_name,sep='\t',quote=F,row.names=F)
    paste0('Write TSV file: ', f_name,'\n') %>% cat
}


gene_traits = function(
    f_gwas   = NULL,
    f_phewas = NULL,
    f_age    = NULL,
    out      = NULL,
    debug    = FASLE
) {
    # Function specific library
    suppressMessages(library(biomaRt))  # download SNP coordiantes (hg19/hg38)

    # Preparing...
    paste0('\n** Run function: utils.r/bash_script... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    'Done\n' %>% cat

    # Read files
    paste0('Read, ',f_gwas,' = ') %>% cat
    gwas = read.delim(f_gwas)
    dim(gwas) %>% print

    paste0('Read, ',f_phewas,' = ') %>% cat
    phewas = read.csv(f_phewas)
    dim(phewas) %>% print

    paste0('Read, ',f_age,' = ') %>% cat
    age = read.delim(f_age)
    dim(age) %>% print

    # Select gwas columns
    paste0('  gwas selection... ') %>% cat
    gwas_sel = gwas %>% dplyr::select(
        'UPSTREAM_GENE_ID','DOWNSTREAM_GENE_ID','SNP_GENE_IDS',
        'SNPS','P.VALUE','MAPPED_TRAIT'
    ) %>% unique
    paste0(ncol(gwas_sel),'.. ') %>% cat

    gwas_li = apply(gwas_sel,1,function(x) {
        ensgids1 = x[1:3][!x[1:3] %in% '']
        Ensgid_vec = lapply(ensgids1,function(x1) {
            strsplit(x1,'\\, ')[[1]]
        }) %>% unlist %>% unique
        a_trait = paste0('[gwas]_',x[6])
        n = length(Ensgid_vec)
        if(n>0) {
            data.frame(
                Rsid    = rep(x[4],n),
                P_value = rep(x[5],n),
                Trait   = rep(a_trait,n),
                Ensgid  = Ensgid_vec
            )
        } else return(NULL)
    }) %>% unique
    gwas_df = data.table::rbindlist(gwas_li)
    dim(gwas_df) %>% print

    # Select phewas columns
    paste0('  phewas selection = ') %>% cat
    phewas_sel = phewas %>% dplyr::select(
        'snp','phewas_phenotype','p_value','gene_name','phewas_code'
    ) %>% unique
    a_trait = paste0(
        '[phecode ',phewas_sel$phewas_code,
        ']_',phewas_sel$phewas_phenotype
    )
    phewas_df = data.frame(
        Rsid    = phewas_sel$snp,
        P_value = phewas_sel$p_value,
        Trait   = a_trait,
        Symbol  = phewas_sel$gene_name
    )
    dim(phewas_df) %>% print


    # Search biomaRt (hg38) from GWAS Ensgid
    paste0('\n  search gwas ensgid biomaRt hg38 = ') %>% cat
    ensgids = gwas_df$Ensgid %>% unique
    ensgid_df = ensgids %>% as.data.frame
    colnames(ensgid_df) = 'Ensgid'
    hg38_ensg = useMart(
        biomart='ENSEMBL_MART_ENSEMBL',
        dataset='hsapiens_gene_ensembl'
    )
    ensgid_attr = c('ensembl_gene_id','hgnc_symbol')
    ensgid_hg38 = getBM(
        attributes = ensgid_attr,
        filters    = 'ensembl_gene_id',
        values     = ensgids,
        mart       = hg38_ensg
    ) %>% unique
    dim(ensgid_hg38) %>% print

    paste0('  Among queried ',length(ensgids),', unfound Ensgid # = ') %>% cat
    ensgid_merge = merge(
        ensgid_df, ensgid_hg38,
        by.x='Ensgid', by.y='ensembl_gene_id', all.x=T
    ) %>% unique
    #which_na = is.na(ensgid_merge$hgnc_symbol) %>% which
    symb = ensgid_merge$hgnc_symbol
    which_nosymbol = symb[symb=='']
    length(which_nosymbol) %>% print

    paste0('  Merge search result = ') %>% cat
    colnames(ensgid_merge) = c('Ensgid','Symbol')
    gwas_merge = merge(gwas_df,ensgid_merge,by='Ensgid',all.x=T)
    dim(gwas_merge) %>% print


    # Search biomaRt (hg38) from Phewas Symbol
    paste0('\n  search phewas symbol biomaRt hg38 = ') %>% cat
    symbols = phewas_df$Symbol %>% unique
    symbol_df = symbols %>% as.data.frame
    colnames(symbol_df) = 'Symbol'
    hg38_symb = useMart(
        biomart='ENSEMBL_MART_ENSEMBL',
        dataset='hsapiens_gene_ensembl'
    )
    symbol_attr = c('ensembl_gene_id','hgnc_symbol')
    symbol_hg38 = getBM(
        attributes = symbol_attr,
        filters    = 'hgnc_symbol',
        values     = symbols,
        mart       = hg38_symb
    ) %>% unique
    dim(symbol_hg38) %>% print

    paste0('  Among queried ',length(symbols),', unfound Symbol # = ') %>% cat
    symbol_merge = merge(
        symbol_df, symbol_hg38,
        by.x='Symbol', by.y='hgnc_symbol', all.x=T
    ) %>% unique
    ensg = symbol_merge$ensembl_gene_id
    which_noensgid = ensg[ensg=='']
    length(which_noensgid) %>% print

    paste0('  Merge search result = ') %>% cat
    colnames(symbol_merge) = c('Symbol','Ensgid')
    phewas_merge = merge(phewas_df,symbol_merge,by='Symbol',all.x=T)
    dim(phewas_merge) %>% print


    # rbind gwas and phewas dfs
    paste0('\n  Bind gwas and phewas = ') %>% cat
    gwas_merge = gwas_merge %>% dplyr::select('Ensgid','Symbol','Rsid','Trait','P_value')
    phewas_merge = phewas_merge %>% dplyr::select('Ensgid','Symbol','Rsid','Trait','P_value')
    gwas_phewas = rbind(gwas_merge,phewas_merge)
    dim(gwas_phewas) %>% print


    # Merge age-related traits
    paste0('  Merge age-realted traits = ') %>% cat
    age_sub = age %>% dplyr::select('Trait','Category')
    gwas_phewas_trait = merge(gwas_phewas,age_sub,by='Trait',all.x=T)
    dim(gwas_phewas_trait) %>% print


    # Write TSV file
    f_name = paste0(out,'/gwas_phewas_ensgid_rsid_trait_pval.tsv')
    write.table(gwas_phewas_trait,f_name,sep='\t',row.names=F,quote=F)
    paste0('Write TSV file, ',f_name,'\n') %>% cat
}


## INIT Function ##
utils = function(
    args = NULL
) {
    # Help message
    if(length(args$help)>0) {     help     = args$help
    } else                        help     = FALSE
    if(help) { cat(help_message); quit() }

    # Global arguments
    if(length(args$base)>0)       b_path   = args$base
    if(length(args$out)>0)        out      = args$out
    if(length(args$debug)>0) {    debug    = args$debug
    } else                        debug    = FALSE

    # Required arguments
    if(length(args$ann)>0)        f_ann    = args$ann
    if(length(args$gwas)>0)       f_gwas   = args$gwas
    if(length(args$phewas)>0)     f_phewas = args$phewas
    if(length(args$trait)>0) {      f_trait    = args$trait
    } else                        debug    = NULL

    # Run function
    source('src/pdtime.r'); t0=Sys.time()
    if(args$utils=='genetraits') {
        gene_traits(f_gwas,f_phewas,f_age,out,debug)
    } else if(args$utils=='gene') {
        gene(b_path,out,f_trait,debug)
    } else if(args$utils=='bash') {
        bash_script(b_path,f_ann,out,debug)
    } else {
        paste0('[Error] There is no such function "',args$utils,'" in utils.\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n\n') %>% cat
}