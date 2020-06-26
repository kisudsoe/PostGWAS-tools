help_message = '
db_download, v2020-02-27
This is a function call for downloading databases
    Roadmap, ENCODE, RegulomeDB, GTEx v8, and lncRNASNP2

Usage:
    Rscript postgwas-exe.r --dbdown roadmap --out <out folder>
    Rscript postgwas-exe.r --dbdown encode --out <out folder>
    Rscript postgwas-exe.r --dbdown regulome --out <out folder>
    Rscript postgwas-exe.r --dbdown gtex --out <out folder>
    Rscript postgwas-exe.r --dbdown lncrna --out <out folder>
    Rscript postgwas-exe.r --dbdown genes --out <out folder> --hg hg19
    Rscript postgwas-exe.r --dbdown genes --out <out folder> --hg hg38


Functions:
    roadmap   Downloading Roadmap data (hg19).
    encode    Downloading ENCODE data (hg19).
    regulome  Downloading RegulomeDB data (≥2b, hg19).
    gtex      Downloading GTEx v8 data (hg38).
    lncrna    Downloading lncRNASNP2 data (hg38).
    genes     Downloading Ensembl Biomart Gene coordinates (hg19/hg38).

Global arguments:
    --out     <out folder>
              Download folder path is mendatory. Default is "db" folder.

Required arguments:
    --hg      <hg19/hg38>
              A required argument for the "genes" function. Choose one human genome version.
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
biomart_gene = function(
    out = 'data', # Download folder path
    hg  = 'hg19'  # Human genome version. Choose either 'hg19' or 'hg38'
) {
    # Function specific library
    suppressMessages(library(biomaRt))

    # Get Ensembl BiomaRt
    paste0('\n** Run function: db_download.r/biomart_gene...\n') %>% cat
    if(hg=='hg19') { # Ensembl biomart (grch37)
        hg_gene = useMart(
            biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
            dataset="hsapiens_gene_ensembl", path="/biomart/martservice")
    } else if(hg=='hg38') { # Ensembl biomart (grch38)
        hg_gene = useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
    } else {
        paste0('  Please choose the correct --hg options: either "hg19" or "hg38".\n') %>% cat
        return()
    }
    gene_attr = c('chromosome_name','start_position','end_position','ensembl_gene_id','external_gene_name')
    genes = getBM(attributes = gene_attr,
                filters    = '',
                values     = '',
                mart       = hg_gene)
    genes = as.data.frame(genes)
    paste0('  BiomaRt table, dim\t= ') %>% cat; dim(genes) %>% print

    # Save as a TSV file
    f_name1 = paste0(out,'/ensembl_gene_ann_',hg,'.tsv')
    write.table(genes,f_name1,row.names=F,col.names=T,quote=F,sep='\t')
    cat(paste0('  File write: ',f_name1,'\n'))

    # Filter genes by chromosomes
    genes_sub = subset(genes,chromosome_name %in% c(1:22,'X','Y'))
    chr_ = paste0('chr',genes_sub[,1])
    genes_ = data.frame(chr_,genes_sub[,2:5])
    colnames(genes_) = c('chr','start','end','ENSGid','Symbol')
    paste0('\n  Filtered table, dim\t= ') %>% cat; dim(genes_) %>% print

    # Save as a BED file
    f_name2 = paste0(out,'/ensembl_gene_',hg,'.bed')
    write.table(genes_[,1:4],f_name2,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('  File write: ',f_name2,'\n\n') %>% cat
}

roadmap_down = function(
    out = 'db' # Download folder path
) {
    # Function specific library
    source('src/saveasrds.r')

    # Download from Roadmap
    paste0('\n** Run function: db_download.r/roadmap_down...\n') %>% cat
    mkdir_out(out)

    cid = as.character(formatC(c(1:129),width=3,flag='0'))
    #urls = paste0('https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E',cid,'_25_imputed12marks_hg38lift_dense.bed.gz')
    urls = paste0('https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E',cid,'_25_imputed12marks_dense.bed.gz')
    urls = c(urls,'https://raw.githubusercontent.com/Bioconductor/BioC2015Introduction/master/inst/extdata/epi_metadata.txt')
    o=lapply(urls,function(url) {
        f_name = paste0(out,'/',basename(url))
        tb = try(download.file(url,destfile=f_name)) # debug no file error
        if("try-error" %in% class(tb)) tb=NULL
        cat(paste0('\n',f_name,'\n'))
        try(R.utils::gunzip(f_name))
        f_name2 = tools::file_path_sans_ext(f_name)
        cat(paste0('  Convert: ',f_name2,'\n'))
        try(bedasrds(f_name2))
    })
}

encode_down = function(
    out = 'db' # Download folder path
) {
    # Download from ENCODE
    paste0('\n** Run function: db_download.r/encode_down...\n') %>% cat
    mkdir_out(out)

    url = paste0('http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz')
    f_name = paste0(out,'/',basename(url))
    tb = try(download.file(url,destfile=f_name)) # debug no file error
    if("try-error" %in% class(tb)) tb=NULL
    cat(paste0('\n',f_name,'\n'))
}

regulome_down = function(
    out = 'db' # Download folder path
) {
    # Function specific library
    source('src/saveasrds.r')

    # Download from RegulomeDB
    paste0('\n** Run function: db_download.r/regulome_down...\n') %>% cat
    mkdir_out(out)

    url_base = 'http://legacy.regulomedb.org/downloads/RegulomeDB.' # url changed 2020-02-27
    url_li = list(
        #total= c(url=paste0(url_base,'dbSNP141.txt.gz'),
        #    path=paste0(dir,'dbSNP141.txt.gz')),
        cat1 = c(url=paste0(url_base,'dbSNP132.Category1.txt.gz'),
            path=paste0(out,'/dbSNP132.Category1.txt.gz')),
        cat2 = c(url=paste0(url_base,'dbSNP132.Category2.txt.gz'),
            path=paste0(out,'/dbSNP132.Category2.txt.gz'))
    )
    paste0('  Download RegulomeDB data\n\n') %>% cat
    for(i in 1:length(url_li)) {
        tb = try(download.file(url_li[[i]][1],destfile=url_li[[i]][2]))
        if('try-error' %in% class(tb)) stop('The file url seems to be changed. Please check the lncRNASNP2 homepasge: http://www.regulomedb.org/')

        # Save as RDS file
        saveasrds(url_li[[i]][2])
    }
}

gtex_down = function(
    out = 'db' # Download folder path
) {
    # Function specific library
    suppressMessages(library(tools))
    source('src/saveasrds.r')
    source('src/pdtime.r'); t0  = Sys.time()

    # Download from GTEx v8
    paste0('\n** Run function: db_download.r/gtex_down...\n') %>% cat
    mkdir_out(out)

    url_base = 'https://storage.googleapis.com/gtex_analysis_v8/'
    url_li = list(
        eqtl = c(
            url=paste0(url_base,'single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar'),
            path=paste0(out,'/GTEx_Analysis_v8_eQTL.tar')),
        snpid = c(
            url=paste0(url_base,'reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz'),
            path=paste0(out,'/GTEx_Analysis_2017-06-05_v8_lookup_table.txt.gz'))
    )
    paste0('  Download GTEx data\n\n') %>% cat
    f_target = url_li[[1]][2] %>% basename %>% file_path_sans_ext
    for(i in 1:length(url_li)) {
        url = url_li[[i]][1]
        name = url_li[[i]][2]
        paste0('  ',name,'\n') %>% cat
        # Please download the files from direct homepage link.
        #tb = try(download.file(url,destfile=name))
        #if('try-error' %in% class(tb)) stop('The file url seems to be changed. Please check the GTEx homepasge: https://gtexportal.org/home/')

        if(i==1) {
            untar(name,exdir=out)
            path = paste0(out,'/',f_target)
            f_li = list.files(path)
            f_name = paste0(out,'/gtex_files.txt')
            write.table(f_li,f_name,row.names=F,col.names=F,quote=F,sep='\t')
            paste0('    File write: ',f_name,'\n') %>% cat
        }
    }

    # Parsing GTEx file names
    cat("\n  Loading GTEx BED files\n")
    path   = read.delim(paste0(out,'/gtex_files.txt'),header=F) %>% unlist
    name   = path %>% file_path_sans_ext %>% file_path_sans_ext # exclude '.txt.gz'
    surfix = sub('.*.v8.','',name)
    f_df   = data.frame(path=paste0(out,'/',f_target,'/',path),surfix)
    f_sig  = subset(f_df,surfix=='signif_variant_gene_pairs')$path
    #f_sig  = subset(f.df,surfix=='egenes')$path # 1:1 match of eGene and the most significant SNP
    f_sig  = unlist(f_sig) %>% as.character

    # Load BED files
    gte_li=list(); n=length(f_sig)
    pb = winProgressBar(title="Loop progress",
        label="Ready to read files..",min=0,max=n,width=500)
    cat('  File reading...\n')
    for(i in 1:n) {
        tissue = sub('.*/(.*).v8..*','\\1',f_sig[i])
        if(file.exists(f_sig[i])) {
            #gte_li[[i]] = read.delim(gzfile(f_sig[i]),header=T)
            tb1 = read.delim(gzfile(f_sig[i]),header=T)
            tb2 = data.frame(tb1,tissue=rep(tissue,nrow(tb1)))
            gte_li[[i]] = tb2
        } else { cat(paste0('No such file: ',f_sig[i]))	}
        cat(paste0('    (',i,'/',n,') ',tissue,'\n'))
        ## Progress time ##
        setWinProgressBar(pb,i,label=paste0(round(i/n*100,0),
                " % (",i,"/",n,") done for ",pdtime(t0,2)))
        ###################
    }
    close(pb)
    gte_df = data.table::rbindlist(gte_li)

    # Summary print
    data.frame(gte_df$pval_nominal) %>% summary %>% print
    cat(paste0('  GTEx table, rows= ',dim(gte_df)[1],' cols= ',dim(gte_df)[2],'\n'))
    cat(paste0('  GTEx data read complete. ',pdtime(t0,2),'\n'))

    # Load annotation
    cat(paste0('  Loading annotation file\n'))
    ann = read.delim(gzfile(url_li[[2]][2]),header=T)
    cat(paste0('  Annotation file read complete. ',pdtime(t0,2),'\n'))
    cat(paste0('  Annotation file, rows= ',dim(ann)[1],' cols= ',dim(ann)[2],'\n'))
    gte_ann = merge(gte_df[,c(1:2,7:9,13)],ann[,c(1:3,7)],by='variant_id')
    cat(paste0('  GTEx annotation, rows= ',dim(gte_ann)[1],' cols= ',dim(gte_ann)[2],'\n'))

    # Save the compiled table as an RDS file
    cat(paste0('  Saving a compiled RDS file..'))
    f_name = paste0(out,'/Gtex_Analysis_v8_eQTL_rsid.rds')
    saveRDS(gte_ann,f_name)
    cat(paste0('  ',f_name,'\n'))
}

lncrna_down = function(
    out = 'db' # Download folder path 
) {
    # Function specific library
    source('src/saveasrds.r')

    # Download from RegulomeDB
    paste0('\n** Run function: db_download.r/regulome_down...\n') %>% cat
    mkdir_out(out)

    url_li = list(
        snplist = c(url='http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/snps_mod.txt',
                path=paste0(out,'/lncRNASNP2_snplist.txt')),
        lncrna  = c(url='http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncrnas.txt',
                path=paste0(out,'/lncrnas.txt')),
        disease = c(url='http://bioinfo.life.hust.edu.cn/static/lncRNASNP2/downloads/lncRNA_associated_disease_experiment.txt',
                path=paste0(out,'/lncrna-diseases_experiment.txt'))
    )
    k1=lapply(url_li,function(url) {
        tb = try(download.file(url[1],destfile=url[2]))
        if("try-error"%in% class(tb)) stop('The file url seems to be changed. Please check the lncRNASNP2 homepasge.')
    })

    # Convert to RDS files
    k2=lapply(url_li,function(url) {
        saveasrds(url[2])
    })
}

mkdir_out = function(dir) {
    if(file.exists(dir)) {
        paste0('  Directory exists: ',dir,'\n') %>% cat
    } else {
        dir.create(file.path(dir))
        paste0('  Directory generated: ',dir,'\n') %>% cat
    }
}

db_download = function(
    args = NULL
) {
    if(length(args$help)>0){ help = args$help
    } else                   help = FALSE
    if(help) {               cat(help_message); quit() }

    if(length(args$out)>0)   out  = args$out
    if(length(args$hg)>0)    hg   = args$hg

    source('src/pdtime.r'); t0=Sys.time()
    if(args$dbdown == 'roadmap') {
        roadmap_down(out)
    } else if(args$dbdown == 'encode') {
        encode_down(out)
    } else if(args$dbdown == 'regulome') {
        regulome_down(out)
    } else if(args$dbdown == 'gtex') {
        gtex_down(out)
    } else if(args$dbdown == 'lncrna') {
        lncrna_down(out)
    } else if(args$dbdown == 'gene') {
        biomart_gene(out,hg)
    } else {
        paste0('[Error] There is no such function in gwas_ldlink: ',
            paste0(args$dbdown,collapse=', '),'\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}