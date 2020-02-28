help_message = '
db_filter, v2020-02-27
This is a function call for filtering data.

Usage: Rscript postgwas-exe.r --dbfilt <function> --base <base file(s)> --out <out folder> <...>

Functions:
    roadmap  Filtering Roadmap data

Global arguments:
    --base   <base file/folder>
             Base file/folder path is mendatory.
    --out    <out folder>
             Out folder path is mendatory. Default is "db" folder.

Required arguments:
    --ctype  <cell type id>
             An optional argument for "roadmap" function.
             See cell-type number information at https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types.
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
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
    f_path = NULL, # Download folder path
    out    = 'db', # Out folder path
    ctype  = NULL, # Optional: Cell type ID 
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/roadmap_filt...\n') %>% cat
    if(length(ctype)>0) {
        cid    = numbers::formatC(ctype,width=3,flag='0') %>% as.character # Convert number to '###' format
        f_name = paste0(out,'/roadmap_',ctype,'_enh.bed')
    } else {
        cid    = numbers::formatC(c(1:129),width=3,flag='0') %>% as.character # 001~129
        f_name = paste0(out,'/roadmap_enh.bed')
    }

    # Read Roadmap files
    paste0('  Reading files..\n') %>% cat
    n=length(cid)
    road_li = lapply(c(1:n), function(i) {
        if(mod(i,10)==0) paste0('    ',i,'/',n,' being processed.\n') %>% cat
        #path = paste0(f_path,'E',cid[i],'_25_imputed12marks_hg38lift_dense.bed.rds')
        path = paste0(f_path,'/E',cid[i],'_25_imputed12marks_dense.bed.rds')
        road = try(readRDS(path))
        if('try-error' %in% class(road)) { # if file is not found.
            paste0('  ',path,' - file not found.\n') %>% cat
            road_enh = NULL
        } else { # If file exist.
            road_enh = subset(road, name %in% c("13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc"))
        }
        if(n<10) {
            paste0('  ',path,' ') %>% cat
            dim(road_enh) %>% print
        }
        return(road_enh)
    })
    paste0('  Finished reading and filtering ',n,' files.\n\n') %>% cat
    road_enh = data.table::rbindlist(road_li)
    rm(road_li)

    # Save file
    road_enh = road_enh[,c(1:3,6)]
    write.table(road_enh,f_name,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('Write file: ',f_name,'\n') %>% cat
}

db_filter = function(
    args = NULL
) {
    if(length(args$help)>0) {  help    = args$help
    } else                     help    = FALSE
    if(help)                   cat(help_message)

    if(length(args$base)>0)    b_path  = args$base
    if(length(args$out)>0)     out     = args$out
    if(length(args$debug)>0) { debug   = args$debug
    } else                     debug   = FALSE

    if(length(args$ctype)>0) { ctype   = args$ctype
    } else                     ctype   = NULL
    if(length(args$pval)>0)    pval    = args$pval

    source('src/pdtime.r'); t0=Sys.time()
    if(args$dbfilt == 'roadmap') {
        roadmap_filt(b_path,out,ctype,debug)
    } else if(args$dbfilt == 'gtex') {
        gtex_filt(b_path,out,pval,debug)
    } else {
        paste0('[Error] There is no such function in db_filter: ',
            paste0(args$ldlink,collapse=', '),'\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}