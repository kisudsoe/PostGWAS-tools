# This file was edited for postgwastools at 2020-02-27 by Seungsoo Kim

# Read file and save as rds file
suppressMessages(library(data.table))
suppressMessages(library(tools))
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))

saveasrds = function(f_paths) {
    n = length(f_paths); df.li = list()
    df.li=lapply(c(1:n),function(i) {
        f_path = f_paths[i]
        if(file_ext(f_path)=='gz') {         # for *.gz file
            try(R.utils::gunzip(f_path))
            f_path_ = file_path_sans_ext(f_path)
            return(fread(f_path_,sep='\t',header=T,stringsAsFactors=F))
        } else if(file_ext(f_path)=='csv') { # for *.csv file
            return(fread(f_path,sep=',',header=T,stringsAsFactors=F))
        } else if(file_ext(f_path) %in% c('tsv','txt')) {
            return(fread(f_path,sep='\t',header=T,stringsAsFactors=F))
        } else {
            print(paste0('[Error] File type ',file_ext(f_paths[i],' is not support yet.')))
        }
    })
    df = data.table::rbindlist(df.li)
    dim(df) %>% print
    f_name = paste0(f_paths[1],'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('File write: ',f_name,'\n'))
}

bedasrds = function(f_path) {
    df = as.data.frame(import(f_path,format='bed'))
    f_name = paste0(f_path,'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('File write: ',f_name,'\n'))
}