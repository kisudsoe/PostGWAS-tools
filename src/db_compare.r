help_message = '
gwas_catalog, v2020-04-25
This is a function for DB comparison with PCA analysis.

Usage: Rscript postgwas-exe.r --dbcomp <functions> --base <base file> --out <out folder> <...>

Functions:
    bedmerge   Generating merged BED file for bedtools merge input.
    uniontb    Generating union table from the dist files.
    pcaplot    Performing PCA analysis and draw plot.

Global arguments:
    --base     Base folder path for input BED files are mendatory.
    --out      <Default: data>
               Out files tartget path is mendatory. Default is "bedmerge.bed" file.

Required arguments:
    --meta     A mendatory argument for the "pcaplot" function.
    --filt     An optional argument for the "pcaplot" function.
               Filtering by group
    --hic      <Default: FALSE>
               An optional argument for the "uniontb" function.
               If input files are hic data, please set the --hic TRUE.
'

## Load libraries ##
suppressMessages(library(dplyr))

## Functions Satrt ##
pcaplot = function(
    f_uniontb = NULL,   # Input uniontb file path
    out       = 'data', # Out folder path
    meta      = NULL,   # Meta file for uniontb file columns
    filt      = NULL,   # Option for uniontb filter by group
    debug
) {
    # Preparing...
    suppressMessages(library(FactoMineR))
    suppressMessages(library(ggplot2))
    suppressMessages(library(plyr))
    paste0('\n** Run function: db_compare.r/pcaplot... \n') %>% cat

    # Read uniontb file
    paste0('  Read file: ',f_uniontb,'... ') %>% cat
    nm = tools::file_path_sans_ext(basename(f_uniontb))
    uniontb = read.csv(f_uniontb,stringsAsFactors=F)
    dim(uniontb) %>% print

    # Convert TRUE/FALSE to 1/0
    paste0('  Convert T/F to 1/0... ') %>% cat
    uniontb = uniontb[order(uniontb$ID),] # reorder row by coord
    uniontb_val = uniontb[,-1]
    cols = sapply(uniontb_val, is.logical)
    uniontb_val[,cols] = lapply(uniontb_val[,cols],as.numeric)
    'done\n' %>% cat

    # Read meta file
    paste0('  Read meta file: ',meta,'... ') %>% cat
    meta = read.delim(meta,stringsAsFactors=F)
    dim(meta) %>% print

    # Filter uniontb column (option)
    if(!is.null(filt)) {
        paste0('  Filter option: ') %>% cat
        paste0(filt,collapse=', ') %>% cat; cat(' >> ')
        colnm = colnames(uniontb_val)
        which_col = which(meta$group %in% filt)
        length(which_col) %>% print
        
        uniontb_val = uniontb_val[,which_col]
        meta = meta[which_col,]
    }
    
    # PCA analysis
    paste0('  PCA analysis... ') %>% cat
    pca = PCA(t(uniontb_val)) # FactoMineR library
    'done\n' %>% cat

    # Save as CSV file
    pca_coord = pca$ind$coord %>% as.data.frame
    pca_coord = cbind(
        pca_coord,
        group = meta$group,
        cell_type = meta$cell_type
    )
    f_name1 = paste0(out,'/',nm,'_pca.csv')
    write.csv(pca_coord,f_name1,quote=F)
    paste0('  Write file: ',f_name1,'\n') %>% cat

    # Set hull
    paste0('  Draw PCA plot... ') %>% cat
    a_hull = function(pca_coord) {
        pca_coord[chull(
            pca_coord$Dim.1,
            pca_coord$Dim.2),]
    }
    pca_hull = ddply(pca_coord,"group",a_hull)
    paste0('hull.. ') %>% cat
    
    # Draw PCA plot
    f_name2 = paste0(out,'/',nm,'_pca.png')
    point_shapes = c(15,16,17,18,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
    n = length(pca_coord$cell_type) %>% unique
    ggplot()+theme_bw()+
        geom_point(pca_coord,mapping=aes(Dim.1,Dim.2,colour=group,shape=cell_type),size=3)+
        geom_polygon(pca_hull,mapping=aes(Dim.1,Dim.2,fill=group),alpha=0.3)+
        scale_shape_manual(values=point_shapes[1:n])+
        labs(title=paste0('PCA for ',nm),
            x=paste0('Dim 1 (',round(pca$eig[1,2],2),'%)'),
            y=paste0('Dim 2 (',round(pca$eig[2,2],2),'%)')
        )
    ggsave(f_name2,width=7,height=5)
    'done\n' %>% cat
    paste0('Save plot: ',f_name2,'\n') %>% cat
}

uniontb = function(
    f_dists = NULL,  # Input dist file dir path
    out     = 'uniontb.csv', # Out file path
    hic     = FALSE, # Whether hic data
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_compare.r/uniontb... \n') %>% cat
    if(hic=='TRUE') hic = TRUE

    # Read files paths
    n = length(f_dists); paths = NULL
    for(i in 1:n) {
        paste0(i,' ',f_dists[i],': ') %>% cat
        paths_i = list.files(f_dists[i],full.names=T)
        length(paths_i) %>% print
        if(length(paths_i)<20) { paste0('  ',paths_i,'\n') %>% cat }
        paths = c(paths,paths_i)
    }
    path_n = length(paths)
    paste0('** ',path_n,' dist file paths are ready.\n') %>% cat

    # Extract names
    nms = tools::file_path_sans_ext(basename(paths))

    # Read data as list
    paste0('\n  Read data... ') %>% cat
    dist_li = lapply(c(1:path_n),function(i) {
        dist = read.delim(paths[i],header=F,)
        k = ncol(dist)
        p = k-3; q = k-2; r = k-1

        if(hic) {
            loops_full  = lapply(dist[,4],function(x)
                strsplit(x%>%as.character,"\\.")[[1]][1]) %>% unlist
            loops  = loops_full %>% unique
            loop_n = length(loops)

            ids_li = lapply(c(1:loop_n),function(i) {
                which_i = which(loops_full==loops[i])
                if(length(which_i)==2) {
                    loop = dist[which_i,]
                    id1 = paste0(loop[1,p],":",loop[1,q],"-",loop[1,r])
                    id2 = paste0(loop[2,p],":",loop[2,q],"-",loop[2,r])
                    return(paste0(id1," - ",id2))
                }
            })
            ids = ids_li %>% unlist %>% unique
        } else {
            ids = paste0(dist[,p],":",dist[,q],"-",dist[,r])
        }
        
        if(path_n>100) {
            if(i%%20==0) { paste0('\n    ',i,'/',path_n,' process done. ') %>% cat }
        } else if(path_n>20) {
            if(i%%10==0) { paste0('\n    ',i,'/',path_n,' process done. ') %>% cat }
        } else if(path_n<=20) {
            if(i%%5==0) { paste0('\n    ',i,'/',path_n,' process done. ') %>% cat }
        }
        return(ids)
    })
    paste0(length(dist_li),'.. ') %>% cat

    # Generate union table
    paste0('\n  Data union... ') %>% cat
    posid = dist_li %>% unlist %>% unique
    m = length(dist_li)
    df = NULL
    for(i in 1:m) {
        df = cbind(df,dist_li[[i]][match(posid,dist_li[[i]])])
    }
    df1 = (df!='')           # Transcform values to TRUE, if ID exists
    df1[is.na(df)] = FALSE   # Transform NA to FALSE value
    df1 = as.data.frame(df1) # Make 'df1' to data.frame
    colnames(df1) = nms

    df2 = data.frame(
        ID = posid,
        df1
    )
    dim(df2) %>% print

    # Save as CSV file
    write.csv(df2,out,row.names=F)
    paste0('  Write file: ',out,'\n') %>% cat
}

bedmerge = function(
    f_beds = NULL, # Input BED file dir path
    out    = 'bedmerge.bed', # Out file path
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_compare.r/bedmerge... \n') %>% cat
    #ifelse(!dir.exists(out), dir.create(out),''); 'ready\n' %>% cat

    # Read files paths
    n = length(f_beds); paths = NULL
    for(i in 1:n) {
        paste0(i,' ',f_beds[i],': ') %>% cat
        paths_i = list.files(f_beds[i],full.names=T)
        length(paths_i) %>% print
        if(length(paths_i)<10) { paste0('  ',paths_i,'\n') %>% cat }
        paths = c(paths,paths_i)
    }
    path_n = length(paths)
    paste0('** ',path_n,' BED file paths are ready.\n') %>% cat

    paste0('\n  Data merge... ') %>% cat
    beds_li = lapply(c(1:path_n),function(i) {
        #paste0('  ',i,' ',paths[i]) %>% cat
        bed = try(read.delim(paths[i],header=F))
        if('try-error' %in% class(bed)) {
            paste0('  [ERROR] ',paths[i],'\n') %>% cat
            return(NULL)
        }
        if(path_n>20) {
            if(i%%20==0) { paste0('\n    ',i,'/',path_n,' process done. ') %>% cat }
        }
        return(bed[,1:4]) # chr, start, end, name
    })
    paste0(length(beds_li),'.. ') %>% cat
    beds = data.table::rbindlist(beds_li)
    paste0('  done: ') %>% cat
    dim(beds) %>% print

    # Write bed file
    write.table(beds,out,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('  Write file: ',out,'\n') %>% cat
}
## Functions End ##

db_compare = function(
    args = NULL
) {
    if(length(args$help)>0) {     help = args$help
    } else                        help = FALSE
    if(help) { cat(help_message); quite() }

    # Global arguments
    if(length(args$base)>0)       b_path   = args$base
    if(length(args$out)>0)        out      = args$out
    if(length(args$debug)>0) {    debug    = args$debug
    } else                        debug    = FALSE

    # Required arguments
    if(length(args$meta)>0) {     meta     = args$meta
    } else                        meta     = NULL
    if(length(args$filt)>0) {     filt     = args$filt
    } else                        filt     = NULL
    if(length(args$hic)>0) {      hic      = args$hic
    } else                        hic      = FALSE

    # Run function
    source('src/pdtime.r'); t0=Sys.time()
    if(args$dbcomp=='bedmerge') {
        bedmerge(b_path,out,debug)
    } else if(args$dbcomp=='uniontb') {
        uniontb(b_path,out,hic,debug)
    } else if(args$dbcomp=='pcaplot') {
        pcaplot(b_path,out,meta,filt,debug)
    } else {
        paste0('[Error] There is no such function "',args$dbfilt,'" in db_compare: ',
            paste0(args$ldlink,collapse=', '),'\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}