## Change log ##
# 1/30/2021 Written by Seungsoo Kim, PhD

## Parsing Arguments ##
suppressMessages(library(argparser))
p = arg_parser("Function to generate pivot table of Roadmap ChromHMM status")

### Example
p = add_argument(p, '--example',flag=T,
    help="See command examples by functions.")
example_msg = '
Rscript src/roadmap_summary.r --pivot \
    --f_roadmap_dist data/roadmap_dist

Rscript src/roadmap_summary.r --heatmap \
    --f_roadmap_summ data/summary_roadmap.tsv \
    --f_meta db/roadmap_meta.tsv \
    --annot BLOOD,PANCREAS,THYMUS \
    --f_snps data/gwas_5e-08_129_hg19.bed,data/snp_364_encode_dist.bed,data/snp_627_roadmap_dist.bed,data/snp_94_regulome2b.bed,data/snp_745_gtex.bed \
    --file_ext png

Rscript src/roadmap_summary.r --heatmap \
    --f_roadmap_summ data/summary_roadmap.tsv \
    --f_snp_filt data/snp_627_roadmap_dist.bed \
    --f_meta db/roadmap_meta.tsv \
    --annot BLOOD,PANCREAS,THYMUS \
    --f_snps data/gwas_5e-08_129_hg19.bed,data/snp_364_encode_dist.bed,data/snp_627_roadmap_dist.bed,data/snp_94_regulome2b.bed,data/snp_745_gtex.bed \
    --file_ext png
'

### Arguments for this function
p = add_argument(p, '--pivot',flag=T,
    help="[Function] Run pivot function to generate ChromHMM summary table.")
p = add_argument(p, '--f_roadmap_dist',
    help="[Path] Directory including TSV files of bedtools closest result of Roadmap ChromHMM data.")
p = add_argument(p, '--heatmap',flag=T,
    help="[Function] Run heatmap function to draw heatmap for ChromHMM summary.")
p = add_argument(p, '--f_roadmap_summ',
    help="[Path] Roadmap ChromHMM summary TSV file")
p = add_argument(p, '--f_snp_filt',
    help="[Path] SNP BED file for filtering ChromHMM summary table.")
p = add_argument(p, '--f_meta',
    help="[Path] Roadmap meta TSV file")
p = add_argument(p, '--annot',
    help="[BLOOD,PANCREAS,THYMUS,...] Choose ANATOMYs of roadmap meta to annotate heatmap.")
p = add_argument(p, '--f_snps',
    help="[Path,Path,...] SNPs BED files for annotation.")
p = add_argument(p, '--file_ext',
    help="[png/svg] Choose output figure format.")
p = add_argument(p, '--out',
    help="[Path] Out folder.")

argv = parse_args(p)


## Load Libraries ##
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))


## Function start ##
pivot = function(
    f_roadmap_dist = NULL,
    out = NULL
) {
    paste0('\n** Run pivot function in roadmap_summary.r **\n\n') %>% cat
    # Get file list from a directory
    f_roadmap_dists = list.files(f_roadmap_dist, full.names=T)
    n = length(f_roadmap_dists)
    paste0('* ',n,' files are found: Read [') %>% cat
    if(n==1) f_roadmap_dists = f_roadmap_dist

    # Read files
    m = n%/%10
    chromhmm_li = lapply(c(1:n), function(i) {
        if(i%%m==0) '.' %>% cat
        dist_df = read.delim(f_roadmap_dists[i],header=F,stringsAsFactors=F)
        k = ncol(dist_df)
        dist_df_col_sub = dist_df[,c(4,8,k)] %>% unique
        colnames(dist_df_col_sub) = c('Rsid','ChromHMM','Dist')
        dist_df_row_sub = subset(dist_df_col_sub,Dist==0)

        f_base = basename(f_roadmap_dists[i])
        #file_base = tools::file_path_sans_ext(f_base)
        file_base = strsplit(f_base,'\\_')[[1]][2]
        data.frame(dist_df_row_sub,Cell=file_base)
    })
    paste0('] = ') %>% cat

    # Merge data by Rsid
    chromhmm_df = data.table::rbindlist(chromhmm_li) %>% unique
    dim(chromhmm_df) %>% print

    # Pivot table
    paste0('* Transform table ') %>% cat
    #chromhmm_pivot = chromhmm_df %>%
    #    select(Rsid,ChromHMM,Cell) %>%
    #    spread(key=Cell,value=ChromHMM)
    chromhmm_pivot = chromhmm_df %>%
        select(Rsid,ChromHMM,Cell) %>%
        dcast(Rsid~Cell,value.var="ChromHMM")
    paste0('-> done\n') %>% cat

    # Write merged table as TSV file
    f_name1 = paste0(out,'/summary_roadmap.tsv')
    write.table(chromhmm_pivot,f_name1,row.names=F,sep='\t')
    paste0('* File write: ',f_name1,'\n') %>% cat

    # Mask data for enhancers
    paste0('\n* Filter enhancers ') %>% cat
    enhs = c("13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc")
    #`%notin%` = Negate(`%in%`) # define %notin% operator
    #chromhmm_pivot[chromhmm_pivot %notin% enhs] = NULL
    paste0('[') %>% cat
    k = nrow(chromhmm_pivot); l = k %/% 5
    enh_annot_li = lapply(c(1:k),function(i) {
        if(i%%l==0) '.' %>% cat
        x = chromhmm_pivot[i,] %>% unlist
        x_enhs = x[which(x %in% enhs)]
        if(length(x_enhs)>0) {
            eids = paste0(names(x_enhs),collapse=',')
        } else eids = '-'
        data.frame(Rsid=x[1],Enh_eid=eids,Enh_eid_num=length(x_enhs))
    })
    paste0('] = ') %>% cat
    enh_annot_df = data.table::rbindlist(enh_annot_li)
    dim(enh_annot_df) %>% print

    f_name2 = paste0(out,'/summary_roadmap_enh.tsv')
    write.table(enh_annot_df,f_name2,row.names=F,sep='\t')
    paste0('* File write: ',f_name2,'\n') %>% cat
}


heatmap = function(
    f_roadmap_summ = NULL,
    f_snp_filt     = NULL,
    f_meta         = NULL,
    f_snps         = NULL,
    annot          = NULL,
    file_ext       = NULL,
    out            = NULL
) {
    paste0('\n** Run heatmap function in roadmap_summary.r **\n\n') %>% cat
    suppressMessages(library(ComplexHeatmap))

    # Read file
    paste0('* Read file: ',f_roadmap_summ,' = ') %>% cat
    summ = read.delim(f_roadmap_summ,stringsAsFactors=F)
    dim(summ) %>% print

    if(!is.na(f_snp_filt)) {
        paste0('  [Option] Filt SNPs by ',f_snp_filt,' = ') %>% cat
        snp_filt = read.delim(f_snp_filt,header=F,stringsAsFactors=F)
        rsids = snp_filt[,4]
        summ = subset(summ,Rsid %in% rsids)
        dim(summ) %>% print
        
        # Set file name
        f_base = basename(f_snp_filt)
        file_base = tools::file_path_sans_ext(f_base)
        f_name = paste0('fig/roadmap_summ_',file_base,'.',file_ext)
    } else {
        f_name = paste0('fig/roadmap_summ','.',file_ext)
    }

    # Read meta-info. file
    paste0('* Meta-info table = ') %>% cat
    meta = read.delim(f_meta,stringsAsFactors=F)
    dim(meta) %>% print

    # Reorder meta-info following summ table EID order
    paste0('  Reordering ') %>% cat
    meta_name = tools::file_path_sans_ext(f_meta %>% basename)
    summ_col_names = colnames(summ)[-1]
    meta = meta[match(summ_col_names,meta$EID),]
    
    paste0('-> extract ha1 ') %>% cat
    summ_col = paste0(meta$EDACC_NAME," (",meta$EID,")")
    ha1 = HeatmapAnnotation(Name=anno_text(summ_col))
    show_column_names = FALSE

    paste0('-> extract ha2 ') %>% cat
    anatomy = meta$ANATOMY
    annots  = strsplit(annot,',')[[1]]
    `%notin%` = Negate(`%in%`) # define %notin% operator
    anatomy[meta$ANATOMY %notin% annots] = 'Other'
    anatomy_uq = anatomy %>% unique %>% sort
    ana_num = length(anatomy_uq)

    if(ana_num<5) { ann_cols = terrain.colors(length(anatomy_uq))
    } else ann_cols = rainbow(length(anatomy_uq))

    # Set column annotation color
    names(ann_cols) = anatomy_uq
    ha2 = HeatmapAnnotation(
        Anatomy=anatomy,
        col=list(Anatomy=ann_cols),
        gp=gpar(col="black",lwd=.5)
    )
    paste0('-> done\n') %>% cat

    # Configuration for heatmap
    paste0('  Configuration for heatmap ') %>% cat
    summ_mat = summ[,-1]
    rownames(summ_mat) = summ[,1]
    summ_mat[summ_mat=="1_TssA"]      = "1"
    summ_mat[summ_mat=="2_PromU"]     = "2"
    summ_mat[summ_mat=="3_PromD1"]    = "3"
    summ_mat[summ_mat=="4_PromD2"]    = "4"
    summ_mat[summ_mat=="5_Tx5'"]      = "5"
    summ_mat[summ_mat=="6_Tx"]        = "6"
    summ_mat[summ_mat=="7_Tx3'"]      = "7"
    summ_mat[summ_mat=="8_TxWk"]      = "8"
    summ_mat[summ_mat=="9_TxReg"]     = "9"
    summ_mat[summ_mat=="10_TxEnh5'"]  = "10"
    summ_mat[summ_mat=="11_TxEnh3'"]  = "11"
    summ_mat[summ_mat=="12_TxEnhW"]   = "12"
    summ_mat[summ_mat=="13_EnhA1"]    = "13"
    summ_mat[summ_mat=="14_EnhA2"]    = "14"
    summ_mat[summ_mat=="15_EnhAF"]    = "15"
    summ_mat[summ_mat=="16_EnhW1"]    = "16"
    summ_mat[summ_mat=="17_EnhW2"]    = "17"
    summ_mat[summ_mat=="18_EnhAc"]    = "18"
    summ_mat[summ_mat=="19_DNase"]    = "19"
    summ_mat[summ_mat=="20_ZNF/Rpts"] = "20"
    summ_mat[summ_mat=="21_Het"]      = "21"
    summ_mat[summ_mat=="22_PromP"]    = "22"
    summ_mat[summ_mat=="23_PromBiv"]  = "23"
    summ_mat[summ_mat=="24_ReprPC"]   = "24"
    summ_mat[summ_mat=="25_Quies"]    = "25"
    summ_mat = data.frame(apply(summ_mat,2,function(x) as.numeric(as.character(x))))
    my_col = c(
        "#FF0000","#FF4500","#FF4500","#FF4500","#008000","#008000","#008000","#009600","#C2E105","#C2E105","#C2E105","#C2E105","#FFC34D","#FFC34D","#FFC34D","#FFFF00","#FFFF00","#FFFF00","#FFFF66","#66CDAA","#8A91D0","#E6B8B7","#7030A0","#808080","#FFFFFF"
    )
    names(my_col) = 1:25
    '-> done\n' %>% cat

    # Read snps files
    paste0('* Read SNPs BED files ')  %>% cat
    f_snps_multi = strsplit(f_snps,'\\,')[[1]]
    n = length(f_snps_multi)
    ' [' %>% cat
    snps_li = lapply(c(1:n),function(i) {
        '.' %>% cat
        f_base = basename(f_snps_multi[i])
        file_base = tools::file_path_sans_ext(f_base)
        file_split = strsplit(file_base,'\\_')[[1]][1:3]
        snps_df = read.delim(f_snps_multi[i],header=F,stringsAsFactors=F)
        annot = paste0(file_split,collapse='_')
        Res = data.frame(Rsid=snps_df[,4],annot)
        colnames(Res) = c('Rsid',annot)
        return(Res)
    })
    paste0('] Annotation merge ') %>% cat
    snps_df = Reduce(function(x,y) merge(x=x,y=y,by='Rsid',all=T), snps_li)
    snps_merge = merge(data.frame(Rsid=summ[,1]),snps_df,by='Rsid',all=T) %>% unique
    if(!is.na(f_snp_filt)) {
        paste0('-> [Option] Filt by ',length(rsids),' snps ') %>% cat
        snps_merge = subset(snps_merge,Rsid %in% rsids)
    }
    rownames(snps_merge) = snps_merge[,1]
    snps_mat = as.matrix(snps_merge[,c(-1)])
    paste0('-> done\n') %>% cat

    # Preparing heatmap
    paste0('* Draw heatmap ') %>% cat
    m = nrow(summ_mat)
    if(m>200) {
        wh=c(22,20); show_row_names = FALSE
    } else { 
        wh=c(22,m*0.11+9); show_row_names = TRUE
    }
    if(file_ext=='png') { png(f_name, width=wh[1], height=wh[2], units='in', res=300)
    } else if(file_ext=='svg') { svg(f_name, width=wh[1], height=wh[2]) }

    summ_mat = as.matrix(summ_mat)
    ht1=Heatmap(
        summ_mat,
        name = 'Chr.\nstatus',
        col  = my_col,
        #rect_gp = gpar(col = "black", lwd = .5),
        column_dend_height = unit(1, "in"),
        row_dend_width     = unit(1, "in"),
        column_names_max_height = unit(10,"in"),
        show_column_names  = show_column_names,
        show_row_names     = FALSE,
        bottom_annotation  = ha1,
        top_annotation     = ha2,
        cluster_rows       = TRUE,
        cluster_columns    = TRUE
    )

    k = ncol(snps_mat)+1
    mycol2 = structure(2:k, names=colnames(snps_mat))
    ht2=Heatmap(
        snps_mat,
        name = 'SNPs',
        col  = mycol2,
        show_row_names = show_row_names
    )
    print(ht1+ht2)
    dev.off()
    paste0('-> done\n') %>% cat
    paste0('\nSave as ',f_name,'\n') %>% cat
}

## Function end ##


## Run function ##
source('src/pdtime.r')
t0 = Sys.time() # timer

if(argv$example) {
    cat(example_msg)
} else if(argv$pivot) {
    pivot(
        f_roadmap_dist = argv$f_roadmap_dist,
        out            = argv$out
    )
} else if(argv$heatmap) {
    heatmap(
        f_roadmap_summ = argv$f_roadmap_summ,
        f_snp_filt     = argv$f_snp_filt,
        f_meta         = argv$f_meta,
        annot          = argv$annot,
        f_snps         = argv$f_snps,
        file_ext       = argv$file_ext,
        out            = argv$out
    )
}

paste0('\n',pdtime(t0,1),'\n') %>% cat