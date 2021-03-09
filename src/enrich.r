## Change log ##
# Written by Seungsoo Kim, PhD
# debug01 21.01.26

## Parsing Arguments ##
suppressMessages(library(argparser))
p = arg_parser("Function for enrichment analysis (Fisher/Permutation)")

### Examples
p = add_argument(p, '--example',flag=T,
    help="See command examples by functions.")
example_msg = '
Rscript src/enrich.r --splittfbs \
    --tfbs db/wgEncodeRegTfbsClusteredWithCellsV3.bed \
    --out db/tfbs_cell

Rscript src/enrich.r --splitgtex \
    --gtex db/gtex_signif_5e-8.tsv.rds \
    --out db/gtex_tsv

Rscript src/enrich.r --roadmap_perm \
    --gwas_snp data/seedSNP_1817_bm.bed \
    --f_roadmap db/roadmap_bed \
    --db_source roadmap_bed \
    --perm_n 1000 \
    --out enrich/roadmap_enh

Rscript src/enrich.r --roadmap_perm \
    --gwas_snp data/gwas_5e-08_129_hg19.bed \
    --f_roadmap db/encode_bed \
    --db_source encode_bed \
    --perm_n 100 \
    --out enrich/encode_tfbs

Rscript src/enrich.r --gtex_perm \
    --gwas_snp data/seedSNP_1817_bm.bed \
    --gtex_base db/gtex_signif_5e-8.db \
    --gtex_median_tpm db/gtex_analysis_v8_rnaseq_gene_median_tpm_ensgid.gct.rds
    --permn 1000 \
    --out enrich

Rscript src/enrich.r --heatmap \
    --pmdata enrich/roadmap_bed-snp_484_roadmap_dist-permn_100-zscore.tsv \
    --meta db/roadmap_meta.tsv \
    --out enrich \
    --range -3,3 \
    --annot BLOOD,PANCREAS,THYMUS \
    --fileext png

Rscript src/enrich.r --conv_dicevcf2db \
    --dice_vcf "db/Schmiedel-DICE" \
    --pval 5e-8 \
    --out db

Rscript src/enrich.r --gtex_perm \
	--gwas_snp data/seedSNP_1817_bm.bed \
    --gtex_base db/dice_eqtl_5e-08.db \
    --gtex_median_tpm db/dice-mean_tpm_merged.rds \
    --perm_n 1000 \
    --out enrich

Rscript src/enrich.r --eqtl_hyper hypergeometric \
    --gwas_snp data/seedSNP_1817_bm.bed \
    --eqtl_db db/dice_eqtl_5e-08.db \
    --out enrich

# Not working yet
Rscript src/enrich.r --eqtl_hyper kstest \
    --gwas_snp data/seedSNP_1817_bm.bed \
    --eqtl_db db/dice_eqtl_5e-08.db \
    --gene_median_tpm db/dice-mean_tpm_merged.rds \
    --out enrich
'

### Shared Arguments
p = add_argument(p,'--gwas_snp',
    help="[Path] GWAS SNP list BED file. Columns: <Chr> <Start> <End> <Rsid> <...>")
p = add_argument(p,'--out',
    help="[Path] Target directory for output files.")
p = add_argument(p,'--verbose',flag=T,
    help="[Option] Show detailed process.")
p = add_argument(p,'--perm_n',default=1000,type="numeric",
    help="[Number] Set permutation number. Default=1000")

### Arguments for perm_test function
p = add_argument(p,'--roadmap_perm',flag=T,
    help="[Function] Run permutation test for Roadmap/ENCODE data.")
p = add_argument(p,'--f_roadmap',
    help="[Path] Directory for Roadmap BED files including chromosome status for permtest functions.")
p = add_argument(p,'--db_source',
    help="[roadmap_bed/encode_bed] Type in one of these options for source of your chromosome status.")

### Arguments for draw_heatmap function
p = add_argument(p,'--heatmap',flag=T,
    help="[Function] Draw a heatmap from roadmap_perm results.")
p = add_argument(p,'--pm_data',
    help="[Path] Z-score table file.")
p = add_argument(p,'--meta',
    help="[Path] Add roadmap meta-info file for heatmap annotation.")
p = add_argument(p,'--range',default="-4,4",
    help="[-4,4] Set coloring Z-score range to display.")
p = add_argument(p,'--annot',
    help="[BLOOD,PANCREAS,THYMUS ...] Choose ANATOMYs of roadmap meta-info to annotate heatmap.")
p = add_argument(p,'--file_ext',
    help="[png/sgv] Choose output figure format.")

### Arguments for split_tfbs function
p = add_argument(p,'--split_tfbs',flag=T,
    help="[Function] Split ENCODE TFBS BED file by cell types.")
p = add_argument(p,'--tfbs',
    help="[Path] ENCODE TFBS BED file.")

### Arguments for split_gtex function
p = add_argument(p,'--split_gtex',flag=T,
    help="[Function] Split GTEx eQTL file by tissue types.")
p = add_argument(p,'--gtex',
    help="[Path] GTEx eQTL RDS file.")

### Arguments for gtex_perm_test function
p = add_argument(p,'--gtex_perm',flag=T,
    help="[Function] Run permutation test for GTEx eQTL data.")
p = add_argument(p,'--gtex_base',
    help="[Path] GTEx eQTL DB file.")
p = add_argument(p,'--gtex_median_tpm',
    help="[Path] GTEx gene median tpm GCT RDS file.")

### Arguments for gtex_pm_plot function
p = add_argument(p,'--gtex_pm_plot',flag=T,
    help="[Function] Draw plot for GTEx permutation test result.")
p = add_argument(p,'--fgsea_rds',
    help="[Path] RDS file for GTEx permutation test result.")

### Arguments for conv_Dicevcf2db function
p = add_argument(p,'--conv_dicevcf2db',flag=T,
    help="[Function] Convert the DICE eQTL VCF files to database for GSEA.")
p = add_argument(p,'--dice_vcf',
    help="[Path] DICE eQTL VCF file or parent directory.")
p = add_argument(p,'--pval',default="5e-8",
    help="[5e-8] DICE eQTL VCF file or parent directory.")

### Arguments for conv_Dicevcf2db function
p = add_argument(p,'--eqtl_hyper',
    help="[Function] Run eqtl_hyper function for enrichment test. Choose function between 'hypergeometric' or 'kstest'.")
p = add_argument(p,'--eqtl_db',
    help="[Path] eQTL DB file.")
p = add_argument(p,'--gene_median_tpm',
    help="[Path] Gene median tpm RDS file.")

argv = parse_args(p)
eqtl_hyper_n = length(argv$eqtl_hyper)


## Load Common Libraries ##
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))


## Functions ##
read_dicevcf = function(
    f_dicevcf = NULL,
    pval      = 5e-8
) {
    # Read file
    paste0('  ',f_dicevcf,' = ') %>% cat
    dicevcf = read.delim(f_dicevcf,comment.char="#",header=F,stringsAsFactors=F)
    colnames(dicevcf) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
    nrow(dicevcf) %>% cat

    # Split and compile the SNP & INFO data
    paste0('; compile INFO ') %>% cat
    dice_info_li = lapply(dicevcf$INFO, function(x) {
        info_split = strsplit(x,"\\;")[[1]]
        info_split_li = strsplit(info_split,"\\=")
        info_df = data.frame(
            info_split_li[[1]][2],
            info_split_li[[2]][2],
            info_split_li[[3]][2],
            info_split_li[[4]][2])
        colnames(info_df) = c(
            info_split_li[[1]][1],
            info_split_li[[2]][1],
            info_split_li[[3]][1],
            info_split_li[[4]][1])
        return(info_df)
    })
    dice_info_df = data.table::rbindlist(dice_info_li)
    paste0('-> preparing ') %>% cat
    file_base = basename(f_dicevcf)
    f_base = tools::file_path_sans_ext(file_base)
    dice_merge = data.frame(Rsid=dicevcf$ID,dice_info_df,Cell=f_base)
    paste0('-> filtering = ') %>% cat
    dice_sub = subset(dice_merge, Pvalue<pval)
    dim(dice_sub) %>% print

    return(dice_merge)
}


conv_dicevcf2db = function(
    f_dicevcf = NULL,
    pval      = 5e-8,
    out       = NULL
) {
    paste0('\n** Run conv_dicevcf2db function in enrich.r **\n\n') %>% cat
    suppressMessages(library(RSQLite))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file/directory
    dicevcf_files = list.files(f_dicevcf, pattern="\\.vcf$", full.names=T)
    f_num = length(dicevcf_files)
    if(f_num==0) {
        dicevcf_files = f_dicevcf
        #paste0('[Error] There is no file! Check the target directory again.\n') %>% cat
    }
    paste0('* ',f_num,' VCF files found in ',f_dicevcf,';\n') %>% cat
    dice_li = lapply(dicevcf_files, function(x) read_dicevcf(x, pval))
    dice_df = data.table::rbindlist(dice_li)
    paste0('* Read done = ') %>% cat; dim(dice_df) %>% print

    # Convert dice table to db format
    f_db = paste0(out,'/dice_eqtl_',pval,'.db')
    conn = dbConnect(RSQLite::SQLite(),f_db)
    dbWriteTable(conn, "dice_eqtl", dice_df)
    paste0('\nWrite db file: ',f_db,'\n') %>% cat
}


split_gtex = function(
    f_gtex = NULL,
    out    = NULL
) {
    paste0('\n** Run split_gtex function in enrich.r **\n\n') %>% cat
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    paste0('* GTEx table = ') %>% cat
    gtex = readRDS(f_gtex)
    dim(gtex) %>% print

    # Get unique tissue types
    tissue_types = gtex$tissue %>% unique %>% sort
    n = length(tissue_types)
    paste0('* ',n,' unique tissue types are found.\n\n') %>% cat

    # Split file by unique tissue types
    source('src/pdtime.r')
    o=lapply(c(1:n), function(i) {
        t1 = Sys.time()
        paste0(i,' ',tissue_types[i],':\t') %>% cat
        gtex_sub = subset(gtex,tissue==tissue_types[i])
        m = nrow(gtex_sub)
        paste0(m,' pairs, SNPs = ') %>% cat
        snps_len = gtex_sub$variant_id %>% unique %>% length
        gene_len = gtex_sub$gene_id %>% unique %>% length
        paste0(snps_len,', genes = ',gene_len,'. ') %>% cat

        # Save file as BED format
        f_name = paste0(out,'/',tissue_types[i],'.tsv')
        write.table(gtex_sub[,c(1:5,9)],f_name,sep='\t',row.names=F,quote=F)
        paste0('Save: ',f_name,'; ',pdtime(t1,2),'\n') %>% cat
        return(NULL)
    })
}


split_tfbs = function(
    f_tfbs = NULL,
    out    = NULL
) {
    paste0('\n** Run split_tfbs function in enrich.r **\n\n') %>% cat
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    paste0('* ENCODE TFBS table = ') %>% cat
    tfbs = read.delim(f_tfbs, stringsAsFactors=F)
    colnames(tfbs) = c('Chr','Start','End','TF','ids','Cells')
    dim(tfbs) %>% print

    # Get unique cell types
    cell_types = strsplit(tfbs$Cells,',') %>% unlist %>% unique %>% sort
    n = length(cell_types)
    paste0('* ',n,' unique cell types are found.\n\n') %>% cat

    # Expand data by unique cell types
    source('src/pdtime.r')
    o=lapply(c(1:n), function(i) {
        t1=Sys.time()

        # Subset ENCODE data by cell type
        paste0(i,' ',cell_types[i],':\t') %>% cat
        tfbs_sub = subset(tfbs, grepl(cell_types[i],tfbs$Cells))
        m = nrow(tfbs_sub)
        progress = m%/%10
        paste0(m,' regions, filtering [') %>% cat

        # Split file by unique cell types
        tfbs_filt_li = lapply(c(1:m), function(j) {
            if(j%%progress==0) { '.' %>% cat }

            tfbs_sub_j = tfbs_sub[j,1:4]
            tfbs_cell = tfbs_sub[j,]$Cells
            cells = strsplit(tfbs_cell,',')[[1]]
            which_cells = which(cells==cell_types[i])
            if(length(which_cells)>0) {
                tfbs_sub_df = data.frame(tfbs_sub_j, Cell=cells[which_cells]) # debug
            } else tfbs_sub_df = NULL
            return(tfbs_sub_df)
        })
        paste0('] row = ') %>% cat
        tfbs_filt_df = data.table::rbindlist(tfbs_filt_li)
        tf_len = tfbs_filt_df$TF %>% unique %>% length
        paste0(nrow(tfbs_filt_df),', TFs = ',tf_len,'. ') %>% cat

        # Save file as BED format
        f_name = paste0(out,'/',cell_types[i],'.bed')
        write.table(tfbs_filt_df[,1:4],f_name,sep='\t',row.names=F,col.names=F,quote=F)
        paste0('Save: ',f_name,'; ',pdtime(t1,2),'\n') %>% cat
        return(NULL)
    })
}


draw_heatmap = function(
    f_pmdata = NULL,
    f_meta   = NULL,
    range    = NULL,
    out      = 'enrich',
    annot    = NULL,
    fileext  = 'png'
) {
    paste0('\n** Run draw_heatmap function in enrich.r **\n\n') %>% cat
    suppressMessages(library(ComplexHeatmap))
    suppressMessages(library(circlize))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    f_pmdata_multi = strsplit(f_pmdata,'\\, ')[[1]]
    n = length(f_pmdata_multi)
    paste0('* Read ',n,' file(s):\n') %>% cat
    pmdata_li = lapply(c(1:n), function(i) {
        paste0(i,' ',basename(f_pmdata_multi[i])) %>% cat
        pmdata = read.delim(f_pmdata_multi[i], stringsAsFactors=F)

        # Extract row order
        paste0(';\t-> Reorder rows by Status ') %>% cat
        row_nums_li = lapply(pmdata[,1],function(x) {
            strsplit(x,'\\_')[[1]][1]
        })
        row_nums = unlist(row_nums_li) %>% as.numeric
        pmdata = pmdata[order(row_nums),]
        '-> done\n' %>% cat

        # Prepare heatmap table
        pmdata_mat = pmdata[,-1]
        rownames(pmdata_mat) = pmdata$Status
        return(pmdata_mat)
    })

    # [Optional] Read meta-info. file
    if(length(f_meta)>0) {
        paste0('\n* [Optional] Add meta-info. table = ') %>% cat
        meta = read.delim(f_meta, stringsAsFactors=F)
        dim(meta) %>% print

        # Reorder meta-info. with original colnames
        paste0('  Reorder meta-info ') %>% cat
        meta_db = tools::file_path_sans_ext(f_meta %>% basename)
        pmdata_colnms = colnames(pmdata_li[[1]])
        meta = meta[match(pmdata_colnms,meta$EID),]
        paste0('-> extract ha1') %>% cat
        if(meta_db=='roadmap_meta') {
            pmdata_col = paste0(meta$EDACC_NAME," (",meta$EID,")")
        } else if(meta_db=='encode_meta') {
            #pmdata_mat = pmdata_mat[,meta$Cell_Name] # Error: undefined columns selected
            pmdata_col = paste0(meta$Cell_Name," (",meta$Symbol,")")
        }
        ha1 = HeatmapAnnotation(Name=anno_text(pmdata_col))
        show_column_names = FALSE

        # [Optional] Add column annotation to heatmap by ANATOMY
        paste0('-> extract ha2 ') %>% cat
        if(length(f_meta)>0 & length(annot)>0) {
            if(meta_db=='roadmap_meta') { anatomy = meta$ANATOMY
            } else if(meta_db=='encode_meta') { anatomy = meta$Tissue }
            annots = strsplit(annot,',')[[1]]
            `%notin%` = Negate(`%in%`) # define %notin% operator
            anatomy[meta$ANATOMY %notin% annots] = 'Other' # A bug to change as NA
            anatomy_uq = anatomy %>% unique %>% sort
            ana_num = length(anatomy_uq)

            if(ana_num<5) { ann_cols = terrain.colors(length(anatomy_uq))
            } else ann_cols = rainbow(length(anatomy_uq))

            ## Set column annotation color
            names(ann_cols) = anatomy_uq
            ha2 = HeatmapAnnotation(
                Anatomy=anatomy,
                col=list(Anatomy=ann_cols),
                gp=gpar(col="black", lwd=.5)
            )
        } else { ha2 = NULL }
        paste0('-> done\n') %>% cat
    } else {
        paste0('* [Notice] Meta-info table is not found.') %>% cat
        show_column_names = TRUE
        ha1 = NULL; ha2 = NULL
    }
    
    # Configuration for heatmap
    paste0('* Configuration for heatmap ') %>% cat
    z_range = strsplit(range,',')[[1]] %>% as.numeric
    my_col = colorRamp2(c(z_range[1],0,z_range[2]),c("#2E86C1","#FEF9E7","#C0392B"))
    if(meta_db=='roadmap_meta') {
        cluster_row_col = c(FALSE,TRUE)
        wh = c(20,5*n+5)
    } else if(meta_db=='encode_meta') {
        cluster_row_col = c(FALSE,FALSE)
        wh = c(22,25)
    }

    # Prepare for draw heatmap
    paste0('-> set file name ') %>% cat
    if(n==1) {
        file_base = basename(f_pmdata)
        file_name = tools::file_path_sans_ext(file_base)
        f_name = paste0(out,'/',file_name,'.',fileext)
    } else {
        f_name = paste0(out,'/merged_',n,'.',fileext)
    }
    if(fileext=='png') { png(f_name,width=wh[1],height=wh[2],units='in',res=300)
    } else if(fileext=='svg') { svg(f_name,width=wh[1],height=wh[2]) }
    paste0('-> done\n') %>% cat

    paste0('* Draw heatmap: ') %>% cat
    ht_list = NULL
    for(i in 1:n) {
        paste0(i,' ') %>% cat
        pmdata_mat = as.matrix(pmdata_li[[i]])
        if(i>1) { ha2_ann = NULL } else ha2_ann = ha2
        if(i<n) { ha1_ann = NULL } else ha1_ann = ha1
        if(n>1) { 
            file_base = basename(f_pmdata_multi[[i]])
            file_name = tools::file_path_sans_ext(file_base)
            ha3_ann = rowAnnotation(SNPs=anno_block(#gp=gpar(fill=i+1),
                labels=paste0(file_name),
                labels_gp=gpar(col="black",fontsize=12)))
        } else ha3_ann = NULL
        ht=Heatmap(
            pmdata_mat,
            name = "Enrich.\nz-score",
            col  = my_col,
            rect_gp = gpar(col = "black", lwd = .5),
            column_dend_height = unit(1, "in"),
            row_dend_width = unit(1, "in"),
            column_names_max_height = unit(10,"in"),
            show_column_names = show_column_names,
            bottom_annotation = ha1_ann,
            top_annotation = ha2_ann,
            left_annotation = ha3_ann,
            cluster_rows = cluster_row_col[1],
            cluster_columns = cluster_row_col[2]
        )
        ht_list = ht_list %v% ht
    }
    print(ht_list)
    dev.off()
    paste0('-> done\n') %>% cat
    paste0('\nSave as ',f_name,'\n') %>% cat
}


read_files_enrich = function(
    f_gwas_snp      = NULL,
    eqtl_db         = NULL,
    gene_median_tpm = NULL,
    statistics      = NULL
) {
    suppressMessages(library(data.table))
    suppressMessages(library(tidyr))
    suppressMessages(library(RSQLite))

    paste0('* Input SNPs = ') %>% cat
    snps = read.delim(f_gwas_snp, header=F, stringsAsFactors=F)
    colnames(snps) = c('Chr','Start','End','Rsid')
    dim(snps) %>% print

    #paste0('* GTEx eQTL pairs = ') %>% cat
    #eqtl_pairs = readRDS(eqtl_db)
    #colnames(eqtl_pairs) = c('Variant_id','Gene_id','Pval','Slope','Slope_se','Tissue','Chr','Pos','Rsid')
    #eqtl_pairs$Ensgid = sapply(eqtl_pairs$Gene_id,function(x) strsplit(x,'\\.')[[1]][1])
    #dim(eqtl_pairs) %>% print

    paste0('* ',eqtl_db,', eQTL pairs = ') %>% cat
    conn = dbConnect(RSQLite::SQLite(),eqtl_db)
    db_name = basename(eqtl_db)
    if(db_name=='gtex_signif_5e-8.db') {
        eqtl_pairs = dbGetQuery(conn,"SELECT Ensgid, Tissue, Rsid FROM gtex") %>% unique
        eqtl_pairs$eqtl_pair = paste0(eqtl_pairs$Rsid,'-',eqtl_pairs$Ensgid)
        eqtl_snp_gene = eqtl_pairs$eqtl_pair %>% unlist %>% unique
        length(eqtl_snp_gene) %>% print

        if(!is.na(gene_median_tpm)) {
            paste0('* GTEx gene median tpm GCT = ') %>% cat
            gct = readRDS(gene_median_tpm)
            nrow(gct) %>% cat
            gct_gather = gct %>% gather("Tissue","TPM",-Ensgid,-Description)
        } else gct_gather = NULL
    } else if(db_name=='dice_eqtl_5e-08.db') {
        eqtl_pairs = dbGetQuery(conn,"SELECT Gene, Cell, Rsid FROM dice_eqtl") %>% unique
        colnames(eqtl_pairs) = c("Ensgid","Tissue","Rsid")
        eqtl_pairs$eqtl_pair = paste0(eqtl_pairs$Rsid,'-',eqtl_pairs$Ensgid)
        eqtl_snp_gene = eqtl_pairs$eqtl_pair %>% unlist %>% unique
        length(eqtl_snp_gene) %>% print

        if(!is.na(gene_median_tpm)) {
            paste0('* DICE gene mean tpm merged = ') %>% cat
            gct = readRDS(gene_median_tpm)
            nrow(gct) %>% cat
            gct_gather = gct %>% gather("Tissue","TPM",-Ensgid)
        } else gct_gather = NULL
    } else {
        paste0('\n[Error] SQL DB file is mendatory; "gtex_signif_5e-8.db" or "dice_eqtl_5e-08.db".\n') %>% cat
    }

    if(!is.na(gene_median_tpm)) { # debug01 21.01.26
        ' -> transform ' %>% cat
        gct_filt_gene = subset(gct_gather, Ensgid %in% eqtl_genes)$Ensgid %>% unique
        paste0('-> Filt gene = ',length(gct_filt_gene),'\n') %>% cat
    }

    # Filter SNP eQTL genes and their tissues
    paste0('* Filtering eQTLs by SNP: ') %>% cat
    eqtl_pairs_sub = subset(eqtl_pairs,Rsid %in% snps$Rsid)
    eqtl_pairs_uniq = eqtl_pairs_sub$eqtl_pair %>% unique
    paste0('paris = ',length(eqtl_pairs_uniq)) %>% cat
    eqtl_genes_sub = eqtl_pairs_sub$Ensgid %>% unique
    paste0(' (genes = ',length(eqtl_genes_sub)) %>% cat
    eqtl_tissue = eqtl_pairs_sub$Tissue %>% unique %>% sort
    n = length(eqtl_tissue)
    paste0(', tissues = ',n,')') %>% cat
    
    # Preparing gene sets by tissue eQTL genes
    if(statistics=='gsea') {
        ' [' %>% cat; m = n%/%5
        geneset_li = lapply(c(1:n),function(i) {
            if(i%%m==0) '.' %>% cat
            subset(eqtl_pairs_sub,Tissue==eqtl_tissue[i])$Ensgid %>% unique
            })
        names(geneset_li) = eqtl_tissue
        '] ready\n\n' %>% cat
    } else if(statistics=='hyper') {
        geneset_li = eqtl_pairs_uniq
        '; ready\n\n' %>% cat
    }

    dbDisconnect(conn)
    return(list(eqtl_pairs,eqtl_tissue,gct_gather,geneset_li))
}


eqtl_hyper_test = function(
    eqtl_hyper         = NULL,
    f_gwas_snp         = NULL,
    f_eqtl_db          = NULL,
    f_gene_median_tpm  = NULL,
    #perm_n             = 1000,
    out                = 'enrich'
) {
    paste0('\n** Run eqtl_hyper_test function in enrich.r **\n\n') %>% cat
    suppressMessages(library(hypeR))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read files
    data_li = read_files_enrich(f_gwas_snp,f_eqtl_db,f_gene_median_tpm,'hyper')
    eqtl_pairs  = data_li[[1]]
    eqtl_tissue = data_li[[2]]
    gct_gather  = data_li[[3]]
    eqtl_query  = data_li[[4]]

    # Prepare genesets by tissue
    n = length(eqtl_tissue)
    background_num = eqtl_pairs$eqtl_pair %>% unique %>% length
    m = n%/%10; 'Prepare eqtl_pair_set [' %>% cat
    eqtl_pair_set = lapply(c(1:n),function(i) {
        if(i%%m==0) '.' %>% cat
        eqtlpair_tissue_all = subset(eqtl_pairs,Tissue==eqtl_tissue[i])$eqtl_pair %>% unique
        return(eqtlpair_tissue_all)
    })
    names(eqtl_pair_set) = eqtl_tissue
    '] done\n' %>% cat

    # Run hypeR
    paste0('Run hypeR.. ') %>% cat
    hyp_obj = hypeR(
        signature  = eqtl_query,
        genesets   = eqtl_pair_set,
        test       = eqtl_hyper,
        background = background_num
    )
    'done\n' %>% cat

    #hyperRes = data.table::rbindlist(hyperRes_li) %>% as.data.frame
    f_base = basename(f_gwas_snp)
    file_base = tools::file_path_sans_ext(f_base)
    db_name = basename(f_eqtl_db)
    if(db_name=='gtex_signif_5e-8.db') {
        f_name1 = paste0(out,'/gtex-',file_base,'-hyper2.tsv')
    } else if(db_name=='dice_eqtl_5e-08.db') {
        f_name1 = paste0(out,'/dice-',file_base,'-hyper2.tsv')
    }

    #write.table(hyperRes[order(hyperRes$pval), ],f_name1,sep='\t',row.names=F)
    hyp_to_table(hyp_obj, file_path=f_name1)
    paste0('\nWrite file: ',f_name1,'\n') %>% cat
}


gtex_perm_test = function(
    f_gwas_snp        = NULL,
    f_gtex_base       = NULL,
    f_gtex_median_tpm = NULL,
    out               = 'enrich',
    perm_n            = 1000
) {
    paste0('\n** Run gtex_perm_test function in enrich.r **\n\n') %>% cat
    suppressMessages(library(fgsea))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read files
    data_li = read_files_enrich(f_gwas_snp,f_gtex_base,f_gtex_median_tpm,'gsea')
    gtex_pairs  = data_li[[1]]
    eqtl_tissue = data_li[[2]]
    gct_gather  = data_li[[3]]
    geneset_li  = data_li[[4]]

    # Run by tissue
    fgseaRes_li = list()
    n = length(eqtl_tissue)
    for(i in 1:n) {
        # Prepare data
        paste0(i,' ',eqtl_tissue[i],': ') %>% cat
        eqtlgene_tissue_all = subset(gtex_pairs,Tissue==eqtl_tissue[i])$Ensgid %>% unique
        gene_n = length(geneset_li[[i]])
        paste0('genes = ',gene_n,' / ',length(eqtlgene_tissue_all),', ') %>% cat
        
        ranked_ts = subset(gct_gather, Tissue==eqtl_tissue[i])
        ranked_ts_gene = subset(ranked_ts, Ensgid %in% eqtlgene_tissue_all)
        ranked = ranked_ts_gene$TPM
        names(ranked) = ranked_ts_gene$Ensgid
        ranked = ranked[order(ranked_ts_gene$TPM)] # sort by rank
        paste0('exp. = ',length(ranked),'\n') %>% cat

        # Run fgsea (fgsea)
        if(length(ranked)>0) {
            result = fgsea(
                pathways = geneset_li[i],
                stats    = ranked,
                nperm    = perm_n
            )
            fgseaRes_li[[i]] = data.frame(
                result[,1:7],
                queried_gene=gene_n,
                total=length(ranked)
            )
        } else {
            fgseaRes_li[[i]] = NULL
        }
    }
    fgseaRes = data.table::rbindlist(fgseaRes_li) %>% as.data.frame
    # Save result as a file
    f_base = basename(f_gwas_snp)
    file_base = tools::file_path_sans_ext(f_base)
    db_name = basename(f_gtex_base)
    if(db_name=='gtex_signif_5e-8.db') {
        f_name1 = paste0(out,'/gtex-',file_base,'-permn_',perm_n,'.tsv')
        f_name2 = paste0(out,'/gtex-',file_base,'-permn_',perm_n,'.rds')
    } else if(db_name=='dice_eqtl_5e-08.db') {
        f_name1 = paste0(out,'/dice-',file_base,'-permn_',perm_n,'.tsv')
        f_name2 = paste0(out,'/dice-',file_base,'-permn_',perm_n,'.rds')
    }
    
    write.table(fgseaRes[order(fgseaRes$pval), ],f_name1,sep='\t',row.names=F)
    paste0('\nWrite file: ',f_name1,'\n') %>% cat

    fgseaRes_li = list(
        geneset_li = geneset_li,
        ranked     = ranked,
        fgseaRes   = fgseaRes[order(fgseaRes$pval), ]
    )
    saveRDS(fgseaRes_li,f_name2)
    paste0('Write file: ',f_name2,'\n') %>% cat
}

gtex_pm_plot = function(
    f_fgsea_rds = NULL,
    out     = NULL
) {
    paste0('\n** Run gtex_pm_plot function in enrich.r **\n\n') %>% cat
    suppressMessages(library(fgsea))
    suppressMessages(library(ggplot2))
    ifelse(!dir.exists(out), dir.create(out), "")
    fgseaRes_li = readRDS(f_fgsea_rds)
    geneset_li = fgseaRes_li[[1]]
    ranked = fgseaRes_li[[2]]
    fgseaRes = fgseaRes_li[[3]]

    # Prepare tissues for plot
    fgseaRes_sig = fgseaRes[fgseaRes$pval<0.05,]
    dim(fgseaRes_sig) %>% print
    topTissueUp_df = fgseaRes_sig[fgseaRes_sig$ES>0,]
    if(nrow(topTissueUp_df)>0) { 
        topTissueUp = topTissueUp_df[order(topTissueUp_df$pval),]$pathway
    } else topTissueUp = NULL
    topTissueDn_df = fgseaRes_sig[fgseaRes_sig$ES<0,]
    if(nrow(topTissueDn_df)>0) { 
        topTissueDn = topTissueDn_df[order(topTissueDn_df$pval),]$pathway
    } else topTissueDn = NULL
    topTissues  = c(topTissueUp,rev(topTissueDn))

    # Draw plot
    f_base = basename(f_fgsea_rds)
    file_base = tools::file_path_sans_ext(f_base)
    f_name = paste0(out,'/',file_base,'.png')
    png(f_name,width=8,height=12,units='in',res=200)
    plotGseaTable(geneset_li[topTissues],ranked,fgseaRes, gseaParam=0.5)
    dev.off()
    paste0('Draw figure: ',f_name,'\n') %>% cat
}


read_status_file = function(path,db_src) {
    #paste0(db_src,' = ') %>% cat
    if(db_src=='roadmap_bed') {
        status_raw = read.delim(path,header=F,skip=1)
        status_raw = status_raw[,1:4]
        colnames(status_raw) = c('Chr','Start','End','Ann')
    } else if(db_src=='encode_bed') {
        status_raw = read.delim(path,header=F)
        status_raw = status_raw[,1:4]
        colnames(status_raw) = c('Chr','Start','End','Ann')
    } else {
        paste0('\n\n[Error] Unknown option is flagged. Please choose one of [roadmap_bed/encode].\n') %>% cat
        stop()
    }
    return(status_raw)
}


perm_test = function(
    f_gwas_snp = NULL,
    f_status   = NULL,
    db_src     = NULL,
    out        = 'enrich',
    perm_n     = 1000,
    verbose    = NULL
) {
    paste0('\n** Run perm_test function in enrich.r **\n\n') %>% cat
    suppressMessages(library(regioneR))
    ifelse(!dir.exists(out), dir.create(out), "")

    # Read file
    paste0('* Gwas snp = ') %>% cat
    gwas_snp = read.delim(f_gwas_snp, header=F)
    colnames(gwas_snp) = c('Chr','Start','End','Rsid')
    gwas_snp_bed = toGRanges(gwas_snp[,1:4], format="BED")
    dim(gwas_snp) %>% print

    # Run perm_test for each status files
    f_status_paths = list.files(f_status, full.names=T, include.dirs=T)
    n = length(f_status_paths)
    if(n==1) { f_status_paths = f_status; n=1 }
    paste0('* ',n,' files were found from ',f_status,'.\n\n') %>% cat
    cell_types=NULL; zscore_li=list(); pval_li=list(); overlap_li=list()
    source('src/pdtime.r')
    for(i in 1:n) {
        t1=Sys.time()
        file_name = basename(f_status_paths[i])
        if(db_src=='roadmap_bed') {
            cell_type = strsplit(file_name,'_')[[1]][1]
        } else if(db_src=='encode_bed') {
            cell_type = tools::file_path_sans_ext(file_name)
        }
        cell_types = c(cell_types,cell_type)

        ## Read file
        paste0(i,' ',cell_type,': ') %>% cat
        status = read_status_file(f_status_paths[i],db_src)
        #paste0(nrow(status),'; ') %>% cat

        ## Run perm_test
        pt_df = perm_test_calc(gwas_snp_bed,status,db_src,perm_n,verbose)
        zscore_li[[i]]  = data.frame(Status=pt_df$Status, Zscore =pt_df$Zscore)
        pval_li[[i]]    = data.frame(Status=pt_df$Status, Pval   =pt_df$Pval)
        overlap_li[[i]] = data.frame(Status=pt_df$Status, Overlap=pt_df$Overlap)
        paste0(pdtime(t1,2),'\n') %>% cat
    }
    pm_merge = function(x,y) { merge(x=x, y=y, by='Status', all=T) }
    zscore_df  = Reduce(pm_merge, zscore_li)
    pval_df    = Reduce(pm_merge, pval_li)
    overlap_df = Reduce(pm_merge, overlap_li)
    colnames(zscore_df)  = c('Status',cell_types)
    colnames(pval_df)    = c('Status',cell_types)
    colnames(overlap_df) = c('Status',cell_types)
    
    # Write result as files
    gwas_base = basename(f_gwas_snp)
    gwas_f_name = tools::file_path_sans_ext(gwas_base)
    
    f_name1 = paste0(out,'/',db_src,'-',gwas_f_name,'-permn_',perm_n,'-zscore.tsv')
    write.table(zscore_df,f_name1,sep='\t',row.names=F,quote=F)
    paste0('\n* Write file: ',f_name1,'\n') %>% cat
    
    f_name2 = paste0(out,'/',db_src,'-',gwas_f_name,'-permn_',perm_n,'-pval.tsv')
    write.table(pval_df,f_name2,sep='\t',row.names=F,quote=F)
    paste0('* Write file: ',f_name2,'\n') %>% cat

    f_name3 = paste0(out,'/',db_src,'-',gwas_f_name,'-permn_',perm_n,'-overlap.tsv')
    write.table(overlap_df,f_name3,sep='\t',row.names=F,quote=F)
    paste0('* Write file: ',f_name3,'\n') %>% cat
}

perm_test_calc = function(
    gwas_snp_bed = NULL,
    status       = NULL,
    db_src       = NULL,
    perm_n       = 1000,
    verbose      = FALSE
) {
    # Calculate permTest to get z-scores and p-values by chromosome status
    #paste0('permTest - ') %>% cat
    status_ann = status$Ann %>% unique %>% sort
    n = length(status_ann)
    paste0(n, ' [') %>% cat

    # Set background status as universe
    status_all_bed = toGRanges(status[,1:4],format="BED")
    
    check_n = 5; progress = n%/%check_n
    pt_li = lapply(c(1:n),function(i) {
        if(n>=check_n & i%%progress==0) { '.' %>% cat } 
        else if(n<check_n) { '.' %>% cat }
        # Convert data.frame to GRanges format
        status_sub = subset(status,Ann==status_ann[i])
        status_sub_bed = toGRanges(status_sub,format="BED")
        
        pt = permTest(
            A                  = status_sub_bed,
            B                  = gwas_snp_bed,
            ntime              = perm_n,
            randomize.function = resampleRegions,
            universe           = status_all_bed,
            evaluate.function  = numOverlaps,
            force.parallel     = T,
            verbose            = verbose
        )
        pt_df = data.frame(
            Status  = status_ann[i],
            Zscore  = pt$numOverlaps$zscore,
            Pval    = pt$numOverlaps$pval,
            Overlap = pt$numOverlaps$observed
        )
        return(pt_df)
    })
    paste0('] ') %>% cat
    pt_df = data.table::rbindlist(pt_li)
    return(pt_df)
}

## Functions End ##


## Run Function ##
source('src/pdtime.r')
t0 = Sys.time()

if(argv$example) {
    cat(example_msg)
} else if(argv$roadmap_perm) {
    perm_test(
        f_gwas_snp = argv$gwas_snp,
        f_status   = argv$f_roadmap,
        db_src     = argv$db_source,
        perm_n     = argv$perm_n,
        out        = argv$out,
        verbose    = argv$verbose
    )
} else if(argv$heatmap) {
    draw_heatmap(
        f_pmdata = argv$pm_data,
        f_meta   = argv$meta,
        out      = argv$out,
        range    = argv$range,
        annot    = argv$annot,
        fileext  = argv$file_ext
    )
} else if(argv$split_tfbs) {
    split_tfbs(
        f_tfbs = argv$tfbs,
        out    = argv$out
    )
} else if(argv$split_gtex) {
    split_gtex(
        f_gtex = argv$gtex,
        out    = argv$out
    )
} else if(argv$gtex_perm) {
    gtex_perm_test(
        f_gwas_snp  = argv$gwas_snp,
        f_gtex_base = argv$gtex_base,
        f_gtex_median_tpm = argv$gtex_median_tpm,
        out         = argv$out,
        perm_n      = argv$perm_n
    )
} else if(argv$gtex_pm_plot) {
    gtex_pm_plot(
        f_fgsea_rds = argv$fgsea_rds,
        out         = argv$out
    )
} else if(argv$conv_dicevcf2db) {
    pval = argv$pval %>% as.numeric
    conv_dicevcf2db(
        f_dicevcf = argv$dice_vcf,
        out       = argv$out,
        pval      = pval
    )
} else if(eqtl_hyper_n > 0) {
    eqtl_hyper_test(
        eqtl_hyper         = argv$eqtl_hyper,
        f_gwas_snp         = argv$gwas_snp,
        f_eqtl_db          = argv$eqtl_db,
        f_gene_median_tpm  = argv$gene_median_tpm,
        #perm_n             = argv$perm_n,
        out                = argv$out
    )
}

paste0('\n',pdtime(t0,1),'\n') %>% cat