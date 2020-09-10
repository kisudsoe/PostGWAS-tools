help_message = '
db_filter, v2020-07-23
This is a function call for filtering data.

Usage:
    Rscript postgwas-exe.r --dbfilt ucsc --base <CDS> <Gene> <Promoter> --out <out folder>
    Rscript postgwas-exe.r --dbfilt roadmap --base <base file(s)> --out <out folder> --enh TRUE
    Rscript postgwas-exe.r --dbfilt roadmap --base <base file(s)> --out <out folder> --enh FALSE
    Rscript postgwas-exe.r --dbfilt roadmap --base <base file(s)> --out <out folder> --enh TRUE --sep TRUE
    ...
    Rscript postgwas-exe.r --dbfilt gtex_ovl --base <base file> --gtex <GTEx RDS file> --out <out folder>
    Rscript postgwas-exe.r --dbfilt gtex_ovl --base <base file> --gtex <GTEx RDS file> --out <out folder> --tissue <optional:GTEx tissue name>


Functions:
    ucsc        Compile the ucsc downloaded promoter/gene/cds region annotations.
    roadmap     Filtering Roadmap data by enhancer tags.
    gtex        Filtering GTEx data by eQTL p-value.
    gtex_ovl    Overlapping the GTEx data with the input GWAS SNPs.
    hic_bed     --base <HiCCUPS files> --out <out folder>
                Converting the HiCCUPS data to the BED format.
    dist        Filtering distance data from Bedtools closest function.
    regulome    Filtering and overlapping by Regulome score ≥2b.
    lnc_ovl     Overlapping the lncRNASNP2 data with the input GWAS SNPs.

Global arguments:
    --base      <base file/folder>
                Base file/folder path is mendatory.
                For ucsc function, you have to input three UCSC downloaded BED files by this order:
                  [1] cds region file, [2] whole gene region file, [3] proximal promoter region file
    --out       <out folder>
                Out folder path is mendatory. Default is "db" folder.

Required arguments:
    --ctype     <cell type id>
                An optional argument for the "roadmap" function.
                See cell-type number information at
                https://github.com/mdozmorov/genomerunner_web/wiki/Roadmap-cell-types.
    --enh       <default: TRUE>
                An optional argument for the "roadmap" function to filter enhancer regions.
    --sep       <default: FALSE>
                An optional argument for the "roadmap" function to generate cell-type seperated results.
    --meta      <roadmap meta file path>
                An optional argument for the "roadmap", "dist" function.
                For "roadmap" function, this argument needs "--sep TRUE" argument.
                Output file will be organized by the cell types.
    --infotype  <default: FALSE>
                An optional argument for the "dist" function.
                If input as "ucsc", three output files will be generated.
                  [1] cds region, [2] whole gene region, [3] proximal promoter region
                If input as "tags", This option allows to generate seperated output files by the tags.
    --pval      <p-value threshold>
                A required argument for the "gtex" function to filter significant eQTLs.
    --gtex      <Filtered GTEx RDS or SQL DB file path>
                A required argument for the "gtex_ovl" function to overlap GWAS SNPs with GTEx data.
    --tissue    <GTEx tissue name>
                An optional argument for the "gtex_ovl" function to filter a specific tissue.
    --regulm    <Regulome data folder>
                A required argument for the "regulome" function to load the data.
    --lncrna    <lncRNASNP2 data folder>
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
hic_bed = function(
    f_hics = NULL,  # Input Hi-C file
    out    = 'data' # Out folder path
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/hic_bed... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    'Ready\n' %>% cat

    # If the base path is folder, get the file list
	paths1 = list.files(f_hics,full.name=T)
	n = length(paths1)
    if(n>0) {
        paste0(n,' Files/folders input.\n') %>% cat
        paths = NULL
        for(i in 1:n) {
            paths2  = list.files(paths1[i],full.name=T)
            if(length(paths2)==0) {
                paths = c(paths,paths1[i])
                paste0('  ',i,' ',paths1[i],'\n') %>% cat
            } else {
                paths = c(paths,paths2)
                paste0('  ',i,' ',length(paths2),' files in the ',
                    paths1[i] %>% basename,'\n') %>% cat
            }
        }
    } else paths = f_hics
    fl_name = basename(paths)
    file_nm = tools::file_path_sans_ext(fl_name)
	paste0('Total ',length(paths),' file(s) is/are input.\n') %>% cat

    # Multiple process
    n = length(paths)
    o = lapply(c(1:n),function(i) {
        paste0('\n  ',i,'\t',fl_name[i],'\t') %>% cat
        # Read Hi-C data file
        hic_full = data.table::fread(paths[i])
        dim(hic_full) %>% print
        hic = hic_full[,1:6]
        colnames(hic) = c('chr1','x1','x2','chr2','y1','y2')

        # Convert to BED format
        m = nrow(hic)
        loops_li = lapply(c(1:m),function(j) {
            hic_row = hic[j,] %>% unlist
            chrs = paste0('chr',c(hic_row[1],hic_row[4]))
            loop = paste0('loop',j,'.',c('x','y'))
            out  = data.frame(
                chr   = chrs,
                start = c(hic_row[2],hic_row[5]),
                end   = c(hic_row[3],hic_row[6]),
                name  = loop
            )
            return(out)
        })
        loops_df = data.table::rbindlist(loops_li)
        #loops_df = loops_df[order(loops_df$chr,loops_df$start),] #<- Doesn't work..
        
        # Save as BED file
        f_name = paste0(out,'/',file_nm[i],'.bed')
        write.table(loops_df,f_name,sep='\t',row.names=F,col.names=F,quote=F)
        paste0('  Write file: ',f_name,'\n') %>% cat
    })
}

ucsc_compile = function(
    f_ucsc = NULL,  # Input UCSC downloaded BED file paths
    out    = 'data' # Out folder path
) {
    # Load function-specific library

    # Preparing...
    paste0('\n** Run function: db_filter.r/ucsc_compile... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    'Ready\n' %>% cat

    # Read UCSC BED files
    paste0('  Read CDS file, dim\t\t= ') %>% cat
    cds  = read.delim(f_ucsc[1],head=F,stringsAsFactors=F); dim(cds) %>% print
    paste0('  Read Gene file, dim\t\t= ') %>% cat
    gene = read.delim(f_ucsc[2],head=F,stringsAsFactors=F); dim(gene) %>% print
    paste0('  Read Promoter file, dim\t= ') %>% cat
    prom = read.delim(f_ucsc[3],head=F,stringsAsFactors=F); dim(prom) %>% print

    # Generate tags
    paste0('\n  Generate tags... ') %>% cat
    cds_tag = apply(cds,1,function(row) {
        row_split = strsplit(row[4],"\\_")[[1]]
        row_tag   = paste0(row_split[1:4],collapse='_')
        return(row_tag)
    }) %>% unlist
    '1.. ' %>% cat

    gene_tag = paste0(gene[,4],'_wholeGene')
    '2.. ' %>% cat

    prom_tag = apply(prom,1,function(row) {
        row_split = strsplit(row[4],"\\_")[[1]]
        row_tag   = paste0(row_split[1],'_proximalPromoter')
        return(row_tag)
    }) %>% unlist
    '3.. done\n' %>% cat

    # Compile the three UCSC data
    paste0('  Compile the three UCSC data... ') %>% cat
    ucsc_bed_li = list(
        cds_bed  = data.frame(cds[,1:3], name=cds_tag),
        gene_bed = data.frame(gene[,1:3],name=gene_tag),
        prom_bed = data.frame(prom[,1:3],name=prom_tag)
    )
    ucsc_bed = data.table::rbindlist(ucsc_bed_li)#,fill=TRUE
    colnames(ucsc_bed) = c('chr','start','end','name')
    dim(ucsc_bed) %>% print

    #na_count = apply(ucsc_bed,1,function(row) {
    #    return((is.na(row)))
    #}) %>% unlist
    #which_na = which(na_count)
    #ucsc_na  = ucsc_bed[which_na,]
    #print(ucsc_na)
    #stop()

    # Write BED file
    f_name = paste0(out,'/ucsc_annot.bed')
    write.table(ucsc_bed,f_name,sep='\t',row.names=F,col.names=F,quote=F)
    paste0('Write a BED file: ',f_name,'\n') %>% cat
}

lncrna_overlap = function(
    snp_path = NULL,  # Input GWAS SNP file path
    lnc_path = NULL,  # lncRNASNP2 data downloaded folder path
    out      = 'data' # Out folder path. Will generate 'summary' subfoleder.
) {
    suppressMessages(library(RSQLite))

    # Preparing...
    paste0('\n** Run function: db_filter.r/lncrna_overlap... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    out_summary = paste0(out,'/summary')
    ifelse(!dir.exists(out_summary), dir.create(out_summary),'')
    'Ready\n' %>% cat

    # Read SNP file
    snp = read.delim(snp_path,header=F)
    colnames(snp) = c('chr','start','end','rsid')
    dbsnp  = gsub('(.*)_.*','\\1',snp[,4])
    snp_df = cbind(snp,dbsnp)
    paste0('Input GWAS SNPs N = ',length(dbsnp),'\n') %>% cat

    # Query SQL db
    f_ext = tools::file_ext(lnc_path)
    if(f_ext=='db') {
        paste0('  Query ',basename(lnc_path),' = ') %>% cat
        conn = dbConnect(RSQLite::SQLite(),lnc_path)
        str_snp = paste0("'",dbsnp,"'", collapse=',')
        rsids_query = sprintf('SELECT * FROM lncrna WHERE dbsnp IN (%s)',str_snp)
        lncsnp = dbGetQuery(conn,rsids_query)
        dim(lncsnp) %>% print
    } else paste0('\n[Error] Input lnc_path seems not SQL db file path. STOP.\n') %>% cat
    
    # Overlapping lncRNA data with input GWAS SNPs
    snp_lnc = merge(snp_df,lncsnp,by='dbsnp')
    paste0('  Summ =\n') %>% cat
    data.frame(
        lncRNA = unique(snp_lnc$lncRNA) %>% length,
        SNPs   = unique(snp_lnc$dbsnp) %>% length
    ) %>% print

    # Save as a BED file
    f_name1 = paste0(out_summary,'/snp_lncrnasnp_',unique(snp_lnc$dbsnp)%>%length,'.bed')
    write.table(snp_lnc[,2:5],f_name1,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('\n  Write file: ',f_name1,'\n') %>% cat

    # Save as a TSV file
    f_name2 = paste0(out,'/lncrnasnp_',unique(snp_lnc$dbsnp)%>%length,'.tsv')
    write.table(snp_lnc,f_name2,row.names=F,quote=F,sep='\t')
    paste0('  Write file: ',f_name2,'\n') %>% cat
    dbDisconnect(conn)
}


regulome_filt = function(
    snp_path = NULL,  # Input GWAS SNP file path
    reg_path = NULL,  # Regulome data downloaded folder path
    out      = 'data' # Out folder path. Will generate 'summary' subfolder.
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/regulome_filt... ') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    out_summary = paste0(out,'/summary')
    ifelse(!dir.exists(out_summary), dir.create(out_summary),'')
    'Ready\n' %>% cat

    # Read a GWAS SNP BED file input
    snp   = read.delim(snp_path,header=F)
    rsids = gsub('(.*)_.*','\\1',snp[,4])
    paste0('Input GWAS SNPs N\t= ') %>% cat; length(rsids) %>% print

    # Read the Regulome data files
    f_reg = list.files(reg_path)
    paths = paste0(reg_path,'/',f_reg)
    n = length(paths)
    paste0(n,' Regulome data load...\n') %>% cat
    reg_li = lapply(paths,function(path) {
        paste0('  Read: ',path,'; ') %>% cat
        f_ext = tools::file_ext(path)
        if(f_ext=='gz') out =  try(read.delim(gzfile(path),header=F))
        else if(f_ext=='rds') out = try(readRDS(path))
        else out = try(read.delim(path,header=F))
        colnames(out) = c('chr','id','rsid','description','level')
        paste0('dim = ') %>% cat; dim(out) %>% print
        return(out)
    })
    reg_df = data.table::rbindlist(reg_li)
    
    # Filter by score ≥2b
    reg_1f_only = subset(reg_df,level == '1f') %>% unique
    reg_2b = subset(reg_df,level %in% c("1a","1b","1c","1d","1e","1f","2a","2b")) %>% unique
    paste0('\n  Regulome score >=2b, SNPs\t\t= ') %>% cat
    nrow(reg_2b) %>% print
    paste0('  Functional motifs (1a~2b - 1f only)\t= ') %>% cat
    nrow(reg_2b)-nrow(reg_1f_only) %>% print

    # GWAS SNP counts
    snp_1f_only = subset(reg_1f_only,rsid %in% rsids) %>% unique
    snp_2b      = subset(reg_2b,rsid %in% rsids) %>% unique
    paste0('\n  Regulome >=2b, GWAS SNPs\t\t= ') %>% cat
    nrow(snp_2b) %>% print
    paste0('  GWAS SNPs occupied in
    functional motifs (1a~2b - 1f only)\t= ') %>% cat
    nrow(snp_2b) - nrow(snp_1f_only) %>% print

    # Save as TSV file
    f_name1 = paste0(out,'/regulome_',nrow(snp_2b),'.tsv')
    write.table(snp_2b,f_name1,row.names=F,quote=F,sep='\t')
    paste0('\nWrite file: ',f_name1,'\n') %>% cat

    # Save as BED file
    snp_rsid = data.frame(snp,rsid=rsids)
    snp_bed  = subset(snp_rsid,rsid %in% snp_2b$rsid) %>% unique
    f_name2  = paste0(out_summary,'/snp_regulome2b_',nrow(snp_bed),'.bed')
    write.table(snp_bed[,1:4],f_name2,row.names=F,col.names=F,quote=F,sep='\t')
    paste0('Write file: ',f_name2,'\n\n') %>% cat
}


distance_filt = function(
    f_path   = NULL,   # Bedtools closest result paths
    out      = 'data', # Out folder path
    infotype = NULL,   # (1) ucsc: Input BED file is UCSC annotation file.
                       # (2) tags: Input BED file has unique tags as names.
                       # (3) NULL: Input BED file has general names. 
    debug
) {
    # Preparing..
    paste0('\n** Run function: db_filter.r/distance_filt... ') %>% cat
    f_name = tools::file_path_sans_ext(basename(f_path))
    'Ready\n' %>% cat

    paste0('  File ',f_name,'... ') %>% cat
    rd = data.table::fread(f_path) %>% as.data.frame
    paste0('nrow= ',nrow(rd),'.. ') %>% cat

    # Extract data
    pos = apply(rd,1,function(row) paste0(row[5],':',row[6],'-',row[7])) %>% unlist
    rd_ = cbind(rd[,1:4],pos,rd[,c(8,ncol(rd))])
    colnames(rd_) = c('chr','start','end','rsid','ann_pos','tag','dist')
    dis = rd_$dist %>% as.character %>% as.numeric
    rd_df = cbind(rd_[,1:6],dis) %>% unique
    paste0('done\n') %>% cat

    # Filtering by distance 0 from the annotations
    rd_df2 = subset(rd_df,dis==0) %>% unique
    paste0('  Annotations occupied by SNPs\t= ') %>% cat; unique(rd_df2$ann_pos) %>% length %>% print
    paste0('  SNPs in annotations\t\t= ') %>% cat; unique(rd_df2$rsid) %>% length %>% print
    
    # Save as a TSV file
    #f_name1 = paste0(out,'/',f_name,'_filt.tsv')
    #write.table(rd_df2,f_name1,row.names=F,quote=F,sep='\t')
    #paste0('  Write file: ',f_name1,'\n') %>% cat

    # Whether dealing with UCSC annotations:
    if(!is.null(infotype)) {
        # List tags
        #x <- "a1~!@#$%^&*(){}_+:\"<>?,./;'[]-=" #or whatever
        tag_alter = stringr::str_replace_all(rd_df2$tag, "[[:punct:]]", ".") # Remove special characters
        if(infotype=='ucsc') {
            paste0('\n  UCSC annotations: ') %>% cat
            tags = lapply(tag_alter,function(tag) strsplit(tag,'\\.')[[1]][2]) %>%
                unlist %>% as.factor
        } else if(infotype=='tags') {
            paste0('\n Annotations: ') %>% cat
            tags = tag_alter %>% as.factor
        } else {
            paste0('\n[ERROR] You chose wrong infotype option: "',infotype,'". 
  You should choose one of these: "ucsc", "tags" or NULL.\n') %>% cat
            quit()
        }

        tag_level = levels(tags)
        n = length(tag_level)
        paste0(n,' tags\n') %>% cat
        o = lapply(c(1:n),function(i) {
            # Extract by each tag
            paste0('    ',i,' ',tag_level[i],':\t') %>% cat
            which_row  = which(tags==tag_level[i])
            rd_df2_tag = rd_df2[which_row,]
            length(which_row) %>% cat

            # Save as a BED file
            snp_bed = rd_df2_tag[,1:4] %>% unique
            snp_n   = unique(snp_bed$rsid) %>% length
            f_name2 = paste0(out,'/snp_',infotype,'_',tag_level[i],'_',snp_n,'.bed')
            write.table(snp_bed,f_name2,row.names=F,col.names=F,quote=F,sep='\t')
            paste0('.. Save at: ',f_name2,'\n') %>% cat
        })
    } else {
        # Save as a BED file
        snp_bed = rd_df2[,1:4] %>% unique
        snp_n   = unique(snp_bed$rsid) %>% length
        if(snp_n>0) {
            ifelse(!dir.exists(out), dir.create(out),'')
            f_name2 = paste0(out,'/snp_',f_name,'_',snp_n,'.bed')
            write.table(snp_bed,f_name2,row.names=F,col.names=F,quote=F,sep='\t')
            paste0('  Write file: ',f_name2,'\n') %>% cat
        } else if(snp_n==0) {
            paste0('  [SKIP] SNP N = 0\n') %>% cat
        }
    }
}


distance_filt_multi = function(
    f_paths  = NULL,   # Input file/folder paths
    out      = 'data', # Out folder path
    meta     = NULL,   # Organizing output files by the designated groups
    infotype = NULL,   # Option for "distance_filt" function. (1) ucsc, (2) tags
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/distance_filt_multi...\n') %>% cat
    paste0('Input file/folder N\t= ') %>% cat; length(f_paths) %>% print
    # If the base path is folder, get the file list
    if(length(f_paths)==1) {
        paths = list.files(f_paths,full.name=T)
        if(length(paths)==0) paths = f_paths
    } else paths = f_paths
    paste0('Input file N\t= ') %>% cat; length(paths) %>% print

    # For metadata,
    if(!is.null(meta)) {
        base_f = basename(paths)
        meta_dat = read.delim(meta,stringsAsFactors=F)
        sub_dir = paste0(out,'/',meta_dat$groups)
        paste0('  Read metadata file dim\t= ') %>% cat; dim(meta_dat) %>% print
    }

    # Run function by each file path
    ifelse(!dir.exists(out), dir.create(out),'')
    n = length(paths)
    o = lapply(c(1:n),function(i) {
        # find metadata by EID
        if(!is.null(meta)) {
            j = which(meta_dat$f_names==base_f[i])
            if(length(j)==0) {
                paste0('[BREAK] ',base_f[i],' is not existing in meta file.\n') %>% cat
                return(NULL)
            }
            # mkdir by groups column in meta_dat file
            #ifelse(!dir.exists(sub_dir[j]), dir.create(sub_dir[j]),'')
            out = sub_dir[j]
        }

        # Run function
        source('src/pdtime.r'); t0=Sys.time()
        f_path = paths[i]
        distance_filt(f_path,out,infotype,debug)

        # Print process
        if(i%%10==0) paste0('  ',i,'/',n,' being processed.\n') %>% cat
        if(n>1&n<10) paste0('  ',i,'/',n,' done: ',f_path,'\n') %>% cat
        if(n>1)      paste0(pdtime(t0,2),'\n\n') %>% cat
    })
}


readgtex = function(
    rsids  = NULL,
    f_gtex = NULL
) {
    f_ext = tools::file_ext(f_gtex)
    if(f_ext %in% c('gz','rds')) {
        paste0('Read, ',basename(f_gtex),' = ') %>% cat
        if(f_ext=='gz') {
            gtex = fread(gzfile(f_gtex),header=T)
        } else if(f_ext=='rds') {
            gtex = readRDS(f_gtex)
        } else {
            paste0('[ERROR] Unknown GTEx file format. gz/rds file format is required.')
            quit
        }
        dim(gtex) %>% print

        # Overlapping eQTLs
        paste0('  Overlapped eQTL-gene pairs = ') %>% cat
        eqtls = subset(gtex, Rsid %in% rsids)
        nrow(eqtls) %>% print
    } else if(f_ext=='db') {
        # Query SNPs from GTEx SQL DB file
        paste0('SNP query to ',basename(f_gtex),' = ') %>% cat
        conn = dbConnect(RSQLite::SQLite(),f_gtex)
        str_snp = paste0("'",rsids,"'",collapse=",")
        rsids_query = sprintf("SELECT * FROM gtex WHERE Rsid IN (%s)",str_snp)
        eqtls = dbGetQuery(conn,rsids_query)
        dim(eqtls) %>% print 
        dbDisconnect(conn)
    }
    return(eqtls)
}


gtex_overlap = function(
    f_path    = NULL,   # Input SNPs BED file path
    f_gtex    = NULL,   # Filtered GTEx RDS file path
    out       = 'data', # Out folder path. Will generate gtex_eqlt subfoler.
    tissue_nm = NULL,   # Optional tissue name
    debug
) {
    suppressMessages(library(data.table))
    suppressMessages(library(RSQLite))

    # Preparing...
    paste0('\n** Run function: db_filter.r/gtex_overlap...\n') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    out_gtex_eqtl = paste0(out,'/gtex_eqtl')
    'ready\n' %>% cat

    # Read SNP file
    snp = read.delim(f_path,header=F)
    colnames(snp) = c('chr','start','end','ann')
    rsids = gsub('(.*)_.*','\\1',snp[,4])
    paste0('Input GWAS SNPs N = ',length(rsids),'\n') %>% cat

    # Load filtered GTEx RDS file
    eqtls = readgtex(rsids,f_gtex)


    # Filter by tissue (optional)
    if(!is.null(tissue_nm)) {
        paste0('\n[Option] ',tissue_nm) %>% cat
        eqtls = subset(eqtls,tissue == tissue_nm)
        f_name1 = paste0(out,'/gtex_signif_',tissue_nm,'_',unique(eqtls$Rsid)%>%length,'.tsv')
        paste0(', dim = ') %>% cat; dim(eqtls) %>% print
    } else {
        f_name1 = paste0(out,'/gtex_signif_',unique(eqtls$Rsid)%>%length,'.tsv')
    }
    paste0('  eQTLs N = ') %>% cat; unique(eqtls$Rsid) %>% length %>% print
    paste0('  Associated eGenes = ') %>% cat; unique(eqtls$Ensgid) %>% length %>% print

    # Save as a TSV file
    if(nrow(eqtls)>0) {
        write.table(eqtls,f_name1,row.names=F,quote=F,sep='\t')
        paste0('\nWrite file: ',f_name1,'\n') %>% cat
    } else {
        paste0('\n[Warning] No overlapped eQTL was found. Please check query SNPs in GTEx homepage.\n\n') %>% cat
        return(NULL)
    }
    

    # Generate as BED file
    snp_rsid = data.frame(snp,rsid=rsids)
    tissue_names = unique(eqtls$tissue)
    tissue_n = length(tissue_names)
    paste0('Generating BED files at "',out_gtex_eqtl,
        '" for ',tissue_n,' tissues.. ') %>% cat
    o = lapply(c(1:tissue_n),function(i) {
        # Filter by tissue
        tissue_name  = tissue_names[i]
        eqtls_tissue = subset(eqtls,tissue==tissue_name)
        snp_bed = subset(snp_rsid,rsid %in% eqtls_tissue$Rsid)[,1:4] %>% unique
        if(debug) {
            paste0('\n',tissue_name,'\n  eQTL BED, dim = ') %>% cat
            dim(snp_bed) %>% print
            paste0('  eQTL SNP N\t\t= ') %>% cat
            unique(snp_bed$ann) %>% length %>% print
        }

        # Save as a BED file
        f_name2 = paste0(out_gtex_eqtl,'/snp_gtex_',tissue_name,
            '_',unique(snp_bed$ann)%>%length,'.bed')
        write.table(snp_bed,f_name2,row.names=F,col.names=F,quote=F,sep='\t')
        if(debug) paste0('\nWrite file: ',f_name2,'\n') %>% cat
    })
    'done\n\n' %>% cat
}

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
    f_path = NULL,  # Download folder path
    out    = 'db',  # Out folder path
    ctype  = NULL,  # Optional: Cell type ID
    enh    = TRUE,  # Optional: Filtering enhancer (Default: TRUE)
    sep    = FALSE, # Optional: Cell type separated results (Defualt: FALSE)
    meta   = NULL,  # Optional: Meta data file
    debug
) {
    # Preparing...
    paste0('\n** Run function: db_filter.r/roadmap_filt...\n') %>% cat
    ifelse(!dir.exists(out), dir.create(out),'')
    
    if(!is.null(ctype)) {
        cid    = formatC(ctype,width=3,flag='0') %>% as.character # Convert number to '###' format
        f_name = paste0(out,'/roadmap_',ctype,'_enh.bed')
    } else {
        cid    = formatC(c(1:129),width=3,flag='0') %>% as.character # 001~129
        if(enh) { f_name = paste0(out,'/roadmap_enh.bed')
        } else    f_name = paste0(out,'/roadmap_total.bed')
    }
    if(!is.null(meta)) {
        meta_dat = read.delim(meta)
        sub_dir = paste0(out,'/',meta_dat$ANATOMY)
        paste0('  Read metadata file dim\t= ') %>% cat; dim(meta_dat) %>% print
    }

    # Read Roadmap files
    paste0('  Reading files..\n') %>% cat
    n = length(cid)
    road_li = lapply(c(1:n), function(i) {
        #path = paste0(f_path,'E',cid[i],'_25_imputed12marks_hg38lift_dense.bed.rds') # hg38 by liftover
        path = paste0(f_path,'/E',cid[i],'_25_imputed12marks_dense.bed.rds') # hg19 original
        road = try(readRDS(path))
        if('try-error' %in% class(road)) { # if file is not found.
            paste0('  ',path,' - file not found.\n') %>% cat
            road_enh = NULL
        } else if(enh) { # If file exist.
            road_enh = subset(road,
                name %in% c("13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc"))
        } else road_enh = road

        # Print process
        if(n<10) {
            paste0('  ',i,' ',path,' ') %>% cat
            dim(road_enh) %>% print
        } else if(i%%10==0) paste0('    ',i,'/',n,' being processed.\n') %>% cat

        # For --sep and --meta arguments
        if(sep) {
            road_enh = road_enh[,c(1:3,6)]
            if(!is.null(meta)) {
                # find metadata by EID
                E_cid = paste0('E',cid[i])
                j = which(meta_dat$EID==E_cid)
                if(length(j)==0) {
                    paste0('[BREAK] ',E_cid,' is not existing in meta file.')
                    return(NULL)
                }
                # mkdir by meta_dat$ANATOMY
                ifelse(!dir.exists(sub_dir[j]), dir.create(sub_dir[j]),'')
                f_name = paste0(sub_dir[j],'/roadmap_',cid[i],'_enh.bed')
            } else {
                f_name = paste0(out,'/roadmap_',cid[i],'_enh.bed')
            }
            write.table(road_enh,f_name,row.names=F,col.names=F,quote=F,sep='\t')
        } else {
            return(road_enh)
        }
    })

    if(sep) {
        paste0('  Finished processing ',n,' files.\n\n') %>% cat
    } else {
        paste0('  Finished reading and filtering ',n,' files.\n\n') %>% cat
        road_enh = data.table::rbindlist(road_li)
        rm(road_li)

        # Save file
        road_enh = road_enh[,c(1:3,6)]
        write.table(road_enh,f_name,row.names=F,col.names=F,quote=F,sep='\t')
        paste0('Write file: ',f_name,'\n') %>% cat
    }
}


db_filter = function(
    args = NULL
) {
    if(length(args$help)>0) {     help     = args$help
    } else                        help     = FALSE
    if(help) { cat(help_message); quit() }

    # Global arguments
    if(length(args$base)>0)       b_path   = args$base
    if(length(args$out)>0)        out      = args$out
    if(length(args$debug)>0) {    debug    = args$debug
    } else                        debug    = FALSE

    # Required arguments
    if(length(args$meta)>0) {     meta     = args$meta
    } else                        meta     = NULL
    if(length(args$ctype)>0) {    ctype    = args$ctype
    } else                        ctype    = NULL
    if(length(args$enh)>0) {      enh      = args$enh
    } else                        enh      = TRUE
    if(length(args$sep)>0) {      sep      = args$sep
    } else                        sep      = FALSE
    if(length(args$pval)>0)       pval     = args$pval
    if(length(args$gtex)>0)       gtex     = args$gtex
    if(length(args$tissue)>0) {   tissue   = args$tissue
    } else                        tissue   = NULL
    if(length(args$regulm)>0)     reg_path = args$regulm
    if(length(args$lncrna)>0)     lnc_path = args$lncrna
    if(length(args$infotype)>0) { infotype = args$infotype
    } else                        infotype = NULL

    # Run function
    source('./src/pdtime.r'); t0=Sys.time()
    if(args$dbfilt == 'ucsc') {
        ucsc_compile(b_path,out)
    } else if(args$dbfilt == 'roadmap') {
        roadmap_filt(b_path,out,ctype,enh,sep,meta,debug)
    } else if(args$dbfilt == 'gtex') {
        gtex_filt(b_path,out,pval,debug)
    } else if(args$dbfilt == 'gtex_ovl') {
        gtex_overlap(b_path,gtex,out,tissue,debug)
    } else if(args$dbfilt == 'dist') {
        distance_filt_multi(b_path,out,meta,infotype,debug)
    } else if(args$dbfilt == 'regulome') {
        regulome_filt(b_path,reg_path,out)
    } else if(args$dbfilt == 'lnc_ovl') {
        lncrna_overlap(b_path,lnc_path,out)
    } else if(args$dbfilt == 'hic_bed') {
        hic_bed(b_path,out)
    } else {
        paste0('[Error] There is no such function "',args$dbfilt,'" in db_filte.\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n\n') %>% cat
}