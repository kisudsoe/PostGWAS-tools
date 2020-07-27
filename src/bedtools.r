help_message = '
bedtools, v2020-07-22
This is a function call for generating bedtools command.

Usage:
    Rscript postgwas-exe.r --bedtools bash --base <base file> --out <out folder>


Function:
    bash    Generating bash command scripts to run bedtools

Global arguments:
    --base  <base file path>
            Mendatory. For bash function.
    --out   <out folder path>
            Mendatory. For bash function.
'


## Load global libraries ##
suppressMessages(library(dplyr))


## Functions Start ##
bash_script = function(
    b_path = NULL,
    out    = NULL,
    debug  = FALSE
) {
    # Preparing...
    paste0('\n** Run function: bedtools.r/bash_script... ') %>% cat
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
        'db_gwas/ensembl_gene_hg19.bed',
        'db_gwas/wgEncodeRegTfbsClusteredV3.bed',
        'db_gwas/ucsc_annot.bed'
    )

    # roadmap files
    roadmap_f1 = list.files('db_gwas/roadmap_enh',full.names=F)
    roadmap_f = list.files('db_gwas/roadmap_enh',full.names=T)
    roadmap_f = c(roadmap_f,'db_gwas/roadmap_enh_merge.bed') # Add roadmap merge file

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
    bash = c(bash1,bash2)
    f_name = paste0(out,'/',out,'.sh')
    out_f = file(f_name,"wb") # Set file as Unix type
    write.table(bash,file=out_f,row.names=F,col.names=F,quote=F,sep='')
    paste0('Write bash file: ',f_name,'\n') %>% cat
}


## INIT Function ##
bed_tools = function(
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

    # Run function
    source('src/pdtime.r'); t0=Sys.time()
    if(args$bedtools=='bash') {
        bash_script(b_path,out,debug)
    } else {
        paste0('[Error] There is no such function "',args$bedtools,'" in bedtools.\n') %>% cat
    }
    paste0(pdtime(t0,1),'\n') %>% cat
}