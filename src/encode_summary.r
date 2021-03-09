## Change log ##
# 1/29/2021 Written by Seungsoo Kim, PhD

## Parsing Arguments ##
suppressMessages(library(argparser))
p = arg_parser("Function to generate summary table of encode TFBS")

### Example
p = add_argument(p, '--example',flag=T,
    help="See command examples by functions.")
example_msg = '
Rscript src/encode_summary.r \
	data/encode_dist_raw.tsv
'

### Arguments for this function
p = add_argument(p, '--f_encode_dist',
    help="[Path] TSV file of bedtools closest result of ENCODE data.")
p = add_argument(p, '--out',
    help="[Path] Out folder.")

argv = parse_args(p)


## Load Libraries ##
suppressMessages(library(dplyr))


## Function start ##
source('src/pdtime.r')
t0 = Sys.time() # timer

paste0('\n* Read file: ',argv$f_encode_dist,' = ') %>% cat
encode_dist = read.delim(argv$f_encode_dist,header=F,stringsAsFactors=F)
colnames(encode_dist) = c('chr1','start1','end1','Rsid','chr2','start2','end2','TF','signal','Cells','Dist')
nrow(encode_dist) %>% cat

paste0(' -> Dist=0; ') %>% cat
encode_over = subset(encode_dist,Dist==0)
dim(encode_over) %>% print

paste0('* Summary: ') %>% cat
rsids = encode_over$Rsid %>% unique
n = length(rsids)
paste0(n,' [') %>% cat; m=n%/%10
summ_li = lapply(c(1:n),function(i) {
    if(i%%m==0) '.' %>% cat
    encode_sub  = subset(encode_over,Rsid==rsids[i])
    rsid_coord  = paste0(encode_sub$chr1[1],':',encode_sub$end1[1])
    merge_start = min(encode_sub$start2)
    merge_end   = max(encode_sub$end2)
    merge_coord = paste0(encode_sub$chr2[1],':',merge_start,'-',merge_end)
    tf_cell     = paste0(encode_sub$TF,'(',encode_sub$Cells,')')
    merge_tfs   = paste0(tf_cell,collapse=',')
    data.frame(
        Rsid=encode_sub$Rsid[1],
        Rsid_pos=rsid_coord,
        Encode_tfbs_merged_pos=merge_coord,
        TF_num=nrow(encode_sub),
        TF_cells=merge_tfs)
})
summ_df = data.table::rbindlist(summ_li)
paste0('] = ') %>% cat
dim(summ_df) %>% print

f_name = paste0(argv$out,'/summary_encode.tsv')
write.table(summ_df,f_name,quote=F,row.names=F,sep='\t')
paste0('* Write file: ',f_name,'\n') %>% cat

paste0('\n',pdtime(t0,1),'\n') %>% cat
## Function end ##