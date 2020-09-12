#!/usr/bin/env Rscript

# 1. Run docker image
## (CMD)> docker run -it -v "C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3:/data" -v "C:\Github\PostGWAS-tools:/git" kisudsoe/postgwas-env
# 2. Run this script
## (bash)# Rscript /git/src/snippet-sqlite.r

# RSQLite tutorial, ref: https://www.datacamp.com/community/tutorials/sqlite-in-r

library(dplyr)
library(RSQLite)

#f_rds   = 'db_gwas/gtex_signif_5e-08_Ensgid_dt.rds'

slnc_path = '/data/db_gwas/lncrna/lncRNASNP2_snplist.txt.rds'
ann_path  = '/data/db_gwas/lncrna/lncrnas.txt.rds'
dis_path  = '/data/db_gwas/lncrna/lncrna-diseases_experiment.txt.rds'
f_db    = '/data/db_gwas/db_gwas.db'
db_name = 'lncrna'

# Read rds files
#gtex = readRDS(f_rds) %>% as.data.frame
#dim(gtex) %>% print
paste0('  Read: ',slnc_path,' = ') %>% cat
snplnc = readRDS(slnc_path)
dim(snplnc) %>% print

paste0('  Read: ',ann_path,' = ') %>% cat
ann = readRDS(ann_path)
colnames(ann)[1] = 'lncRNA'
dim(ann) %>% print

paste0('  Read: ',dis_path,' = ') %>% cat
dis = readRDS(dis_path)
dim(dis) %>% print

snp_lnc_ann = merge(snplnc,ann,by='lncRNA') %>% unique
data_in = merge(snp_lnc_ann,dis,by='lncRNA',all.x=T) %>% unique


# Generate gtex.db file
conn = dbConnect(RSQLite::SQLite(),f_db)
dbWriteTable(conn, db_name, data_in)

# Check db list and data
dbListTables(conn) %>% print
my_query = paste0('SELECT * FROM ',db_name,' LIMIT 10')
dbGetQuery(conn, my_query) %>% print

# Close the database connection
dbDisconnect(conn)