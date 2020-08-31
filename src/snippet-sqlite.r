# setwd('/data')

library(dplyr)
library(RSQLite)

# Read gtex file
gtex = readRDS('db_gwas/gtex_signif_5e-08_Ensgid_dt.rds')
dim(gtex) %>% print

# Convert data.table to df
gtex  = gtex %>% as.data.frame

# Generate gtex.db file
conn = dbConnect(RSQLite::SQLite(),'db_gwas/db_gwas.db')
dbWriteTable(conn,"gtex",gtex)

# Check db list
dbListTables(conn)