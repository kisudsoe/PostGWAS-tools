# setwd('/data')
# ref: https://www.datacamp.com/community/tutorials/sqlite-in-r

library(dplyr)
library(RSQLite)

f_rds = 'db_gwas/gtex_signif_5e-08_Ensgid_dt.rds'

# Read gtex file
gtex = readRDS(f_rds) %>% as.data.frame
dim(gtex) %>% print

# Generate gtex.db file
conn = dbConnect(RSQLite::SQLite(),'db_gwas/db_gwas.db')
dbWriteTable(conn,"gtex",gtex)

# Check db list
dbListTables(conn)
dbGetQuery(conn, "SELECT * FROM gtex LIMIT 10")

# Close the database connection
dbDisconnect(conn)