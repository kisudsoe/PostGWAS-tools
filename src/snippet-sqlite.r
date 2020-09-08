# setwd('/data')
# ref: https://www.datacamp.com/community/tutorials/sqlite-in-r

library(dplyr)
library(RSQLite)

f_rds = 'db_gwas/gtex_signif_5e-08_Ensgid_dt.rds'
<<<<<<< HEAD
f_db = 'db_gwas/db_gwas.db'
=======
>>>>>>> 6cc62cfa7e5126607873a865a511af29b258b319

# Read gtex file
gtex = readRDS(f_rds) %>% as.data.frame
dim(gtex) %>% print

# Generate gtex.db file
conn = dbConnect(SQLite(),f_db)
dbWriteTable(conn,"gtex",gtex)

# Check db list
dbListTables(conn)
<<<<<<< HEAD

# Disconnect DB
=======
dbGetQuery(conn, "SELECT * FROM gtex LIMIT 10")

# Close the database connection
>>>>>>> 6cc62cfa7e5126607873a865a511af29b258b319
dbDisconnect(conn)