library(tidyverse)
gse_id <- "GSE102349"

## read in gene expression data, there shoube be 24530 genes and 113 samples
load_data <- function(){
  dat <-
    as.data.frame(
      read.table(
        "analysis/data/raw_data/GSE102349_NPC_mRNA_processed.txt",
        header = TRUE,
        sep = "\t",
        dec = ".",
        na.strings = "NA",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  
  row.names(dat) <- dat[,1]
  dat<-dat[,-1]
  head(dat)
  dim(dat)
  return(dat)
}

dat <- load_data()


load_clincal <- function(){
  #check missing value
  f<-function(x) sum(is.na(x))
  #1 = summarize by row; 2= summarize by column
  d<-as.data.frame(apply(dat,2,f))
  d$id=row.names(d)
  head(d[order(d[,1],decreasing=TRUE),])
  
  #convert missing value to 0
  dat[is.na(dat)] <- 0
  
  # read in sample information 
  # Sample ID, "event", "time to event" and "clinical stage" gain
  # see file GSE102349_series_matrix_survival.csv, clinical file should involved Sample ID, "event", "time to event" and "clinical_stage" extracted from file GSE102349_series_matrix.txt.gz. Here we also group samples based on clinical stage.
  
  # GSE102349_series_matrix_survival.csv
  
  
  info <-
    as.data.frame(
      read.table(
        "analysis/data/raw_data/GSE102349_series_matrix_survival.csv",
        header = TRUE,
        sep = ",",
        dec = ".",
        na.strings = "NA",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  row.names(info) <- info[, 1]
  info <- info[, -1]
  dim(info)
  
  #remove missing data and filter out clinical stage I and II samples
  info=info[!info[, "time"] == "N/A",]
  # info=info[info[, "clinical_stage"] == "III" | info[, "clinical_stage"] == "IV",]
  
  #set time to numeric datatype
  info$time=as.numeric(as.character(info$time))
  
  dim(info)
  head(info)
  table(info$`clinical_stage`)
  return(info)
}

clin <- load_clincal()

# clean clin -------------------------------------------------------------
clin %>% glimpse()
clin %>%   count(event)

clean_clin <-
  clin %>%
  rownames_to_column("sample") %>%
  # mutate(time = time * 30) %>% 
  mutate(
    days_to_death = if_else(event == "Last follow-up", time, 0),
    days_to_last_follow_up = if_else(event == "Disease progression", time, 0)
  ) %>%
  mutate_at(c("days_to_death", "days_to_last_follow_up"),
            manuscriptsR::na_x2na) %>%
  mutate(sample_type = "PrimaryTumor") %>%
  select(c(
    "sample_type",
    "sample",
    "days_to_death",
    "days_to_last_follow_up"
  ))

dim(clean_clin)
# clean expression ----------------------------------------------------------
dat %>% glimpse()
dat %>% head()

clean_expression <- 
  dat %>% 
  select(clean_clin$sample) %>% 
  as.matrix() 
dim(clean_expression)

# clean_expression <- head(clean_expression)
# rownames(clean_expression) <- rownames(rnaExpr)
# 生存分析 -------------------------------------------------------

library(GDCRNATools)
# Last follow-up  ➡  event
# Disease progression  ➡  alive
#  mutate(sample_type = "PrimaryTumor") 
# sample_type 一定要是 PrimaryTumor

survOutput <-
  gdcSurvivalAnalysis(
    gene = rownames(clean_expression) %>% head() %>% c("LXN"),
    rna.expr = clean_expression,
    metadata = clean_clin
) %>% 
  select(-symbol) %>% 
  rownames_to_column("symbol")

gse_id <- "GSE102349"
data_path  <- "analysis/data/derived_data/{gse_id}" %>% stringr::str_glue()
fs::dir_create(data_path)

rio::export(clin, str_glue("{data_path}/{gse_id}_clincal_info.csv"))
dat %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  select(gene_id, everything()) %>%
  rio::export(., str_glue("{data_path}/{gse_id}_expresioon_matrix.csv"))

save.image(str_glue("{data_path}/{gse_id}.Rdata"))

