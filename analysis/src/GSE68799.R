
library(tidyverse)
gse_id <- "GSE68799"
data_path  <- "analysis/data/derived_data/{gse_id}" %>% stringr::str_glue()
fs::dir_create(data_path)
# load data  -------------------------------------------------------------
gse_68799 <- rio::import("analysis/data/raw_data/GSE68799_NPC_FPKM.txt.gz") 

library(rvest)
library(curl)
library(tidyverse)
library(robotstxt)
url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68799"
page <- read_html(curl(url, handle = curl::new_handle("useragent" = "Mozilla/5.0")))
#====🔥find where is the table lives on this webpage====
table_path='tr:nth-child(23) tr td'
meta <- 
page %>%
  html_nodes(table_path) %>%
  html_text() %>% 
  matrix(ncol = 2, byrow = TRUE) %>% 
  as.data.frame() %>% 
  set_names(c("sample_id", "group")) %>% 
  separate(col = group, into = c("group", "sample"), sep = "-") %>% 
  mutate(sample = str_c("X", sample))
meta %>% head()

gse_68799_matrix <- 
gse_68799 %>% 
  as_tibble() %>% 
  column_to_rownames("GeneID")
# DEG -------------------------------------------------------------
library(GDCRNATools)

DEGAll <- 
  gdcDEAnalysis(
  counts = gse_68799_matrix,
  group = meta$group,
  comparison = 'NPC-nonNPC',
  #比较组信息
  method = 'limma'
) %>% 
  rownames_to_column("gene_id") %>% 
  select(gene_id, everything())

write.table(DEGAll, file=stringr::str_glue("{data_path}/{gse_id}_DEG_all_gene.csv"), row.names=F, sep=",")

gse_68799 %>% 
  rio::export(., str_glue("{data_path}/{gse_id}_expression_matrix.csv"))


test <- 
gse_68799_matrix %>% 
  filter(rownames(.) == "LXN_chr3") %>% 
  manuscriptsR::df_turn() %>% 
  rownames_to_column("sample") %>% 
  left_join(meta) 

ggplot(test) +
  aes(x = "", y = LXN_chr3, fill = group) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_hue() +
  theme_minimal()

# save image -------------------------------------------------------------
save.image(str_glue("{data_path}/{gse_id}.Rdata"))
