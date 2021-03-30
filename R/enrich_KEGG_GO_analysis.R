#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gse_id PARAM_DESCRIPTION
#' @param trend PARAM_DESCRIPTION, Default: 'Up'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stringr]{str_glue}}
#'  \code{\link[readr]{read_delim}}
#'  \code{\link[dplyr]{filter}},\code{\link[dplyr]{pull}}
#'  \code{\link[enrichR]{enrichr}},\code{\link[enrichR]{enrichR}}
#' @rdname enrich_KEGG_GO_analysis
#' @export 
#' @importFrom stringr str_glue
#' @importFrom readr read_csv
#' @importFrom dplyr filter pull
#' @importFrom enrichR enrichr
enrich_KEGG_GO_analysis <- function(gse, trend = "Up") {
  #数据库设置
  #install.packages("enrichR")
  #library(enrichR)
  #dbs <- listEnrichrDbs() #列出168个库
  #dbs[1:4,1:4] #从中选择你要富集的库
  #dbs$libraryName #查看库名
  dbs <-c("KEGG_2019_Human", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")
  #
  # deg_genes_path <-
  #   stringr::str_glue("analysis/data/derived_data/{gse_id}/{gse_id}_Deg_Genes.csv")
  # 设置要分析的分子, 上 or 下
  if(trend == "Up"){
    genes <- gse$up_genes %>% toupper()
  } else {
    genes <- gse$down_genes %>% toupper()
  }
  # 富集分析
  enrichr_res <- enrichR::enrichr(genes, dbs) #翻墙会快很多
  
  p_kegg <-
    enrichr_res$KEGG_2019_Human %>%
    plot_enrich(title = str_glue("KEGG pathway({trend})"))
  p_go_BP <-
    enrichr_res$GO_Biological_Process_2018 %>%
    plot_enrich(title = str_glue("GO Biological Process({trend})"))
  p_go_CC <-
    enrichr_res$GO_Cellular_Component_2018 %>%
    plot_enrich(title = str_glue("GO Cellular Component({trend})"))
  p_go_MF <-
    enrichr_res$GO_Molecular_Function_2018 %>%
    plot_enrich(title = str_glue("GO Molecular Function({trend})"))
  res <- list(enrichr_res, p_kegg, p_go_BP, p_go_CC, p_go_MF)
  names(res) <- c("enrichr_res", "p_kegg", "p_go_BP", "p_go_CC", "p_go_MF")
  return(res)
}
