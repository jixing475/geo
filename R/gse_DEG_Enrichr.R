#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gse_id PARAM_DESCRIPTION
#' @param pvalue PARAM_DESCRIPTION, Default: 'P.Value'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  gse <- gse_DEG_Enrichr("GSE12452", pvalue = "P.Value")
#   class(gse)
#'  }
#' }
#' @seealso 
#'  \code{\link[readr]{read_delim}}
#'  \code{\link[stringr]{str_glue}},\code{\link[stringr]{str_remove}}
#'  \code{\link[dplyr]{filter}},\code{\link[dplyr]{pull}}
#'  \code{\link[rlang]{parse_expr}}
#'  \code{\link[purrr]{set_names}}
#' @rdname gse_DEG_Enrichr
#' @export 
#' @importFrom readr read_csv
#' @importFrom stringr str_glue str_remove
#' @importFrom dplyr filter pull
#' @importFrom rlang parse_expr
#' @importFrom purrr set_names
## 修改获取基因 list的标准 P.Value or adj.P.Val(default)
# P.Value or  adj.P.Val
# input: gse id and 
# output: list result

gse_DEG_Enrichr <- function(gse_id, pvalue = "P.Value") {
  gse <- list()
  gse$id <- gse_id
  
  gse$expression_matrix_table <-
    readr::read_csv(
      stringr::str_glue(
        "analysis/data/raw_data/{gse$id}/{gse$id}_Deg_data.csv"
      )
    )
  
  gse$DEG_table <-
    readr::read_csv(
      stringr::str_glue(
        "analysis/data/raw_data/{gse$id}/{gse$id}_Deg_result.csv"
      )
    )
  
  ## 获取 gene list
  gse$up_genes <-
    gse %>%
    .[["DEG_table"]] %>%
    dplyr::filter(logFC >= 1) %>%
    dplyr::filter(!!rlang::parse_expr(pvalue)<= 0.05) %>%
    dplyr::pull(X1) 
  
  gse$down_genes <-
    gse %>%
    .[["DEG_table"]] %>%
    dplyr::filter(logFC <= -1) %>%
    dplyr::filter(!!rlang::parse_expr(pvalue)<= 0.05) %>%
    dplyr::pull(X1) 
  ## 火山图
  gse$volcano_plot <- plot_Volvano(gse = gse)
  ## 富集分析: Up
  gse$enrichr_up <-
    enrich_KEGG_GO_analysis(gse = gse, trend = "Up") %>%
    purrr::set_names( ~ stringr::str_remove(.x, "_2018|_2019"))
  gse$enrichr_down <-
    enrich_KEGG_GO_analysis(gse = gse, trend = "Down") %>% 
    purrr::set_names( ~ stringr::str_remove(.x, "_2018|_2019"))
  return(gse)
}

