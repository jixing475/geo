#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gse PARAM_DESCRIPTION
#' @param db PARAM_DESCRIPTION, Default: 'KEGG_2019_Human'
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
#'  \code{\link[enrichR]{enrichr}},\code{\link[enrichR]{enrichR}}
#'  \code{\link[stringr]{str_glue}}
#' @rdname enrich_db_analysis
#' @export 
#' @importFrom enrichR enrichr
#' @importFrom stringr str_glue
enrich_db_analysis <-
  function(gse, db = "KEGG_2019_Human", trend = "Up") {
    # 设置要分析的分子, 上 or 下
    if (trend == "Up") {
      genes <- gse$up_genes %>% toupper()
    } else {
      genes <- gse$down_genes %>% toupper()
    }
    # 富集分析
    enrichr_res <- enrichR::enrichr(genes, db) #翻墙会快很多
    names(enrichr_res) <- paste0(db, "_table")
    
    enrichr_res[[paste0(db, "_plot")]] <- 
      enrichr_res[[1]] %>%
      plot_enrich(title = stringr::str_glue("{db}({trend})"))
    return(enrichr_res)
  }

#test <- all_res$GSE12452 %>% enrich_db_analysis()
