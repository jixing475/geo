#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gse_id PARAM_DESCRIPTION
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
#'  \code{\link[dplyr]{mutate}},\code{\link[dplyr]{select}}
#'  \code{\link[purrr]{set_names}}
#'  \code{\link[ggplot2]{labs}}
#' @rdname plot_Volvano
#' @export 
#' @importFrom stringr str_glue
#' @importFrom readr read_csv
#' @importFrom dplyr mutate select
#' @importFrom purrr set_names
#' @importFrom ggplot2 ylab
plot_Volvano <- function(gse){
  deg_res <-
    gse$DEG_table %>% 
    dplyr::mutate(group = "protein_coding") %>%
    dplyr::select(c("X1", "group", "logFC", "P.Value")) %>%
    purrr::set_names(c("symbol", "group", "logFC", "FDR"))
  
  p <- my_gdcVolcanoPlot(deg.all = deg_res, pval = 0.05) +
    labs(x = "Log<sub>2</sub>(fold change)",
         y = "-Log<sub>10</sub>P−value") +
    theme(
      axis.title.x = element_markdown(),
      axis.title.y = element_markdown()
    )
  return(p)
}
