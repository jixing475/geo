#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname load_pkg
#' @export 

load_pkg <- function(){
  library(ggtext)
  library(GDCRNATools)
  library(enrichR)
  library(janitor)
  # survival analysis
  library("survival")
  library("survminer")
  library("gtsummary")
}
