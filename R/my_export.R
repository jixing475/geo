
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param res PARAM_DESCRIPTION
#' @param path PARAM_DESCRIPTION, Default: '~/Desktop'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[manuscriptsR]{my_word}}
#'  \code{\link[stringr]{str_glue}}
#'  \code{\link[export]{graph2office}},\code{\link[export]{character(0)}}
#'  \code{\link[rio]{export}}
#'  \code{\link[readr]{read_lines}}
#'  \code{\link[fs]{create}}
#' @rdname my_export
#' @export 
#' @importFrom manuscriptsR my_word
#' @importFrom stringr str_glue
#' @importFrom export graph2ppt export
#' @importFrom rio export
#' @importFrom readr write_lines
#' @importFrom fs dir_create
#' 

my_export <- function(res, path = "~/Desktop"){
  for(it in names(res)){
    # word
    if("flextable" %in% class(res[[it]])){
      print(it)
      manuscriptsR::my_word(res[[it]],
                            para = it,
                            open = FALSE,
                            docx_path = stringr::str_glue("{path}/Tables.docx"))
    } 
    # ppt
    else if(
      intersect(c("ggplot", "ggsurv", "ggsurvplot"), class(res[[it]])) %>% length() %>% `>`(0)
    ){
      print(it)
      print(res[[it]])
      Sys.sleep(1)
      export::graph2ppt(
        file = paste0(path, "/Figures.ppt"),
        append = TRUE,
        width=9, aspectr=1.5
      )
    }
    # csv
    else if("data.frame" %in% class(res[[it]])){
      print(it)
      rio::export(res[[it]], file = stringr::str_glue("{path}/{it}.csv"))
    }
    # character to txt
    else if("character" %in% class(res[[it]])){
      print(it)
      readr::write_lines(res[[it]], file = stringr::str_glue("{path}/{it}.txt"))
    }
    # list recursive
    else if("list" == class(res[[it]])){
      print(it)
      sub_path <- stringr::str_glue("{path}/{it}")
      fs::dir_create(sub_path)
      my_export(res[[it]], path = sub_path)
    }
  }
}
