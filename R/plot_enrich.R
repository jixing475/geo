#' @title plot enrichment analysis
#' @description FUNCTION_DESCRIPTION
#' @param enrich_df PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{arrange}},\code{\link[dplyr]{mutate}},\code{\link[dplyr]{rename}},\code{\link[dplyr]{select}}
#'  \code{\link[janitor]{clean_names}}
#'  \code{\link[utils]{head}}
#'  \code{\link[tidyr]{separate}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{scale_colour_gradient}},\code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{scale_x_discrete}},\code{\link[ggplot2]{theme}}
#'  \code{\link[forcats]{fct_reorder}}
#' @rdname plot_enrich
#' @export 
#' @importFrom dplyr arrange mutate rename select
#' @importFrom janitor clean_names
#' @importFrom utils head
#' @importFrom tidyr separate
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient theme_bw labs scale_y_discrete theme
#' @importFrom forcats fct_reorder2
#' 
plot_enrich <- function(enrich_df, title = "Statistics of Enrichment") {
  data4plot <-
    enrich_df %>%
    dplyr::arrange(Adjusted.P.value) %>% 
    janitor::clean_names() %>%
    utils::head(10) %>%
    tidyr::separate(overlap,
                    into = c("input_number", "background_number"),
                    sep = "/") %>%
    dplyr::mutate(rich_factor = as.numeric(input_number) / as.numeric(background_number)) %>%
    dplyr::mutate(number_term = as.factor(term)) %>%
    dplyr::mutate(PValue = adjusted_p_value, `Count` = as.integer(input_number)) %>%
    dplyr::mutate(PValue = -log10(PValue)) %>% 
    dplyr::select(
      c(
        "number_term",
        "rich_factor",
        "Count",
        "background_number",
        "PValue",
        "genes"
      )
    )
  # set variable
  low_p  <- min(data4plot$PValue)
  high_p <- max(data4plot$PValue)
  
  p <-
    data4plot %>%
    ggplot2::ggplot(ggplot2::aes(
      x = rich_factor,
      y = forcats::fct_reorder2(number_term, `Count`, rich_factor)
    )) +
    ggplot2::geom_point(ggplot2::aes(size = `Count`, colour = PValue)) +
    ggplot2::scale_size(range = c(2, 9)) +
    ggplot2::scale_color_gradient(
      "-log10(p.adjust)",
      limits = c(low_p, high_p),
                                  low = "#417BB4",
                                  high = "#AF5546") +
    ggplot2::theme_bw(base_size = 14, base_family='Times') +
    ggplot2::labs(title = title, x = "Enrichment Ratio", y = NULL) +
    ggplot2::scale_y_discrete(labels = get_wraper(40)) +
    # when term is 20, use 3/1
    ggplot2::theme(aspect.ratio=4/3) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}
