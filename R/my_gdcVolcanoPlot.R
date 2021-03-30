#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param deg.all PARAM_DESCRIPTION
#' @param fc PARAM_DESCRIPTION, Default: 2
#' @param pval PARAM_DESCRIPTION, Default: 0.01
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes}},\code{\link[ggplot2]{geom_point}},\code{\link[ggplot2]{labs}},\code{\link[ggplot2]{lims}},\code{\link[ggplot2]{geom_abline}},\code{\link[ggplot2]{ggtheme}},\code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}}
#' @rdname my_gdcVolcanoPlot
#' @export 
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab xlim geom_vline geom_hline theme_bw theme element_line element_blank element_rect element_text
my_gdcVolcanoPlot <-
  function (deg.all, fc = 2, pval = 0.05) {
    geneList <- deg.all
    geneList$threshold <- c()
    geneList$threshold[geneList$logFC > log(fc, 2) & geneList$FDR <
                         pval] <- "Up−regulation"
    geneList$threshold[geneList$logFC >= -log(fc, 2) &
                         geneList$logFC <=
                         log(fc, 2) | geneList$FDR >= pval] <- "None"
    geneList$threshold[geneList$logFC < -log(fc, 2) & geneList$FDR <
                         pval] <- "Down−regulation"
    geneList$threshold <- as.factor(geneList$threshold)
    lim <- max(max(geneList$logFC), abs(min(geneList$logFC))) +
      0.5
    ggplot2::ggplot(data = geneList, ggplot2::aes(x = logFC, y = -log10(FDR))) +
      ggplot2::geom_point(aes(color = threshold), 
                          alpha = 0.4, size = 2) +
      ggplot2::scale_color_manual(name = "", 
                                  values = c("#0374B4", "#BEBEBE", "#BB3C28"))+
      # ggplot2::xlab("log2(Fold Change)") +
      # ggplot2::ylab("-log10(FDR)") +
      ggplot2::xlim(c(-lim, lim)) +
      ggplot2::geom_vline(
        xintercept = c(-log(fc,2), log(fc, 2)),
        color = "darkgreen",
        linetype = 3
      ) +
      ggplot2::geom_hline(
        yintercept = -log(pval, 10),
        color = "darkgreen",
        linetype = 3
      ) +
      ggplot2::theme_bw(base_family='Times') +
      ggplot2::theme(
        axis.line = ggplot2::element_line(colour = "black"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "black"),
        panel.background = ggplot2::element_blank(),
        legend.position = c(0.1, 0.9) 
      ) +
       ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                      axis.title = ggplot2::element_text(size = 16))
    
  }


