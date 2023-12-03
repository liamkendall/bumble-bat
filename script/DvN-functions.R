###functions

#basic summary for binary data
#can handle non-integers  (by rounding) 

binary.summary=function (x,n){
  
  x=round(n*x)
  out=c(rep(1,x),rep(0,n-x))
  mean=mean(out)
  sd=sd(out)
  se=sd(out)/sqrt(n)
  
  return(cbind(n,mean,sd,se))
  
}

orchard_sans_pi=function (object, mod = "1", group, xlab, N = NULL, alpha = 0.5, 
                          angle = 90, cb = TRUE, k = TRUE, g = TRUE, tree.order = NULL, 
                          trunk.size = 3, branch.size = 1.2, twig.size = 0.5, transfm = c("none", 
                                                                                          "tanh"), condition.lab = "Condition", legend.pos = c("bottom.right", 
                                                                                                                                               "bottom.left", "top.right", "top.left", "top.out", "bottom.out", 
                                                                                                                                               "none"), k.pos = c("right", "left", "none"), colour = FALSE, 
                          fill = TRUE, weights = "prop", by = NULL, at = NULL, upper = TRUE, 
                          flip = TRUE) 
{
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  if (any(class(object) %in% c("robust.rma", "rma.mv", "rma", 
                               "rma.uni"))) {
    if (mod != "1") {
      results <- orchaRd::mod_results(object, mod, group, 
                                      N, by = by, at = at, weights = weights, upper = upper)
    }
    else {
      results <- orchaRd::mod_results(object, mod = "1", 
                                      group, N, by = by, at = at, weights = weights, 
                                      upper = upper)
    }
  }
  if (any(class(object) %in% c("orchard"))) {
    results <- object
  }
  mod_table <- results$mod_table
  data_trim <- results$data
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, 
                                labels = mod_table$name)
  data_trim$scale <- (1/sqrt(data_trim[, "vi"]))
  legend <- "Precision (1/SE)"
  if (!is.null(tree.order) & length(tree.order) != nlevels(data_trim[, 
                                                                     "moderator"])) {
    stop("Length of 'tree.order' does not equal number of categories in moderator")
  }
  if (!is.null(tree.order)) {
    data_trim$moderator <- factor(data_trim$moderator, levels = tree.order, 
                                  labels = tree.order)
    mod_table <- mod_table %>% dplyr::arrange(factor(name, 
                                                     levels = tree.order))
  }
  if (is.null(N) == FALSE) {
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)")
  }
  if (transfm == "tanh") {
    cols <- sapply(mod_table, is.numeric)
    mod_table[, cols] <- Zr_to_r(mod_table[, cols])
    data_trim$yi <- Zr_to_r(data_trim$yi)
    label <- xlab
  }
  else {
    label <- xlab
  }
  mod_table$K <- as.vector(by(data_trim, data_trim[, "moderator"], 
                              function(x) length(x[, "yi"])))
  mod_table$g <- as.vector(num_studies(data_trim, moderator, 
                                       stdy)[, 2])
  group_no <- length(unique(mod_table[, "name"]))
  cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
            "#AA4499", "#44AA99", "#999933", "#882255", "#661100", 
            "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", 
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  if (colour == TRUE) {
    color <- as.factor(data_trim$stdy)
    color2 <- NULL
  }
  else {
    color <- data_trim$mod
    color2 <- mod_table$name
  }
  if (fill == TRUE) {
    fill <- color
  }
  else {
    fill <- NULL
  }
  if (names(mod_table)[2] == "condition") {
    condition_no <- length(unique(mod_table[, "condition"]))
    plot <- ggplot2::ggplot() + ggbeeswarm::geom_quasirandom(data = data_trim, 
                                                             ggplot2::aes(y = yi, x = moderator, size = scale, 
                                                                          colour = color, fill = fill), alpha = alpha, 
                                                             shape = 21) + ggplot2::geom_hline(yintercept = 0, 
                                                                                               linetype = 2, colour = "black", alpha = alpha) + 
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, 
                                                             ymin = lowerCL, ymax = upperCL), size = branch.size, 
                              position = ggplot2::position_dodge2(width = 0.3)) + 
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, 
                                                              x = name, ymin = lowerCL, ymax = upperCL, shape = as.factor(condition), 
                                                              fill = color2), size = twig.size, position = ggplot2::position_dodge2(width = 0.3), 
                               fatten = trunk.size) + ggplot2::scale_shape_manual(values = 20 + 
                                                                                    (1:condition_no)) + ggplot2::theme_bw() + ggplot2::guides(fill = "none", 
                                                                                                                                              colour = "none") + ggplot2::theme(legend.position = c(0, 
                                                                                                                                                                                                    1), legend.justification = c(0, 1)) + ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) + 
      ggplot2::theme(legend.direction = "horizontal") + 
      ggplot2::theme(legend.background = ggplot2::element_blank()) + 
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) + 
      ggplot2::labs(shape = condition.lab) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, 
                                                                                                colour = "black", hjust = 0.5, angle = angle))
    if (flip) {
      plot <- plot + ggplot2::coord_flip()
    }
  }
  else {
    plot <- ggplot2::ggplot() + ggbeeswarm::geom_quasirandom(data = data_trim, 
                                                             ggplot2::aes(y = yi, x = moderator, size = scale, 
                                                                          colour = color, fill = fill), alpha = alpha, 
                                                             shape = 21) + ggplot2::geom_hline(yintercept = 0, 
                                                                                               linetype = 2, colour = "black", alpha = alpha) + 
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, 
                                                             ymin = lowerCL, ymax = upperCL), size = branch.size) + 
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, 
                                                              x = name, ymin = lowerCL, ymax = upperCL, fill = color2), 
                               size = twig.size, fatten = trunk.size, shape = 21) + 
      ggplot2::theme_bw() + ggplot2::guides(fill = "none", 
                                            colour = "none") + ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) + 
      ggplot2::theme(legend.direction = "horizontal") + 
      ggplot2::theme(legend.background = ggplot2::element_blank()) + 
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) + 
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, 
                                                         colour = "black", hjust = 0.5, angle = angle))
    if (flip) {
      plot <- plot + ggplot2::coord_flip()
    }
  }
  if (legend.pos == "bottom.right") {
    plot <- plot + ggplot2::theme(legend.position = c(1, 
                                                      0), legend.justification = c(1, 0))
  }
  else if (legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position = c(0, 
                                                      0), legend.justification = c(0, 0))
  }
  else if (legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position = c(1, 
                                                      1), legend.justification = c(1, 1))
  }
  else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position = c(0, 
                                                      1), legend.justification = c(0, 1))
  }
  else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position = "top")
  }
  else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position = "bottom")
  }
  else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }
  if (cb == TRUE) {
    plot <- plot + ggplot2::scale_fill_manual(values = cbpl) + 
      ggplot2::scale_colour_manual(values = cbpl)
  }
  if (k == TRUE && g == FALSE && k.pos == "right") {
    plot <- plot + ggplot2::annotate("text", y = (max(data_trim$yi) + 
                                                    (max(data_trim$yi) * 0.1)), x = (seq(1, group_no, 
                                                                                         1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no]), 
                                     parse = TRUE, hjust = "right", size = 3.5)
  }
  else if (k == TRUE && g == FALSE && k.pos == "left") {
    plot <- plot + ggplot2::annotate("text", y = (min(data_trim$yi) + 
                                                    (min(data_trim$yi) * 0.1)), x = (seq(1, group_no, 
                                                                                         1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no]), 
                                     parse = TRUE, hjust = "left", size = 3.5)
  }
  else if (k == TRUE && g == TRUE && k.pos == "right") {
    plot <- plot + ggplot2::annotate("text", y = (max(data_trim$yi) + 
                                                    (max(data_trim$yi) * 0.1)), x = (seq(1, group_no, 
                                                                                         1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no], 
                                                                                                                  "~", "(", mod_table$g[1:group_no], ")"), parse = TRUE, 
                                     hjust = "right", size = 3.5)
  }
  else if (k == TRUE && g == TRUE && k.pos == "left") {
    plot <- plot + ggplot2::annotate("text", y = (min(data_trim$yi) + 
                                                    (min(data_trim$yi) * 0.1)), x = (seq(1, group_no, 
                                                                                         1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no], 
                                                                                                                  "~", "(", mod_table$g[1:group_no], ")"), parse = TRUE, 
                                     hjust = "left", size = 3.5)
  }
  else if (k == TRUE && g == FALSE && k.pos %in% c("right", 
                                                   "left", "none") == FALSE) {
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, 
                                                                 group_no, 1) + 0.3), label = paste("italic(k)==", 
                                                                                                    mod_table$K[1:group_no]), parse = TRUE, size = 3.5)
  }
  else if (k == TRUE && g == TRUE && k.pos %in% c("right", 
                                                  "left", "none") == FALSE) {
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, 
                                                                 group_no, 1) + 0.3), label = paste("italic(k)==", 
                                                                                                    mod_table$K[1:group_no], "~", "(", mod_table$g[1:group_no], 
                                                                                                    ")"), parse = TRUE, size = 3.5)
  }
  return(plot)
}

