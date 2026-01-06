plot_tornado <- function(dsa_results, base_case_icer, wtp, param_labels) {
  
  library(ggplot2)
  
  # Merge the labels with dsa_results
  dsa_results <- merge(dsa_results, param_labels, by = "parameter", all.x = TRUE)
  dsa_results$display_label <- dsa_results$description 
  
  # Calculate parameter ranges and order by impact
  param_ranges <- aggregate(icer ~ display_label, data = dsa_results, 
                            function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  colnames(param_ranges)[2] <- "range"
  param_ranges <- param_ranges[order(-param_ranges$range), ]
  
  # Merge and set factor levels
  plot_data <- merge(dsa_results, param_ranges, by = "display_label")
  plot_data$display_label <- factor(plot_data$display_label, 
                                    levels = rev(param_ranges$display_label))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = icer, y = display_label)) +
    geom_line(aes(group = display_label), linewidth = 1.5, color = "steelblue") +
    geom_point(aes(shape = bound), size = 3, color = "steelblue") +
    geom_vline(xintercept = base_case_icer, linetype = "dashed", color = "black") +
    geom_vline(xintercept = wtp, linetype = "dotted", color = "red") +
    scale_shape_manual(values = c("lower" = 16, "upper" = 17)) +
    labs(
      title = "Tornado Diagram: Deterministic Sensitivity Analysis",
      subtitle = paste0("Base case ICER: £", format(round(base_case_icer), big.mark = ",")),
      x = "ICER (£/QALY)",
      y = NULL
    ) +
    theme_minimal()
  
  return(p)
}