# Functions for determining significant fold and plotting #
# ----------------------------------------------------------------
# Function to add columns of significance based a cutoff, direction of changes, axis limits for plotting
add_sig_fc <- function(df, pval_col = "padj", pval_cutoff = 0.01, log2FC_col = "log2FoldChange", log2FC_cutoff = 0, 
                       plot_fc_cutoff = 5, plot_pval_cutoff = 10^-50) {
  # Significance column
  df$sig <- "ns"
  df[which(df[, pval_col] < pval_cutoff), ]$sig <- paste("p < ", pval_cutoff)
  # Fold change column
  df$changes <- "Unchanged"
  df[which(df$sig != "ns" & df[, log2FC_col] > log2FC_cutoff), ]$changes <- "Up"
  df[which(df$sig != "ns" & df[, log2FC_col] < -log2FC_cutoff), ]$changes <- "Down"
  # For plotting points that are out-of-bounds
  df$log2FC_lim <- df[, log2FC_col]
  tryCatch({df[which(df[, log2FC_col] > plot_fc_cutoff), ]$log2FC_lim <- plot_fc_cutoff}, error = function(e){})
  tryCatch({df[which(df[, log2FC_col] < -plot_fc_cutoff), ]$log2FC_lim <- -plot_fc_cutoff}, error = function(e){})
  df$pval_lim <- df[, pval_col]
  tryCatch({df[which(is.na(df$pval_lim)), ]$pval_lim <- 1}, error = function(e){})
  tryCatch({df[which(df[, pval_col] < plot_pval_cutoff), ]$pval_lim <- plot_pval_cutoff}, error = function(e){})
  df$plot_lim <- "in"
  tryCatch({df[which(df$pval_lim == plot_pval_cutoff), ]$plot_lim <- "top"}, error = function(e){})
  tryCatch({df[which(df$pval_lim > plot_pval_cutoff & abs(df$log2FC_lim) == plot_fc_cutoff), ]$plot_lim <- "side"}, error = function(e){})
  tryCatch({df[which(df$pval_lim == plot_pval_cutoff & abs(df$log2FC_lim) == plot_fc_cutoff), ]$plot_lim <- "both"}, error = function(e){})
  return(df)
}
# ----------------------------------------------------------------
# Function to make volcano plot
plot_volcano <- function(res, nn, xlabel = NULL, point_label = NULL) {
  require(ggplot2)
  require(ggrepel)
  res$repel <- ""
  tryCatch({res[which(res$Protein %in% point_label), ]$repel <- res[which(res$Protein %in% point_label), ]$Protein}, error = function(e){})
  p <- ggplot(res, aes(x = log2FC_lim, y = -log10(pval_lim), color = changes, shape = plot_lim)) +
    geom_point(alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_text(data = nn, aes(x = min(res$log2FC_lim), y = ypos, label = paste0("n = ", Freq), color = changes), 
              size = 7/.pt, hjust = 0, show.legend = FALSE, inherit.aes = FALSE) +
    geom_label_repel(aes(label = repel, color = changes), size = 7/.pt, 
                     max.overlaps = Inf, label.padding = 0.1, min.segment.length = 0, show.legend = FALSE) +
    facet_wrap(~pair, nrow = 1) +
    scale_color_manual(name = "Significant change: ", values = c(Up = "orange", Down = "purple", Unchanged = "grey50")) +
    scale_shape_manual(name = "", values = c(`in` = 1, top = 2, side = 0, both = 14)) +
    xlab(bquote(log[2] ~ "fold change" ~ .(xlabel))) + ylab(expression("-log"[10]*"(adjusted p-value)")) +
    theme_bw(base_size = 10) + 
    theme(strip.background = element_rect(fill = "white"), panel.grid = element_blank(),
          legend.position = "top", strip.text = element_text(face = "italic")) +
    guides(color = guide_legend(override.aes = list(alpha = 1)), shape = "none")
  return(p)
}
# ----------------------------------------------------------------
# Function to make scatter plot of mRNA vs protein log2 fold changes
plot_scatter <- function(combined_filter, facet_group = "pair", nnc) {
  p <- ggplot(combined_filter, aes(x = log2FoldChange, y = logFC)) +
    geom_point(aes(color = sig_both), alpha = 0.25) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    geom_label_repel(aes(label = repel, color = sig_both), size = 7/.pt, 
                     max.overlaps = Inf, label.padding = 0.1, min.segment.length = 0, show.legend = FALSE) +
    scale_color_manual(name = "mRNA vs. protein: ", 
                       values = c(Neither = "grey50", `Protein only` = "blue", `mRNA only` = "red", Both = "forestgreen")) +
    xlab(expression("log"[2]*" fold change in mRNA abundance")) + ylab(expression("log"[2]*" fold change in protein abundance")) +
    theme_bw(base_size = 10) + 
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), strip.text = element_text(face = "italic"),
          legend.position = "top") +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  if (facet_group == "pair") {
    p <- p + facet_wrap(~pair, nrow = 1) + 
      geom_text(data = nnc, aes(x = -12.5, y = ypos, label = paste0("n = ", Freq), color = sig_both), 
                size = 7/.pt, hjust = 0, show.legend = FALSE, inherit.aes = FALSE) +
      stat_cor(method = "spearman", size = 7/.pt, label.x = -12.5, label.y = 1.8, cor.coef.name = "rho") +
      coord_cartesian(xlim = c(-12, 5), ylim = c(-2, 2))
  }
  if (facet_group == "prot_group") {
    p <- p + facet_wrap(~prot_group)
  }
  return(p)
}
# ----------------------------------------------------------------