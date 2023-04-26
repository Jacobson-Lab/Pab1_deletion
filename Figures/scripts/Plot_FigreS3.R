# --------------------------------------------------
# Figures S3
# --------------------------------------------------
library(dplyr)
library(ggplot2)
set3 <- RColorBrewer::brewer.pal(n = 12, name = "Set3") # Color palette for plotting

# Random forest performance metrics ----------------
nrmse <- read.table("../data/FigS3A_pab1_rf_reg_5foldCV_5repeats_nrmse.txt", header = TRUE, sep = "\t")
nrmse$sample <- recode_factor(nrmse$sample, WT = "WT", pbp1d = "pbp1Δ", pab1d = "pab1Δ\npbp1Δ")
nrmse_testing <- nrmse[which(nrmse$dataset == "testing"), ]
p_nrmse <- ggplot(nrmse_testing, aes(x = sample, y = nrmse_range)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = sample)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.9), size = 0.2) +
  scale_fill_manual(values = set3[c(4, 7, 5)]) +
  ylab("NRMSE") + xlab("") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), strip.background = element_blank(), axis.text.x = element_text(face = "italic"),
        legend.position = "none")

auroc <- read.table("../data/FigS3B_pab1_rf_class_5foldCV_5repeats_auroc.txt", header = TRUE, sep = "\t")
auroc$sample <- recode_factor(auroc$sample, WT = "WT", pbp1d = "pbp1Δ", pab1d = "pab1Δ\npbp1Δ")
auroc_testing <- auroc[which(auroc$dataset == "testing"), ]
p_auroc <- ggplot(auroc_testing, aes(x = sample, y = auc)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = sample)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.9), size = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(values = set3[c(4, 7, 5)]) +
  ylab("AUROC") + xlab("") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), strip.background = element_blank(), axis.text.x = element_text(face = "italic"),
        legend.position = "none")

# Stop codon and surrounding nucleotides -----------
dfs <- read.table("../data/FigS3C_pab1_umi_f0_0_reg_ntm6-stop-ntp9_padj_BH.txt", header = TRUE, sep = "\t")
dfs$Sample <- recode_factor(dfs$Sample, WT = "WT", pbp1d = "pbp1Δ", pab1d = "pab1Δpbp1Δ")
dfs$feature <- recode_factor(dfs$feature, random_factor = "Random", nt_m06 = "-6",  nt_m05 = "-5", nt_m04 = "-4", 
                             nt_m03 = "-3", nt_m02 = "-2", nt_m01 = "-1", stop_codon = "Stop", 
                             nt_p04 = "+4", nt_p05 = "+5", nt_p06 = "+6", nt_p07 = "+7", nt_p08 = "+8", nt_p09 = "+9")
p_stop <- ggplot(dfs) +
  geom_tile(aes(x = Var, y = Sample, fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~feature, scales = "free", space = "free_x") +
  scale_fill_gradient2(name = "Group's median\nrelative to sample median     ",
                       low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-1, 1), oob = scales::squish,
                       breaks = seq(-1, 1, length = 5), labels = c("< -1.0", "-0.5\nLower", "0.0", "0.5\nHigher", "> 1.0")) +
  scale_size_manual(name = "Significance  ", values = c(`p < 0.05` = 0, ns = 3), limits = c("p < 0.05", "ns")) +
  scale_y_discrete(limits = rev(levels(dfs$Sample))) +
  xlab("") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(face = "italic"), axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# P-site codon -------------------------------------
dfpc <- read.table("../data/FigS3D_pab1_umi_f0_0_reg_codon_m01_padj_BH.txt", header = TRUE, sep = "\t")
dfpc$Sample <- recode_factor(dfpc$Sample, WT = "WT", pbp1d = "pbp1Δ", pab1d = "pab1Δpbp1Δ")
dfpc$aa <- factor(dfpc$aa, levels = c("F", "S", "Y", "*", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"))
p_psite <- ggplot(dfpc) +
  geom_tile(aes(x = Var, y = Sample, fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~aa, scales = "free", space = "free_x") +
  scale_fill_gradient2(name = "Group's median\nrelative to sample median     ",
                       low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-3, 3), oob = scales::squish,
                       breaks = seq(-3, 3, length = 5), labels = c("< -3.0", "-1.5\nLower", "0.0", "1.5\nHigher", "> 3.0")) +
  scale_size_manual(name = "Significance  ", values = c(`p < 0.05` = 0, ns = 3), limits = c("p < 0.05", "ns")) +
  scale_y_discrete(limits = rev(levels(dfs$Sample))) +
  xlab("P-site codon") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(face = "italic"), axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# Combine plots
library(patchwork)
design = "
AB#
CCC
DDD
"
pA <- p_nrmse + theme(axis.title.y = element_text(margin = margin(r = -45, unit = "pt")))
p <- pA + p_auroc + p_stop + p_psite +
  plot_layout(design = design, heights = c(4.5, 1, 1)) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft")

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "FigS3.pdf", family = "Arial", width = 7.5, height = 6.5) 
p
dev.off()

