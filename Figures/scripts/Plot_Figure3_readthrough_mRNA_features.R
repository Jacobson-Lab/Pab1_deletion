# --------------------------------------------------
# Figures 3
# --------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggh4x)

# Random forest feature importance -----------------
df_reg <- read.table("../data/Fig3A_pab1_rf_reg_5foldCV_5repeats_feature_importance.txt", sep = "\t", header = TRUE)

  # Function to rename features and prepare data for plotting
prep_imp <- function(all.imp_df) {
  # Recode sample
  all.imp_df$sample <- recode_factor(all.imp_df$sample, WT = "WT", pbp1d = "pbp1Δ", pab1d = "pab1Δpbp1Δ")
  # Rename features for x-axis tick
  all.imp_df$xtick <- as.factor(all.imp_df$feature)
  all.imp_df$xtick <- recode_factor(all.imp_df$xtick,
                                    tunnel_lower_H_count = "Helical", tunnel_lower_neg_charge = "- charged", tunnel_lower_no_charge = "No charge", tunnel_lower_pos_charge = "+ charged", tunnel_lower_aromatic = "Aromatic", tunnel_lower_polar = "Polar", tunnel_lower_nonpolar = "Nonpolar", tunnel_lower_hydrophilic = "Hydrophylic", tunnel_lower_neutral = "Neutral", tunnel_lower_hydrophobic = "Hydrophobic", tunnel_lower_v_hydrophobic = "Very hydrophobic", 
                                    tunnel_central_neg_charge = " - charged", tunnel_central_no_charge = " No charge", tunnel_central_pos_charge = " + charged", tunnel_central_aromatic = " Aromatic", tunnel_central_polar = " Polar", tunnel_central_nonpolar = " Nonpolar", tunnel_central_hydrophilic = " Hydrophylic", tunnel_central_neutral = " Neutral", tunnel_central_hydrophobic = " Hydrophobic", tunnel_central_v_hydrophobic = " Very hydrophobic",
                                    tunnel_constriction_neg_charge = "  - charged", tunnel_constriction_no_charge = "  No charge", tunnel_constriction_pos_charge = "  + charged", tunnel_constriction_aromatic = "  Aromatic", tunnel_constriction_polar = "  Polar", tunnel_constriction_nonpolar = "  Nonpolar", tunnel_constriction_hydrophilic = "  Hydrophylic", tunnel_constriction_neutral = "  Neutral", tunnel_constriction_hydrophobic = "  Hydrophobic", tunnel_constriction_v_hydrophobic = "  Very hydrophobic",
                                    tunnel_upper_neg_charge = "   - charged", tunnel_upper_no_charge = "   No charge", tunnel_upper_pos_charge = "   + charged", tunnel_upper_aromatic = "   Aromatic", tunnel_upper_polar = "   Polar", tunnel_upper_nonpolar = "   Nonpolar", tunnel_upper_hydrophilic = "   Hydrophylic", tunnel_upper_neutral = "   Neutral", tunnel_upper_hydrophobic = "   Hydrophobic", tunnel_upper_v_hydrophobic = "   Very hydrophobic",
                                    aa_m02 = "E-site aa", aa_m01 = "P-site aa", nt_m18 = "-18", nt_m17 = "-17", nt_m16 = "-16", nt_m15 = "-15", nt_m14 = "-14", nt_m13 = "-13", nt_m12 = "-12", nt_m11 = "-11", nt_m10 = "-10", nt_m09 = "-9", nt_m08 = "-8", nt_m07 = "-7", nt_m06 = "-6", nt_m05 = "-5", nt_m04 = "-4", nt_m03 = "-3", nt_m02 = "-2", nt_m01 = "-1", stop_codon = "Stop codon", nt_p04 = "+4", nt_p05 = "+5", nt_p06 = "+6", nt_p07 = "+7", nt_p08 = "+8", nt_p09 = "+9", nt_p10 = "+10", nt_p11 = "+11", nt_p12 = "+12", nt_p13 = "+13", nt_p14 = "+14", nt_p15 = "+15", nt_p16 = "+16",
                                    nis_stop = "1st 3'UTR stop", l_utr3 = "3'-UTR length", 
                                    utr3_fraction_A = "3'-UTR A fraction", utr3_fraction_C = "3'-UTR C fraction", utr3_fraction_G = "3'-UTR G fraction", utr3_fraction_T = "3'-UTR U fraction",
                                    last_nt_15 = " -15", last_nt_14 = " -14", last_nt_13 = " -13", last_nt_12 = " -12", last_nt_11 = " -11", last_nt_10 = " -10", last_nt_09 = " -9", last_nt_08 = " -8", last_nt_07 = " -7", last_nt_06 = " -6", last_nt_05 = " -5", last_nt_04 = " -4", last_nt_03 = " -3", last_nt_02 = " -2", last_nt_01 = " -1", last_nt_00 = " 0",
                                    random_factor = "Random factor", random_num = "Random number")
  # Group features for faceting
  all.imp_df$xgroup <- NA
  all.imp_df$xgroup <- gsub("tunnel_lower_.*", "aa 20-30\nfrom PTC", all.imp_df$feature)
  all.imp_df$xgroup <- gsub("tunnel_central_.*", "aa 13-19\nfrom PTC", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("tunnel_constriction_.*", "aa 10-12\nfrom PTC", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("tunnel_upper_.*", "aa 3-9\nfrom PTC", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("aa_m.*", "", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("nt_m.*", " nt from stop ", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("stop_codon", " ", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("nt_p.*", "nt from stop", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("l_utr3", "  ", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("utr3_fraction_.*", "  ", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("last_nt.*", "nt from 3' end", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("nis_stop", "   ", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("dist_bp", "2°", all.imp_df$xgroup)
  all.imp_df$xgroup <- gsub("random_.*", "NC", all.imp_df$xgroup)
  all.imp_df$xgroup <- factor(all.imp_df$xgroup, levels = c("aa 20-30\nfrom PTC", "aa 13-19\nfrom PTC", "aa 10-12\nfrom PTC", "aa 3-9\nfrom PTC", "", " nt from stop ", " ", "nt from stop", "  ", "   ", "nt from 3' end", "NC"))
  # Group features for faceting nascent peptide
  all.imp_df$xgroup2 <- NA
  all.imp_df$xgroup2 <- gsub("tunnel.*", "Nascent peptide in the exit tunnel", all.imp_df$feature)
  all.imp_df$xgroup2 <- gsub("aa_m.*", "Nascent peptide in the exit tunnel", all.imp_df$xgroup2)
  all.imp_df$xgroup2 <- gsub("nt_m.*", "", all.imp_df$xgroup2)  # 0
  all.imp_df$xgroup2 <- gsub("stop_codon", " ", all.imp_df$xgroup2) # 1
  all.imp_df$xgroup2 <- gsub("nt_p.*", "  ", all.imp_df$xgroup2) # 2
  all.imp_df$xgroup2 <- gsub("l_utr3", "   ", all.imp_df$xgroup2) # 3
  all.imp_df$xgroup2 <- gsub("utr3_fraction_.*", "   ", all.imp_df$xgroup2) # 3
  all.imp_df$xgroup2 <- gsub("nis_stop", "    ", all.imp_df$xgroup2) # 4
  all.imp_df$xgroup2 <- gsub("random_.*", "     ", all.imp_df$xgroup2) # 5
  all.imp_df$xgroup2 <- factor(all.imp_df$xgroup2, levels = c("Nascent peptide in the exit tunnel", "", " ", "  ", "   ", "    ", "     "))
  # Tile size for significant level
  all.imp_df$hw <- ifelse(all.imp_df$sig == "ns", yes = 0.6, no = 0.9)
  return(all.imp_df)
}
  # Function to plot heatmap
plot_imp <- function(dat, legend_name = NULL, legend_lim = c(-8, 8), legend_break_length = 9) {
  legend_breaks <- seq(ceiling(legend_lim[1]), floor(legend_lim[2]), length = legend_break_length)
  legend_labs <- paste0(c("< ", rep("", times = legend_break_length-2), "> "), legend_breaks)
  p <- ggplot() +
    geom_tile(data = dat, 
              mapping = aes(x = xtick, y = sample, fill = imp_ave, height = hw, width = hw)) +
    facet_nested(.~xgroup2 + xgroup, scales = "free", space = "free", switch = "y", nest_line = TRUE) +
    scale_y_discrete(limits = rev) +
    scale_fill_gradient2(name = legend_name,
                         low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                         limits = legend_lim, oob = scales::squish,
                         breaks = legend_breaks, labels = legend_labs) +
    xlab("") + ylab("") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "top", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.25, "cm"), 
          legend.box.spacing = unit(0.5, "cm"), legend.title.align = 0,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), axis.text.y = element_text(face = "italic"), axis.title = element_blank(),
          strip.background = element_blank(), strip.text.y = element_blank(), strip.placement = "outside")
  return(p)
}

  # Plot
df_reg <- prep_imp(df_reg)
p_reg <- plot_imp(dat = df_reg, legend_name = "%IncMSE   ", legend_lim = c(-8, 8), legend_break_length = 9)

# Readthrough vs 3'-UTR length ---------------------
df_plot <- read.table("../data/Fig3B_pab1_re_reg_l_utr3_spearman.txt", header = TRUE, sep = "\t")
df_plot$strain <- recode_factor(df_plot$strain, WT = "WT", pbp1d = "pbp1\u0394", pab1d = "pab1\u0394\npbp1\u0394")

set3 <- RColorBrewer::brewer.pal(n = 12, name = "Set3") # Color palette for plotting
p_lutr3 <- ggplot(df_plot, aes(fill = strain, y = cor, x = strain)) + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_point(shape = 21, aes(size = Freq), color = "white") + scale_size_continuous(range = c(2, 10)) +
  geom_text(aes(label = Freq), size = 7/.pt, color = "black") +
  facet_wrap(~stop_codon, nrow = 1) +
  scale_fill_manual(values = set3[c(4, 7, 5)]) +
  xlab("") + ylab("Spearman's correlation coefficient\nReadthrough efficiency vs.\n3'-UTR length (nt)") +
  theme_bw(base_size = 10) + 
  theme(strip.background = element_rect(fill = "white"), strip.placement = "outside", panel.grid = element_blank(),
        legend.position = "none", strip.text.x = element_text(face = "italic"), 
        axis.text.x = element_text(face = "italic"))

# Combine plots
library(patchwork)
p <- p_reg / (p_lutr3 + theme(axis.title.y = element_text(margin = margin(r = -30, unit = "pt")))) + 
  plot_layout(heights = c(1.5, 6), guides = 'keep') + 
  plot_annotation(tag_levels = "A") & 
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
cairo_pdf(filename = "Fig3_readthrough_mRNA_features.pdf", family = "Arial", width = 7.5, height = 6) 
p
dev.off()
