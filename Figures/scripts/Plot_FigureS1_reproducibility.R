# --------------------------------------------------
# Figure S1
# --------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggh4x)

# RNA-seq, Ribo-seq libraries
tab <- read.csv("../../Data/RNAseq-Riboseq/pab1_strains_RSEM_FPKM_genes.results.csv", row.names = 1)
tab <- tab[, colnames(tab)[grepl(pattern = "^R.*", x = colnames(tab))]]
tab <- log10(tab)
tab <- do.call(data.frame, lapply(tab, function(x) replace(x, is.infinite(x), NA)))

cor_rna <- cor(tab, method = "pearson", use = "pairwise.complete.obs")
cor_rna <- reshape2::melt(cor_rna)
cor_rna$strain1 <- sub("_rep.*", "", cor_rna$Var1)
cor_rna$strain1 <- sub(".*[AF]_", "", cor_rna$strain1)
cor_rna$type1 <- sub("_.*", "", cor_rna$Var1)
cor_rna$rep1 <- sub(".*rep", "", cor_rna$Var1)
cor_rna$strain2 <- sub("_rep.*", "", cor_rna$Var2)
cor_rna$strain2 <- sub(".*[AF]_", "", cor_rna$strain2)
cor_rna$type2 <- sub("_.*", "", cor_rna$Var2)
cor_rna$rep2 <- sub(".*rep", "", cor_rna$Var2)
cor_rna$strain1 <- recode_factor(cor_rna$strain1, WT = "WT", pbp1d = "pbp1\u0394", pab1d_pbp1d = "pab1\u0394pbp1\u0394")
cor_rna$strain2 <- recode_factor(cor_rna$strain2, WT = "WT", pbp1d = "pbp1\u0394", pab1d_pbp1d = "pab1\u0394pbp1\u0394")

p_rna <- ggplot(cor_rna, aes(x = rep1, y = rep2)) +
  geom_tile(aes(fill = value), height = 0.9, width = 0.9, color = NA) +
  geom_text(aes(label = round(value, digits = 2)), size = 7/.pt) +
  scale_fill_distiller(name = "Pearson's r  ", palette = "Spectral", direction = -1, 
                       limits = c(0.7, 1)) +
  scale_y_discrete(limits = rev) +
  facet_nested(type2+strain2~type1+strain1, scales = "free", space = "free") +
  xlab("") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), 
        panel.border = element_rect(linewidth = 0.2),
        strip.background = element_rect(fill = "white"), strip.text = element_text(face = "italic"), 
        legend.position = "top", legend.key.height = unit(0.25, "cm"), legend.key.width = unit(1.5, "cm"), aspect.ratio = 1)

# Mass spec
df <- readxl::read_xlsx("../../Data/Mass spec/Mass_spec_TMT10plex_offlinefractionated_log2_normalized_intensity_median_normalization.xlsx")
TMTcols <- colnames(df)[grepl(pattern = "^TMT.*", x = colnames(df))]
df <- df[, TMTcols]
colnames(df) <- c("WT_rep1", "WT_rep3", "WT_rep2", "pab1d_pbp1d_rep1", "pab1d_pbp1d_rep3", "pab1d_pbp1d_rep2", "pbp1d_rep2", "pbp1d_rep1", "pbp1d_rep3")
df <- df[-which(df$WT_rep1 == "No data"), ]   # Remove rows with no data
df[df == "Value missing"] <- NA               # Replace "Value missing" with NA
df <- sapply(df, as.numeric)                  # Convert data to numeric
df <- df[rowSums(is.na(df)) != ncol(df), ]    # Remove rows with all NA

cor_protein <- cor(df, method = "pearson", use = "pairwise.complete.obs")
cor_protein <- reshape2::melt(cor_protein)
cor_protein$strain1 <- sub("_rep.*", "", cor_protein$Var1)
cor_protein$rep1 <- sub(".*rep", "", cor_protein$Var1)
cor_protein$type1 <- "Protein"
cor_protein$strain2 <- sub("_rep.*", "", cor_protein$Var2)
cor_protein$rep2 <- sub(".*rep", "", cor_protein$Var2)
cor_protein$type2 <- "Protein"
cor_protein$strain1 <- recode_factor(cor_protein$strain1, WT = "WT", pbp1d = "pbp1\u0394", pab1d_pbp1d = "pab1\u0394pbp1\u0394")
cor_protein$strain2 <- recode_factor(cor_protein$strain2, WT = "WT", pbp1d = "pbp1\u0394", pab1d_pbp1d = "pab1\u0394pbp1\u0394")

p_prot <- ggplot(cor_protein, aes(x = rep1, y = rep2)) +
  geom_tile(aes(fill = value), height = 0.9, width = 0.9, color = NA) +
  geom_text(aes(label = round(value, digits = 2)), size = 7/.pt) +
  scale_fill_distiller(name = "Pearson's r  ", palette = "Spectral", direction = -1, 
                       limits = c(0.9, 1)) +
  scale_y_discrete(limits = rev) +
  facet_nested(type2+strain2~type1+strain1, scales = "free", space = "free") +
  xlab("") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), 
        panel.border = element_rect(linewidth = 0.2),
        strip.background = element_rect(fill = "white"), strip.text = element_text(face = "italic"), 
        legend.position = "top", legend.key.height = unit(0.25, "cm"), legend.key.width = unit(1.5, "cm"),  aspect.ratio = 1)

library(patchwork)
pp <- p_rna + p_prot + plot_layout(nrow = 2, heights = c(2, 1)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft")

library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "FigS1_reproducibility.pdf", family = "Arial", width = 6, height = 10) 
pp
dev.off()
