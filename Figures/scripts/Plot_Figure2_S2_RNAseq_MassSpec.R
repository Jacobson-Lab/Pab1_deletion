# --------------------------------------------------
# Figures 2 and S2
# --------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
source("../../Analysis scripts/functions_misc.R")

# RNA-seq results
res_rna <- read.csv("../../Analysis scripts/DE RNA-seq/pab1_RNAseq_DESeq2.csv")
res_rna <- add_sig_fc(df = res_rna, pval_col = "padj", pval_cutoff = 0.01, log2FC_col = "log2FoldChange", log2FC_cutoff = 0,
                      plot_fc_cutoff = 5, plot_pval_cutoff = 10^-50)
res_rna$pair <- factor(res_rna$pair, levels = sort(unique(res_rna$pair))[c(3, 2, 1)])
res_rna$changes <- factor(res_rna$changes, levels = c("Up", "Down", "Unchanged"))
nn <- data.frame(table(pair = res_rna$pair, changes = res_rna$changes))
nn$ypos <- rep(c(30, 25, 20), each = 3)

p_rna_main <- plot_volcano(res = res_rna[which(res_rna$pair == "pab1Δpbp1Δ / pbp1Δ"), ], 
                           nn = nn[which(nn$pair == "pab1Δpbp1Δ / pbp1Δ"), ], xlabel = "in mRNA abundance", 
                           point_label = c("PAB1", "PBP1"))
p_rna_supp <- plot_volcano(res = res_rna[which(res_rna$pair != "pab1Δpbp1Δ / pbp1Δ"), ], 
                           nn = nn[which(nn$pair != "pab1Δpbp1Δ / pbp1Δ"), ], xlabel = "in mRNA abundance")

# Mass spec results
res_prot <- read.csv("../../Analysis scripts/DE Mass spec/pab1_mass_spec_limma.csv")
res_prot <- res_prot[!is.na(res_prot$logFC), ]
res_prot <- add_sig_fc(df = res_prot, pval_col = "adj.P.Val", pval_cutoff = 0.015, log2FC_col = "logFC", log2FC_cutoff = 0,
                       plot_fc_cutoff = 2, plot_pval_cutoff = 10^-50)
res_prot$pair <- factor(res_prot$pair, levels = sort(unique(res_prot$pair))[c(3, 2, 1)])
res_prot$changes <- factor(res_prot$changes, levels = c("Up", "Down", "Unchanged"))
nnp <- data.frame(table(pair = res_prot$pair, changes = res_prot$changes))
nnp$ypos <- rep(c(1.5, 1, 0.5), each = 3)

p_prot_main <- plot_volcano(res = res_prot[which(res_prot$pair == "pab1Δpbp1Δ / pbp1Δ"), ], 
                            nn = nnp[which(nnp$pair == "pab1Δpbp1Δ / pbp1Δ"), ], xlabel = "in protein abundance",
                            point_label = c("PAB1", "PBP1"))
p_prot_supp <- plot_volcano(res = res_prot[which(res_prot$pair != "pab1Δpbp1Δ / pbp1Δ"), ], 
                            nn = nnp[which(nnp$pair != "pab1Δpbp1Δ / pbp1Δ"), ], xlabel = "in protein abundance",
                            point_label = c("PAB1", "PBP1"))

# RNA-seq vs. Mass spec
ID <- read.table("../../Data/yeast_transcript_gene_IDs.txt", header = TRUE)
colnames(ID) <- c("row", "Protein")

res_rna_cluster <- read.csv("../../Analysis scripts/DE RNA-seq/pab1_RNAseq_DESeq2_combine_for_protein_cluster.csv")
res_rna_cluster <- add_sig_fc(df = res_rna_cluster, pval_col = "padj", pval_cutoff = 0.01, log2FC_col = "log2FoldChange", log2FC_cutoff = 0,
                              plot_fc_cutoff = 5, plot_pval_cutoff = 10^-50)
res_rna_cluster <- left_join(res_rna_cluster, ID, by = "row")
res_rna_cluster$mrna <- ifelse(test = grepl(pattern = "pre-mRNA", x = res_rna_cluster$row), yes = "pre-mRNA", no = "mRNA")
res_rna_cluster <- res_rna_cluster[which(res_rna_cluster$mrna == "mRNA"), ] # Protein data don't have pre-mRNA equivalent
res_rna_cluster$pair <- factor(res_rna_cluster$pair, levels = sort(unique(res_rna_cluster$pair))[c(3, 2, 1)])

combined <- full_join(res_rna_cluster, res_prot, by = c("pair", "Protein"), suffix = c("_mRNA", "_protein"))
combined$sig_both <- "Neither"
combined[which(combined$sig_mRNA == "ns" & combined$sig_protein != "ns"), ]$sig_both <- "Protein only"
combined[which(combined$sig_mRNA != "ns" & combined$sig_protein == "ns"), ]$sig_both <- "mRNA only"
combined[which(combined$sig_mRNA != "ns" & combined$sig_protein != "ns"), ]$sig_both <- "Both"
combined$sig_both <- factor(combined$sig_both, levels = c("Neither", "Protein only", "mRNA only", "Both"))

combined_filter <- combined[which(!is.na(combined$logFC) & !is.na(combined$log2FoldChange)), ]
combined_filter$repel <- ""
point_label = c("PAB1", "PBP1")
combined_filter[which(combined_filter$Protein %in% point_label), ]$repel <- combined_filter[which(combined_filter$Protein %in% point_label), ]$Protein
nnc <- data.frame(table(pair = combined_filter$pair, sig_both = combined_filter$sig_both))
nnc$ypos <- rep(c(1.3, 1.0, 0.7, 0.4), each = 3)

plot_scatter_main <- plot_scatter(combined_filter = combined_filter[which(combined_filter$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                                  nnc = nnc[which(nnc$pair == "pab1Δpbp1Δ / pbp1Δ"), ], facet_group = "pair")
plot_scatter_supp <- plot_scatter(combined_filter = combined_filter[which(combined_filter$pair != "pab1Δpbp1Δ / pbp1Δ"), ],
                                  nnc = nnc[which(nnc$pair != "pab1Δpbp1Δ / pbp1Δ"), ], facet_group = "pair")

  # Facet plot by groups of translation-related proteins
prot_group <- list()
prot_group[["Initiation factors"]] <- scan("../../Data/protein_group_list/translation_initiation_factor_genes.txt", character())
prot_group[["Elongation factors"]] <- scan("../../Data/protein_group_list/translation_elongation_factor_genes.txt", character())
prot_group[["Ribosomal proteins"]] <- scan("../../Data/protein_group_list/ribosomal_protein_genes.txt", character())
prot_group[["Release factors and recycling factors"]] <- scan("../../Data/protein_group_list/translation_termination_recycling_factor_genes.txt", character())
prot_group[["Decapping/5'-3' Decay"]] <- scan("../../Data/protein_group_list/decay5.txt", character())
prot_group[["Deadenylation/3'-5' Decay"]] <- scan("../../Data/protein_group_list/decay3.txt", character())

combined_filter$prot_group <- ""
for (x in names(prot_group)) {
  combined_filter[which(combined_filter$Protein %in% prot_group[[x]]), ]$prot_group <- x
}
prot_group_lev <- c("", names(prot_group))
combined_filter$prot_group <- factor(combined_filter$prot_group, levels = prot_group_lev)
combined_filter$repel <- ""

point_label <- data.frame(gene_name = unlist(prot_group[c(1:2, 4:6)]))
protein_name <- readxl::read_xlsx("../../Data/protein_name.xlsx")
point_label <- left_join(point_label, protein_name, by = "gene_name")
point_label[is.na(point_label$protein_name), ]$protein_name <- point_label[is.na(point_label$protein_name), ]$gene_name
combined_filter <- left_join(combined_filter, point_label, by = c("Protein" = "gene_name"))
combined_filter[which(combined_filter$Protein %in% point_label$gene_name), ]$repel <- combined_filter[which(combined_filter$Protein %in% point_label$gene_name), ]$protein_name

plot_prot_group_translation <- plot_scatter(combined_filter = combined_filter[which(combined_filter$prot_group %in% prot_group_lev[2:5] & combined_filter$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                                            facet_group = "prot_group")
plot_prot_group_translation <- plot_prot_group_translation + 
  scale_x_continuous(limits = c(-1.55, 1.55)) + scale_y_continuous(limits = c(-1.55, 1.55)) + 
  theme(legend.position = "none") +
  xlab(expression(atop("log"[2]*" fold change in mRNA abundance", italic("pab1Δpbp1Δ / pbp1Δ")))) +
  ylab(expression(atop("log"[2]*" fold change in protein abundance", italic("pab1Δpbp1Δ / pbp1Δ"))))

  # (For Fig 6B)
plot_prot_group_decay <- plot_scatter(combined_filter = combined_filter[which(combined_filter$prot_group %in% prot_group_lev[6:7] & combined_filter$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                                      facet_group = "prot_group")
plot_prot_group_decay <- plot_prot_group_decay + 
  xlab(expression(atop("log"[2]*" fold change in mRNA abundance", italic("pab1Δpbp1Δ / pbp1Δ")))) +
  ylab(expression(atop("log"[2]*" fold change in protein abundance", italic("pab1Δpbp1Δ / pbp1Δ"))))

# Combine plots
library(patchwork)
p_main <- (p_rna_main + p_prot_main + plot_scatter_main) + guide_area() + 
  plot_layout(guides = 'collect', ncol = 2, nrow = 2) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")
p_main2 <- (p_main / plot_prot_group_translation) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft")
p_supp <- (p_rna_supp / p_prot_supp / plot_scatter_supp) + plot_layout(guides = 'keep') +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")

# Export plots
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Fig2_RNAseq_protein_DE_analyses.pdf", family = "Arial", width = 6, height = 10) 
p_main2
dev.off()
cairo_pdf(filename = "FigS2_RNAseq_protein_DE_analyses.pdf", family = "Arial", width = 7, height = 7.5) 
p_supp
dev.off()
