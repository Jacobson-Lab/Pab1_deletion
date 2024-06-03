# --------------------------------------------------
# Figures 1, S2, and 3
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

# RNA-seq or Ribo-seq vs. Mass spec
combine_data <- function(input_sequencing_file, res_prot = res_prot) {
  ID <- read.table("../../Data/yeast_transcript_gene_IDs.txt", header = TRUE)
  colnames(ID) <- c("row", "Protein")
  
  res_cluster <- read.csv(input_sequencing_file)
  res_cluster <- add_sig_fc(df = res_cluster, pval_col = "padj", pval_cutoff = 0.01, log2FC_col = "log2FoldChange", log2FC_cutoff = 0,
                            plot_fc_cutoff = 5, plot_pval_cutoff = 10^-50)
  res_cluster <- left_join(res_cluster, ID, by = "row")
  res_cluster$mrna <- ifelse(test = grepl(pattern = "pre-mRNA", x = res_cluster$row), yes = "pre-mRNA", no = "mRNA")
  res_cluster <- res_cluster[which(res_cluster$mrna == "mRNA"), ] # Protein data don't have pre-mRNA equivalent
  res_cluster$pair <- factor(res_cluster$pair, levels = sort(unique(res_cluster$pair))[c(3, 2, 1)])
  
  combined <- full_join(res_cluster, res_prot, by = c("pair", "Protein"), suffix = c("_mRNA", "_protein"))
  combined$sig_both <- "Neither"
  combined[which(combined$sig_mRNA == "ns" & combined$sig_protein != "ns"), ]$sig_both <- "Protein only"
  combined[which(combined$sig_mRNA != "ns" & combined$sig_protein == "ns"), ]$sig_both <- "mRNA only"
  combined[which(combined$sig_mRNA != "ns" & combined$sig_protein != "ns"), ]$sig_both <- "Both"
  combined$sig_both <- factor(combined$sig_both, levels = c("Neither", "Protein only", "mRNA only", "Both"))
  
  combined_filter <- combined[which(!is.na(combined$logFC) & !is.na(combined$log2FoldChange)), ]
  combined_filter$repel <- ""
  point_label = c("PAB1", "PBP1")
  combined_filter[which(combined_filter$Protein %in% point_label), ]$repel <- combined_filter[which(combined_filter$Protein %in% point_label), ]$Protein
  return(combined_filter)
}

  # RNA-seq vs. Mass spec
combined_filter <- combine_data(input_sequencing_file = "../../Analysis scripts/DE RNA-seq/pab1_RNAseq_DESeq2_combine_for_protein_cluster.csv", res_prot = res_prot)
nnc <- data.frame(table(pair = combined_filter$pair, sig_both = combined_filter$sig_both))
nnc$ypos <- rep(c(1.3, 1.0, 0.7, 0.4), each = 3)

plot_scatter_main <- plot_scatter(combined_filter = combined_filter[which(combined_filter$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                                  nnc = nnc[which(nnc$pair == "pab1Δpbp1Δ / pbp1Δ"), ], facet_group = "pair", xlabel = "in mRNA abundance")
plot_scatter_supp <- plot_scatter(combined_filter = combined_filter[which(combined_filter$pair != "pab1Δpbp1Δ / pbp1Δ"), ],
                                  nnc = nnc[which(nnc$pair != "pab1Δpbp1Δ / pbp1Δ"), ], facet_group = "pair", xlabel = "in mRNA abundance")

  # Ribo-seq vs. Mass spec
combined_filter2 <- combine_data(input_sequencing_file = "../../Analysis scripts/DE Ribo-seq/pab1_Riboseq_DESeq2_combine_for_protein_cluster.csv", res_prot = res_prot)
nnc2 <- data.frame(table(pair = combined_filter2$pair, sig_both = combined_filter2$sig_both))
nnc2$ypos <- rep(c(1.3, 1.0, 0.7, 0.4), each = 3)

plot_scatter_main2 <- plot_scatter(combined_filter = combined_filter2[which(combined_filter2$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                                  nnc = nnc2[which(nnc2$pair == "pab1Δpbp1Δ / pbp1Δ"), ], facet_group = "pair", xlabel = "Ribo-seq reads")
plot_scatter_supp2 <- plot_scatter(combined_filter = combined_filter2[which(combined_filter2$pair != "pab1Δpbp1Δ / pbp1Δ"), ],
                                  nnc = nnc2[which(nnc2$pair != "pab1Δpbp1Δ / pbp1Δ"), ], facet_group = "pair", xlabel = "Ribo-seq reads")

# Facet plot by groups of related proteins
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

  # Translation-related proteins
plot_prot_group_translation <- plot_scatter(combined_filter = combined_filter[which(combined_filter$prot_group %in% prot_group_lev[2:5] & combined_filter$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                                            facet_group = "prot_group")
plot_prot_group_translation <- plot_prot_group_translation + 
  scale_x_continuous(limits = c(-1.55, 1.55)) + scale_y_continuous(limits = c(-1.55, 1.55)) + 
  theme(legend.position = "none") +
  xlab(expression(atop("log"[2]*" fold change in mRNA abundance", italic("pab1Δpbp1Δ / pbp1Δ")))) +
  ylab(expression(atop("log"[2]*" fold change in protein abundance", italic("pab1Δpbp1Δ / pbp1Δ"))))

  # mRNA Decay-related proteins
plot_prot_group_decay <- plot_scatter(combined_filter = combined_filter[which(combined_filter$prot_group %in% prot_group_lev[6:7] & combined_filter$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                                      facet_group = "prot_group")
plot_prot_group_decay <- plot_prot_group_decay + 
  xlab(expression(atop("log"[2]*" fold change in mRNA abundance", italic("pab1Δpbp1Δ / pbp1Δ")))) +
  ylab(expression(atop("log"[2]*" fold change in protein abundance", italic("pab1Δpbp1Δ / pbp1Δ"))))

# Facet by TE
  # First, run code in translation_efficiency_analysis_DESeq2.Rmd to get res_TE
combined_filter_rna <- left_join(combined_filter, res_TE, by = c("pair", "row"), suffix = c("", "_TE"))
combined_filter_rna$changes2 <- ifelse(test = combined_filter_rna$changes == "Unchanged", yes = "TE unchanged", no = "TE changed")
p_rna_te <- plot_scatter(combined_filter = combined_filter_rna[which(combined_filter_rna$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                         facet_group = "changes2", xlabel = "in mRNA abundance")

combined_filter_ribo <- left_join(combined_filter2, res_TE, by = c("pair", "row"), suffix = c("", "_TE"))
combined_filter_ribo$changes2 <- ifelse(test = combined_filter_ribo$changes == "Unchanged", yes = "TE unchanged", no = "TE changed")
p_ribo_te <- plot_scatter(combined_filter = combined_filter_ribo[which(combined_filter_ribo$pair == "pab1Δpbp1Δ / pbp1Δ"), ],
                         facet_group = "changes2", xlabel = "Ribo-seq reads")

# Combine and export plots
library(patchwork)
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)

p_main1 <- (p_prot_main + p_rna_main) + 
  plot_layout(guides = 'collect', ncol = 2, nrow = 1) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")
p_main2 <- (plot_scatter_main + plot_scatter_main2) + 
  plot_layout(guides = 'collect', ncol = 2, nrow = 1) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")
p_main <- (p_main1 / p_main2) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")
cairo_pdf(filename = "Fig1_RNAseq_Riboseq_protein_DE_analyses.pdf", family = "Arial", width = 7.5, height = 5) 
p_main
dev.off()

p_supp <- (p_prot_supp / p_rna_supp / plot_scatter_supp / plot_scatter_supp2) + plot_layout(guides = 'keep') +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")
cairo_pdf(filename = "FigS2_RNAseq_Riboseq_protein_DE_analyses.pdf", family = "Arial", width = 7.5, height = 10) 
p_supp
dev.off()

p_groups <- (plot_prot_group_translation / plot_prot_group_decay) + plot_layout(guides = 'keep', heights = c(2, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")
cairo_pdf(filename = "Fig2_RNAseq_protein_DE_groups.pdf", family = "Arial", width = 7.5, height = 8) 
p_groups
dev.off()

p_te <- (p_rna_te + p_ribo_te) +
  plot_layout(guides = 'collect', ncol = 1, nrow = 2) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "right")
cairo_pdf(filename = "FigS5_RNAseq_protein_DE_TE.pdf", family = "Arial", width = 7.5, height = 5) 
p_te
dev.off()
