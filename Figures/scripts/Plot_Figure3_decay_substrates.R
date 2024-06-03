# --------------------------------------------------
# Figures 3
# --------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(rstatix)
source("../../Analysis scripts/functions_misc.R")

# RNA-seq results
res_rna <- read.csv("../../Analysis scripts/DE RNA-seq/pab1_RNAseq_DESeq2.csv")
res_rna <- add_sig_fc(df = res_rna, pval_col = "padj", pval_cutoff = 0.01, log2FC_col = "log2FoldChange", log2FC_cutoff = 0,
                      plot_fc_cutoff = 5, plot_pval_cutoff = 10^-50)
res_rna$pair <- factor(res_rna$pair, levels = sort(unique(res_rna$pair))[c(3, 2, 1)])
res_rna$changes <- factor(res_rna$changes, levels = c("Up", "Down", "Unchanged"))

# Decay substrates list
substrates <- list()
substrates$nmd_up <- scan("../../Data/decay_substrates_list/nmd_up_907.txt", character())
substrates$lpd_up <- scan("../../Data/decay_substrates_list/lsm1_pat1_dhh1_up_482.txt", character())
substrates$lp_up <- scan("../../Data/decay_substrates_list/lsm1_pat1_up_382.txt", character())
substrates$d_up <- scan("../../Data/decay_substrates_list/dhh1_up_556.txt", character())
library(ggVennDiagram)
venn <- venn_region(process_data(Venn(substrates)))
venn <- as.list(venn)
venn2 <- venn$item
names(venn2) <- venn$name

res_rna$substrates <- "Non-substrates"
res_rna[res_rna$row %in% venn2$nmd_up, ]$substrates <- "upf1Δ/2Δ/3Δ"
res_rna[res_rna$row %in% venn2$lpd_up, ]$substrates <- "dhh1Δ & pat1Δ/lsm1Δ"
res_rna[res_rna$row %in% venn2$lp_up, ]$substrates <- "pat1Δ/lsm1Δ"
res_rna[res_rna$row %in% venn2$d_up, ]$substrates <- "dhh1Δ"
res_rna[res_rna$row %in% venn2$nmd_up..lpd_up, ]$substrates <- "upf1Δ/2Δ/3Δ & dhh1Δ & pat1Δ/lsm1Δ"
res_rna[res_rna$row %in% venn2$nmd_up..lp_up, ]$substrates <- "upf1Δ/2Δ/3Δ & pat1Δ/lsm1Δ"
res_rna[res_rna$row %in% venn2$nmd_up..d_up, ]$substrates <- "upf1Δ/2Δ/3Δ & dhh1Δ"
res_rna$substrates <- factor(res_rna$substrates, levels = sort(unique(res_rna$substrates))[c(3, 1, 4, 5, 2, 6, 8, 7)])

# Plot substrate status for each column
substatus <- data.frame(substrates = unique(res_rna$substrates))
substatus$Dhh1 <- ifelse(test = grepl(pattern = "dhh1", x = substatus$substrates), yes = TRUE, no = FALSE)
substatus$`Pat1/Lsm1` <- ifelse(test = grepl(pattern = "lsm1", x = substatus$substrates), yes = TRUE, no = FALSE)
substatus$NMD <- ifelse(test = grepl(pattern = "upf", x = substatus$substrates), yes = TRUE, no = FALSE)

substatus2 <- reshape2::melt(substatus, id.vars = "substrates")
substatus2$substrates <- factor(substatus2$substrates, levels = levels(res_rna$substrates))
substatus2$variable <- factor(substatus2$variable, levels = c("NMD", "Pat1/Lsm1", "Dhh1"))

substatus_n <- data.frame(table(substrates = res_rna[which(res_rna$pair == "pab1Δpbp1Δ / pbp1Δ"), ]$substrates))
substatus2 <- left_join(substatus2, substatus_n, by = "substrates")
substatus2$Freq_label <- paste0("n = ", substatus2$Freq)
substatus2$Freq_label <- factor(substatus2$Freq_label, levels = unique(substatus2[order(substatus2$substrates), "Freq_label"]))

ps <- ggplot(substatus2, aes(x = 1, y = variable, fill = value)) +
  geom_tile(width = 0.8, height = 0.8) +
  facet_grid(.~Freq_label, switch = "x") +
  scale_fill_manual(name = "", values = c("white", "forestgreen")) +
  ylab("Substrate") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), #aspect.ratio = 0.4,
        #strip.background = element_blank(), strip.text = element_blank(), 
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing = unit(0, "cm"))

# Plot proportions of substrates Up, Down, Unchanged
nns <- data.frame(table(pair = res_rna$pair, changes = res_rna$changes, substrates = res_rna$substrates))
nns <- nns %>% group_by(pair, substrates) %>% mutate(percent = 100*Freq/sum(Freq))
psr <- ggplot(nns[which(nns$pair == "pab1Δpbp1Δ / pbp1Δ" & nns$Freq > 0), ], 
              aes(x = substrates, y = percent, fill = changes)) + 
  geom_col(position = "stack") + 
  geom_text(aes(label = paste0(Freq, "\n(", round(percent), "%)")), 
            position = position_stack(vjust = 0.5), size = 7/.pt) +
  facet_grid(.~substrates, scales = "free_x", space = "free_x") +
  scale_fill_manual(name = expression(paste("mRNA abundance change in ", italic("pab1Δpbp1Δ / pbp1Δ:"))), 
                    values = c(Up = "orange", Down = "purple", Unchanged = "grey50")) +
  xlab("") + ylab("Percentage") +
  theme_bw(base_size = 10) + 
  theme(strip.background = element_blank(), strip.text = element_blank(), panel.grid = element_blank(), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.spacing = unit(0, "cm")) +
  guides(fill = guide_legend(title.position = "top"))

# TE changes results
res_TE <- read.csv("../../Analysis scripts/TE/pab1_TE_DESeq2.csv")
res_TE <- add_sig_fc(df = res_TE, pval_col = "padj", pval_cutoff = 0.05, log2FC_col = "log2FoldChange", log2FC_cutoff = 0,
                     plot_fc_cutoff = 5, plot_pval_cutoff = 10^-15)
res_TE$changes <- factor(res_TE$changes, levels = c("Up", "Down", "Unchanged"))
PPvP <- left_join(res_TE[which(res_TE$pair == "pab1Δpbp1Δ / pbp1Δ"), ], 
                  res_rna[which(res_rna$pair == "pab1Δpbp1Δ / pbp1Δ"), c("row", "substrates", "changes"), ], 
                  by = "row", suffix = c("_TE", "_rna"))
PPvP2 <- PPvP[which(!is.na(PPvP$substrates)), ]
PPvP2$group <- paste0(PPvP2$substrates, "_", PPvP2$changes_rna)

stat.test <- PPvP2 %>% group_by(group) %>% wilcox_test(log2FoldChange~1, mu = 0, p.adjust.method = "BH") %>% # compare to mu = 0
  adjust_pvalue(method = "BH") %>% add_significance()
# Create placeholders for group where n = 0
missing <- data.frame(group = c("upf1Δ/2Δ/3Δ & dhh1Δ_Down", "upf1Δ/2Δ/3Δ & pat1Δ/lsm1Δ_Down", "upf1Δ/2Δ/3Δ & dhh1Δ & pat1Δ/lsm1Δ_Down"), 
                      .y. = "log2FoldChange", group1 = 1, group2 = "null model", n = 0, statistic = 0, p = 1, p.adj = 1, p.adj.signif = "")
stat.test <- rbind(stat.test, missing)
stat.test$substrates <- sub("_.*", "", stat.test$group)
stat.test$changes_rna <- sub(".*_", "", stat.test$group)
stat.test$substrates <- factor(stat.test$substrates, level = levels(PPvP2$substrates))
stat.test$changes_rna <- factor(stat.test$changes_rna, levels = levels(PPvP2$changes_rna))

nn <- data.frame(table(substrates = PPvP2$substrates, changes_rna = PPvP2$changes_rna))

pte <- ggplot(PPvP2, aes(x = changes_rna, y = log2FoldChange)) + 
  geom_boxplot(aes(fill = changes_rna), outlier.shape = NA) + 
  geom_text(data = nn, aes(x = changes_rna, y = -2.5, label = Freq), size = 7/.pt) +
  geom_text(data = stat.test, aes(x = changes_rna, y = 2.9, label = p.adj.signif), size = 7/.pt) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_grid(.~substrates, scales = "free_x", space = "free_x") +
  scale_fill_manual(name = expression(paste("mRNA abundance change in ", italic("pab1Δpbp1Δ / pbp1Δ:"))), 
                    values = c(Up = "orange", Down = "purple", Unchanged = "grey50")) +
  scale_x_discrete(drop = FALSE) + # For making x-axis order following levels specification for panel with missing level
  xlab(expression(paste("mRNA abundance change in ", italic("pab1Δpbp1Δ / pbp1Δ")))) +
  ylab(expression(atop("log"[2]*" fold change in TE", italic("pab1Δpbp1Δ / pbp1Δ")))) +
  coord_cartesian(ylim = c(-2.5, 3)) +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank(),  
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        panel.spacing = unit(0, "cm"), legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Codon optimality
feature_termination <- read.csv("../../Data/mRNA features/mRNA_features_termination.csv")
PPvP2 <- left_join(PPvP2, feature_termination[, c("transcript", "tAI_cds")], by = c("row" = "transcript"))
PPvP3 <- PPvP2[which(!is.na(PPvP2$tAI_cds)), ] # Remove non-coding RNAs or pre-mRNAs that don't have a CDS to calculate tAI

stat.test2 <- PPvP3 %>% group_by(group) %>% wilcox_test(tAI_cds~1, mu = mean(PPvP3$tAI_cds), p.adjust.method = "BH") %>% # compare to mu = mean of the entire data set
  adjust_pvalue(method = "BH") %>% add_significance()
# Create placeholders for group where n = 0
missing <- data.frame(group = c("upf1Δ/2Δ/3Δ & dhh1Δ_Down", "upf1Δ/2Δ/3Δ & pat1Δ/lsm1Δ_Down", "upf1Δ/2Δ/3Δ & dhh1Δ & pat1Δ/lsm1Δ_Down"), 
                      .y. = "log2FoldChange", group1 = 1, group2 = "null model", n = 0, statistic = 0, p = 1, p.adj = 1, p.adj.signif = "")
stat.test2 <- rbind(stat.test2, missing)
stat.test2$substrates <- sub("_.*", "", stat.test2$group)
stat.test2$changes_rna <- sub(".*_", "", stat.test2$group)
stat.test2$substrates <- factor(stat.test2$substrates, level = levels(PPvP3$substrates))
stat.test2$changes_rna <- factor(stat.test2$changes_rna, levels = levels(PPvP3$changes_rna))

no <- data.frame(table(substrates = PPvP3$substrates, changes_rna = PPvP3$changes_rna))

pco <- ggplot(PPvP3, aes(x = changes_rna, y = tAI_cds)) + 
  geom_boxplot(aes(fill = changes_rna), outlier.shape = NA) + 
  geom_text(data = no, aes(x = changes_rna, y = 0.2, label = Freq), size = 7/.pt) +
  geom_text(data = stat.test2, aes(x = changes_rna, y = 0.6, label = p.adj.signif), size = 7/.pt) +
  geom_hline(yintercept = mean(PPvP3$tAI_cds), color = "red", linetype = "dashed") +
  facet_grid(.~substrates, scales = "free_x", space = "free_x") +
  scale_fill_manual(name = expression(paste("mRNA abundance change in ", italic("pab1Δpbp1Δ / pbp1Δ:"))), 
                    values = c(Up = "orange", Down = "purple", Unchanged = "grey50")) +
  scale_x_discrete(drop = FALSE) + # For making x-axis order following levels specification for panel with missing level
  xlab(expression(paste("mRNA abundance change in ", italic("pab1Δpbp1Δ / pbp1Δ")))) +
  ylab("Codon optimality") +
  coord_cartesian(ylim = c(0.2, 0.6)) +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_blank(),  
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing = unit(0, "cm"), legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"))


# Combine plots
library(patchwork)
pA <- (ps + psr) + 
  plot_layout(nrow = 2, heights = c(0.2, 0.8)) & 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "none")
p <- (pA / pte / pco) +
  plot_layout(heights = c(0.4, 0.3, 0.3)) +
  plot_annotation()

# Export plots
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Fig3_decay_substrates.pdf", family = "Arial", width = 7.5, height = 7.5) 
p
dev.off()
