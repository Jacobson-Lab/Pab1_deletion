# --------------------------------------------------
# Figure 4
# --------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(rstatix)
library(patchwork)
set3 <- RColorBrewer::brewer.pal(n = 12, name = "Set3") # Color palette for plotting

# A: Metaprofile
df <- read.table("../data/metaprofile.txt", header = TRUE, sep = "\t")
df$strain <- factor(df$strain, levels = sort(unique(df$strain))[c(3, 2, 1)])
  # Define vertical line for every in-frame nt
vline <- data.frame(region = c(rep("Distance from start (nt)", times = 22), rep("Distance from stop (nt)", times = 40)), 
                    xin = c(seq(-24, 39, 3), seq(-38, 79, 3)))
  # Define vertical line for the fist nucleotides of start and stop codons
rline <- data.frame(region = c("Distance from start (nt)", "Distance from stop (nt)"), xin = c(0, 1))

p_main <- ggplot() +
  geom_vline(data = vline, aes(xintercept = xin), linetype = "dotted", color = "grey50", linewidth = 0.2) +
  geom_vline(data = rline, aes(xintercept = xin), color = "grey50", linewidth = 0.2) +
  geom_line(data = df, aes(x = distance, y = fraction, color = strain), linewidth = 0.5) +
  facet_grid(.~region, scales = "free_x", space = "free_x", switch = "x") +
  scale_color_manual(name = "", values = set3[c(4, 7, 5)]) +
  scale_x_continuous(breaks = seq(-20, 80, 20), labels = seq(-20, 80, 20), expand = expansion(mult = 0.02)) +
  xlab("") + ylab("Normalized read count") + 
  coord_cartesian(ylim = c(0, 0.12)) +
  theme_bw(base_size = 10) + 
  theme(strip.background = element_rect(fill = "white", color = NA), panel.grid = element_blank(), strip.placement = "outside",
        legend.position = "top", legend.text = element_text(face = "italic")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1)))
p_inset <- ggplot() +
  geom_vline(xintercept = c(seq(-44, 85, 3)), linetype = "dotted", color = "grey", linewidth = 0.2) +
  geom_line(data = df[df$region == "Distance from stop (nt)", ], aes(x = distance, y = fraction, color = strain), linewidth = 0.5) +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.01))) +
  scale_color_manual(name = "", values = set3[c(4, 7, 5)]) +
  coord_cartesian(ylim = c(0, 0.001), xlim = c(1, 80)) +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.ticks.length = unit(0, "pt"),
        legend.position = "none", strip.text = element_blank(), panel.spacing.y = unit(0, "cm"), 
        #plot.margin = grid::unit(c(0, 0, 0, 0), "mm"), 
        plot.background = element_rect(color = "black"))
pA <- p_main + inset_element(p_inset, left = 0.6, bottom = 0.1, right = 0.98, top = 0.95, ignore_tag = TRUE)

# B-C: Percentage of reads in mRNA regions
fr <- data.table(read.csv("../data/psite_region_frame.csv"))
fr[, total := sum(count), by = list(sample)]
fr[, total_region := sum(count), by = list(sample, region)]
fr$percentage <- 100*(fr$count/fr$total)
fr$percentage_region <- 100*(fr$total_region/fr$total)
fr$percentage_frame <- 100*(fr$count/fr$total_region)
fr$strain <- sub("_[0-9].*", "", fr$sample)
fr$rep <- sub(".*[d|T]_", "rep_", fr$sample)
fr$strain <- recode_factor(fr$strain, WT = "WT", pbp1d = "pbp1\u0394", pab1d = "pab1\u0394pbp1\u0394")
fr$region <- recode_factor(fr$region, `5utr` = "5'-UTR", start = "Start", cds = "CDS", stop = "Stop", extension = "Extension", `distal_3utr` = "Distal 3'-UTR")
fr$frame <- paste0("Frame ", fr$frame)

  # Percentage of reads in each mRNA region
fr_filter <- unique(fr[, c("sample", "region", "percentage_region", "strain", "rep")])
stat.test <- fr_filter %>%
  group_by(region) %>%
  t_test(percentage_region~strain, p.adjust.method = "BH", paired = FALSE,
         comparisons = list(c("WT", "pbp1\u0394"), c("pbp1\u0394", "pab1\u0394pbp1\u0394"), c("WT", "pab1\u0394pbp1\u0394"))) %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "strain", scales = "free_y")

pB <- ggplot(fr_filter, aes(x = strain, y = percentage_region)) +
  stat_summary(fun = mean, geom = "bar", aes(fill = strain), color = NA) + # Add mean bar
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), size = 6/.pt, vjust = -0.3, color = "black") + # Add mean number above bar
  geom_point(color = "grey25", position = position_dodge2(width = 0.75), size = 0.75, alpha = 0.5, show.legend = FALSE) +
  stat_pvalue_manual(data = stat.test, label = "p.adj.signif", tip.length = 0.01, size = 7/.pt, vjust = 0.5, hide.ns = TRUE) +
  facet_wrap(~region, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "", values = set3[c(4, 7, 5)]) +
  xlab("") + ylab("% Footprint") +
  theme_bw(base_size = 10) + 
  theme(strip.background = element_rect(fill = "white"), strip.placement = "outside",  panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.position = "none", legend.text = element_text(face = "italic"))

# Percentage of reads in each reading frame in each mRNA region
stat.test2 <- fr[fr$frame == "Frame 0", ] %>%
  group_by(region) %>%
  t_test(percentage_frame~strain, p.adjust.method = "BH", paired = FALSE,
         comparisons = list(c("WT", "pbp1\u0394"), c("pbp1\u0394", "pab1\u0394pbp1\u0394"), c("WT", "pab1\u0394pbp1\u0394"))) %>%
  add_significance()
stat.test2 <- stat.test2 %>% add_xy_position(x = "strain")
stat.test2$y.position <- stat.test2$y.position + 5

pC <- ggplot(fr[fr$frame == "Frame 0", ], aes(x = strain, y = percentage_frame)) +
  stat_summary(fun = mean, geom = "bar", aes(fill = strain), color = NA) + # Add mean bar
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), size = 6/.pt, vjust = -0.5, color = "black") + # Add mean number above bar
  geom_point(color = "grey25", position = position_dodge2(width = 0.75), size = 0.75, alpha = 0.5, show.legend = FALSE) +
  stat_pvalue_manual(data = stat.test2, label = "p.adj.signif", tip.length = 0.01, size = 7/.pt, vjust = 0.5, hide.ns = TRUE) +
  geom_hline(yintercept = (1/3)*100, color = "gray50", linetype = "dashed", size = 0.2) +
  facet_nested_wrap(~region, nrow = 1, axes = "y") +
  scale_fill_manual(name = "", values = set3[c(4, 7, 5)]) +
  xlab("") + ylab("% Frame 0 footprint\nin mRNA region") +
  theme_bw(base_size = 10) + 
  theme(strip.background = element_rect(fill = "white"), strip.placement = "outside", panel.grid = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.position = "none", legend.text = element_text(face = "italic"))

# Combine plots
library(patchwork)
p <- (pA + pB + pC) + plot_layout(guides = 'keep', ncol = 1, nrow = 3) +
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
cairo_pdf(filename = "Fig4_metaprofile.pdf", family = "Arial", width = 7.5, height = 7.5) 
p
dev.off()
