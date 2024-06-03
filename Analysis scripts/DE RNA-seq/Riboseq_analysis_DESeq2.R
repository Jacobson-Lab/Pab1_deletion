# Analysis of mRibo abundance changes from Ribo-Seq experiment

library(data.table)
library(dplyr)
library(DESeq2)

wrapper_deseq2 <- function(input_file_path) {
  # countData
  tab <- read.csv(input_file_path)
  rownames(tab) <- tab$transcript_id
  tab <- tab[, colnames(tab)[grepl(pattern = "^RPF.*", x = colnames(tab))]]
  tab <- round(x = tab, digits = 0)
  tab <- tab[rowSums(tab) > 0, ] # Remove rows with no reads in any samples
  
  # colData
  cdt <- data.frame(fullname = colnames(tab))
  cdt$strain <- sub("RPF_", "", cdt$fullname)
  cdt$strain <- sub("_rep[0-9]", "", cdt$strain)
  
  #DESeq
  dds <- DESeqDataSetFromMatrix(countData = tab, colData = cdt, design = ~strain)
  dds$strain <- relevel(dds$strain, ref = "WT")
  DGE <- DESeq(dds)
  
  # results
  res_t <- list()
  res_t[["pbp1\u0394 / WT"]] <- results(DGE, contrast = c("strain", "pbp1d", "WT"), pAdjustMethod = "fdr", tidy = TRUE, independentFiltering = TRUE, alpha = 0.01)
  res_t[["pab1\u0394pbp1\u0394 / WT"]] <- results(DGE, contrast = c("strain", "pab1d_pbp1d", "WT"), pAdjustMethod = "fdr", tidy = TRUE, independentFiltering = TRUE, alpha = 0.01)
  res_t[["pab1\u0394pbp1\u0394 / pbp1\u0394"]] <- results(DGE, contrast = c("strain","pab1d_pbp1d", "pbp1d"), pAdjustMethod = "fdr", tidy = TRUE, independentFiltering = TRUE, alpha = 0.01)
  
  return(res_t)
}

# isoforms results
res_t <- wrapper_deseq2(input_file_path = "../../Data/RNAseq-Riboseq/pab1_strains_RSEM_expected_count_isoforms.results.csv")
openxlsx::write.xlsx(list('pbp1d vs WT' = res_t$`pbp1Δ / WT`, 'pab1d_pbp1d vs WT' = res_t$`pab1Δpbp1Δ / WT`, "pab1d_pbp1d vs pbp1d" = res_t$`pab1Δpbp1Δ / pbp1Δ`), 
                     file = "pab1_Riboseq_DESeq2.xlsx")
res_t2 <- bind_rows(res_t, .id = "pair")
write.csv(res_t2, file = "pab1_Riboseq_DESeq2.csv", row.names = FALSE, quote = FALSE)

# isoforms results - combine reads from protein cluster to compare with mass spec results
res_t <- wrapper_deseq2(input_file_path = "../../Data/RNAseq-Riboseq/pab1_strains_RSEM_expected_count_isoforms.results_combine_for_protein_cluster.csv")
openxlsx::write.xlsx(list('pbp1d vs WT' = res_t$`pbp1Δ / WT`, 'pab1d_pbp1d vs WT' = res_t$`pab1Δpbp1Δ / WT`, "pab1d_pbp1d vs pbp1d" = res_t$`pab1Δpbp1Δ / pbp1Δ`), 
                     file = "pab1_Riboseq_DESeq2_combine_for_protein_cluster.xlsx")
res_t2 <- bind_rows(res_t, .id = "pair")
write.csv(res_t2, file = "pab1_Riboseq_DESeq2_combine_for_protein_cluster.csv", row.names = FALSE, quote = FALSE)
