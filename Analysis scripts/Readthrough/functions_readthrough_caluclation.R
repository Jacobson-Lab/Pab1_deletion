#--------------------------------------------------
# Functions to calculate readthrough efficiency
#--------------------------------------------------
# Count reads in each mRNA's CDS and extension region
  # dt: data table from bam file, read into R by riboWaltz
  # annotation: file contaning mRNA region length and position of first in-frame stop codon in the 3'-UTR
  # cds_m5: the number of nucleotides from the start codon (5' end of CDS, start codon inclusive) to EXCLUDE from the CDS count
  # cds_m3: the number of nucleotides from the stop codon (3' end of CDS, stop codon exclusive) to EXCLUDE from the CDS count
  # utr3_p: the number of nucleotides to INCLUDE in the 3'UTR read count, in the case that a fixed distance to calculate readthrough is desired. 
  #   Default (Inf) means no constraint, and reads will be counted until "nis_pos", which the next in-frame stop codon or end of 3'UTR. 
  #   If a number is specified, the count will be until the minimum of the specified number, nis_pos, or l_utr3
  # frame_: reading frame to include in the count. Options are "all", 0, 1, 2, or a vector of pairs. Default is "all"
read_count <- function(dt, annotation, cds_m5 = 15, cds_m3 = 33, frame_ = "all", utr3_p = Inf) {
  require(data.table)
  dt <- merge(dt, annotation[, c("transcript", "nis_pos")], by = "transcript", all.x = TRUE)
  
  # Subset data
  message("Subsetting data")
  if (frame_ == "all") {
    message("\tFrame: 0, 1, 2")
    cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3, ]
    utr3_sub <- dt[psite_region == "3utr", ]
  } else if (frame_ < 0 | frame_ > 2) {
    message("\tInvalid frame number. Please put either 0, 1, 2, or all")
  } else {
    message(paste0("\tFrame: ", frame_))
    cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3 & frame %in% frame_, ]
    utr3_sub <- dt[psite_region == "3utr" & frame %in% frame_, ]
  }
  message(paste0("\t3'UTR constraint: Will count reads up until the first in-frame stop codon or up until ", utr3_p, " nt"))
  ext_sub <-  utr3_sub[psite_from_stop < nis_pos, ]
  utr3p_sub <-  ext_sub[psite_from_stop <= utr3_p, ]
  
  # Count reads
  message("Counting reads")
  cds_tab <- cds_sub[, list(c_cds = .N), by = list(transcript)]       # Read count in CDS
  utr3_tab <- utr3_sub[, list(c_utr3 = .N), by = list(transcript)]    # Read count in entire 3'-UTR
  ext_tab <- ext_sub[, list(c_ext = .N), by = list(transcript)]       # Read count in extension
  utr3p_tab <- utr3p_sub[, list(c_utr3p = .N), by = list(transcript)] # Read count in fixed length of 3'-UTR
  
  # Combine data
  message("Combining read count tables")
  cu_tab <- Reduce(function(df1, df2) merge(df1, df2, by = "transcript", all.x = TRUE), list(annotation, cds_tab, utr3_tab, ext_tab, utr3p_tab))
  message("Replacing NA count with 0")
  cu_tab[is.na(cu_tab)] <- 0
  
  # Calculate length of extension and distal 3'UTR
  cu_tab$l_ext <- pmin(cu_tab$nis_pos-1, cu_tab$l_utr3)
  cu_tab$l_utr3p <- pmin(cu_tab$nis_pos-1, cu_tab$l_utr3, utr3_p)
  
  cu_tab <- data.table(mutate_if(cu_tab, is.character, as.factor))
  message("Done\n")
  return(cu_tab)
}
#--------------------------------------------------
# Pool replicate read count 
  # data_list: list of tables (outputs of read count function) to be pooled
composite <- function(data_list) {
  require(data.table)
  data_comp <- rbindlist(data_list, use.names = TRUE)[,lapply(.SD, sum), by = list(transcript, stop_codon, l_tr, l_utr5, l_cds, l_utr3, nis_pos, nis_stop, l_ext, l_utr3p)]
  return(data_comp)
}
#--------------------------------------------------
# Calculate RPKM and readthrough efficiency
  # dt: data table containing read count information (output of read_count or composite function)
  # cds_m5: the number of nucleotides from the start codon (5' end of CDS, start codon inclusive) to EXCLUDE from the CDS count
  # cds_m3: the number of nucleotides from the stop codon (3' end of CDS, stop codon exclusive) to EXCLUDE from the CDS count
rt_efficiency <- function(dt, cds_m5 = 15, cds_m3 = 33) {
  require(data.table)
  # Calculate RPKM
  lib_size <- (sum(dt$c_cds) + sum(dt$c_utr3))/10^6
  dt$rpkm_cds <- dt$c_cds/(lib_size * ((dt$l_cds - cds_m5 - cds_m3)/10^3))
  dt$rpkm_utr3 <- dt$c_utr3/(lib_size * ((dt$l_utr3)/10^3))
  dt$rpkm_ext <- dt$c_ext/(lib_size * ((dt$l_ext)/10^3))
  dt$rpkm_utr3p <- dt$c_utr3p/(lib_size * ((dt$l_utr3p)/10^3))
  # Calculate readthrough efficiency
  dt$rte_utr3 <- dt$rpkm_utr3/dt$rpkm_cds
  dt$rte_ext <- dt$rpkm_ext/dt$rpkm_cds
  dt$rte_utr3p <- dt$rpkm_utr3p/dt$rpkm_cds
  return(dt)
}
#--------------------------------------------------