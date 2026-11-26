library(dplyr)
library(GenomicRanges)
library(readr)

getwd()
setwd("C:/Users/anima/Documents/KELLEHER/MOTIF PROJECT/FIMO align")

FIMO_OUTPUT <- read_tsv("FIMO_OUTPUT.tsv", comment = "#")
FIMO_OUTPUT$sequence_name <- "chr2L"

FIMO_unique_IDs<- pull(FIMO_OUTPUT, motif_alt_id) %>% sort() %>% unique()
df_FIMO_unique_IDs <-data.frame(Items = unlist(FIMO_unique_IDs))

# ROI: 12250870  12269023
range(FIMO_OUTPUT$start, na.rm = TRUE) # [1] 12251050  12268850
length(unique(FIMO_OUTPUT$motif_alt_id)) # 201 unique motif_ids



# (1) grouping hits by motif ID & eliminating overlap. -------------
FIMO_GR <- GRanges(
  seqnames     = FIMO_OUTPUT$sequence_name,
  ranges       = IRanges(start = FIMO_OUTPUT$start, end = FIMO_OUTPUT$stop),
  strand       = FIMO_OUTPUT$strand,
  motif_id     = FIMO_OUTPUT$motif_id,
  motif_alt_id = FIMO_OUTPUT$motif_alt_id,
  score        = FIMO_OUTPUT$score,
  pvalue       = FIMO_OUTPUT$`p-value`,
  qvalue       = FIMO_OUTPUT$`q-value`,
  matched_seq  = FIMO_OUTPUT$matched_sequence
)

FIMO_DF <- as.data.frame(FIMO_GR) 
FIMO_GR_split <- split(FIMO_GR, FIMO_GR$motif_alt_id) 
reduced_list    <- lapply(FIMO_GR_split, GenomicRanges::reduce)
FIMO_GR_reduced <- unlist(GRangesList(reduced_list), use.names = FALSE)

motif_counts <- lengths(reduced_list)
FIMO_GR_reduced$motif_alt_id <- rep(names(reduced_list), motif_counts)

length(FIMO_GR)          # before: 1437 hits
length(FIMO_GR_reduced)  # no repeats: 1092 unique hits



# (2) Extracting reduced hits -------------
hits <- findOverlaps(FIMO_GR_reduced, FIMO_GR)
qval_by_query <- tapply(FIMO_GR$qvalue[subjectHits(hits)], queryHits(hits), min)
best_qval <- setNames(as.numeric(qval_by_query), names(qval_by_query))

FIMO_GR_reduced$best_qvalue <- unname(best_qval[as.character(seq_along(FIMO_GR_reduced))])
fimo_reduced_df <- as.data.frame(FIMO_GR_reduced)



# (3) Extracting ATAC seq peak coordinates -------------
RRs <- read.csv("ATAC seq data/regulatory_regions_only.csv")

atac_peaks <- RRs %>%
  dplyr::filter(Type == "ATAC Peak") %>%
  dplyr::group_by(Element) %>%
  dplyr::summarise(start = min(x), end = max(x), .groups = "drop") %>%
  dplyr::distinct(start, end, .keep_all = TRUE)   # no duplicates

atac_peaks_gr <- GRanges(
  seqnames = "chr2L", 
  ranges = IRanges(atac_peaks$start, atac_peaks$end),
  peak_id = atac_peaks$Element
)


# (4) Position intersect (T/F and coordinates)-------------
FIMO_GR_reduced$in_atac_peak <- overlapsAny(FIMO_GR_reduced, atac_peaks_gr, type = "within") 
table(FIMO_GR_reduced$in_atac_peak) #   581 FALSE     511 TRUE


ov <- findOverlaps(FIMO_GR_reduced, atac_peaks_gr, type = "within")
fimo_in_peaks <- FIMO_GR_reduced[queryHits(ov)]
fimo_in_peaks$atac_peak <- atac_peaks_gr$peak_id[subjectHits(ov)]
length(fimo_in_peaks)

fimo_in_peaks_df <- as.data.frame(fimo_in_peaks)
pull(fimo_in_peaks_df, motif_alt_id) %>% sort() %>% unique() # 178 unique motifs


# (5) Examining hits in all peaks -------------
nhits_in_peaks_df <- fimo_in_peaks_df %>%
  dplyr::group_by(atac_peak) %>%
  dplyr::summarise(n_hits = dplyr::n(),
                   n_distinct_motifs = dplyr::n_distinct(motif_alt_id),
                   best_q = min(best_qvalue),
                   .groups = "drop")

atac_peaks$n_hits <- nhits_in_peaks_df$n_hits
atac_peaks$n_unique <- nhits_in_peaks_df$n_distinct_motifs
atac_peaks$best_q <- nhits_in_peaks_df$best_q

summary(atac_peaks)


# (6) Examining hits in significant peak -------------
#     significant peak ID: 12254805
fimo_in_sig_df <- fimo_in_peaks_df %>% 
  dplyr::filter(atac_peak == 12254805) %>%
  dplyr::select(-strand, -in_atac_peak) %>%
  dplyr::distinct(motif_alt_id, .keep_all = TRUE)

fimo_in_sig_df %>%
  dplyr::arrange((abs(best_qvalue))) %>%
  dplyr::select(dplyr::any_of(c("start", "end", "motif_alt_id", "best_qvalue", "atac_peak"))) %>%
  head(20)
  
write.csv(fimo_in_sig_df, "sig_hits_11-24.csv", row.names = FALSE)
