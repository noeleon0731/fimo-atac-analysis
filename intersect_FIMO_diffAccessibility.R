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



# (2) Extracting Reduced hits -------------
hits <- findOverlaps(FIMO_GR_reduced, FIMO_GR)
qval_by_query <- tapply(FIMO_GR$qvalue[subjectHits(hits)], queryHits(hits), min)
best_qval <- setNames(as.numeric(qval_by_query), names(qval_by_query))

FIMO_GR_reduced$best_qvalue <- unname(best_qval[as.character(seq_along(FIMO_GR_reduced))])
fimo_reduced_df <- as.data.frame(FIMO_GR_reduced)



# (3) Extracting ATAC seq peak coverage -------------
library(IRanges)

A4_atac <- read.table("ATAC seq data/ovary atac seq A4.txt", 
                      col.names = c("chrom", "x1", "x2", "y"))
A4_atac$genome <- "A4"
A5_atac <- read.table("ATAC seq data/ovary atac seq A5.txt", 
                      col.names = c("chrom", "x1", "x2", "y"))
A5_atac$genome <- "A5"

ALL_atac <- rbind(A4_atac, A5_atac)
ALL_atac <- ALL_atac %>% filter(x1 >= 12250870, x1 <= 12269023)


atac_diff <- A4_atac %>%
  inner_join(A5_atac, by = c("chrom", "x1", "x2"), suffix = c("_A4", "_A5")) %>%
  filter(x1 >= 12250870, x1 <= 12269023) %>%
  arrange(x1)

call_peaks <- function(df, y_col, cutoff_quantile = 0.75, window = 150) {
  df <- df %>% arrange(x1)
  y  <- df[[y_col]]
  is_max <- y > lag(y, default = -Inf) & y > lead(y, default = -Inf)
  cutoff <- quantile(y, cutoff_quantile)
  summits <- df$x1[is_max & y >= cutoff]
  gr <- GRanges("chr2L", IRanges(start = summits - window, end = summits + window))
  GenomicRanges::reduce(gr)
}

peaks_A4 <- call_peaks(atac_diff, "y_A4")
peaks_A5 <- call_peaks(atac_diff, "y_A5")

cov_A4 <- GenomicRanges::coverage(
  GRanges("chr2L", IRanges(atac_diff$x1 + 1, atac_diff$x2), score = atac_diff$y_A4),
  weight = "score")[["chr2L"]]

cov_A5 <- GenomicRanges::coverage(
  GRanges("chr2L", IRanges(atac_diff$x1 + 1, atac_diff$x2), score = atac_diff$y_A5),
  weight = "score")[["chr2L"]]



# (4) position intersecting hits with A4 and A5 coverage -------------
fimo_reduced_df$atac_A4_mean <- viewMeans(Views(cov_A4, start = start(FIMO_GR_reduced), end = end(FIMO_GR_reduced)))
fimo_reduced_df$atac_A5_mean <- viewMeans(Views(cov_A5, start = start(FIMO_GR_reduced), end = end(FIMO_GR_reduced)))

fimo_reduced_df$in_peak_A4 <- overlapsAny(FIMO_GR_reduced, peaks_A4)
fimo_reduced_df$in_peak_A5 <- overlapsAny(FIMO_GR_reduced, peaks_A5)

# (5) sort by lowest qvalue -------------
fimo_reduced_df %>%
  arrange(abs(best_qvalue)) %>%
  select(start, end, motif_alt_id, atac_A4_mean, atac_A5_mean, best_qvalue) %>%
  head(15)