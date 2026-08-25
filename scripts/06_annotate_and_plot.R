#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Annotate per-region G4 scores with gene symbols and plot:
#   1. anchor meta-profile (from 05_g4_metaprofile.py output)
#   2. Distributions of per-gene sum and max G4 scores
#
# Usage:
#   Rscript 06_annotate_and_plot.R \
#       <scores.bed>        # output of 04_g4hunter_to_bed.py \
#       <profile.tab>       # output of 05_g4_metaprofile.py \
#       <annotation.tab>    # gene.id <tab> ... <tab> gene.symbol (see README) \
#       <output_prefix>
#
# Author: Nejc Haberman
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 06_annotate_and_plot.R <scores.bed> <profile.tab> ",
       "<annotation.tab> <output_prefix>")
}
scores_file     <- args[1]
profile_file    <- args[2]
annotation_file <- args[3]
out_prefix      <- args[4]

# --- per-region scores (BED6+3 from 04_g4hunter_to_bed.py) ---------
scores <- read.table(scores_file, header = FALSE, sep = "\t",
                     col.names = c("chr", "start", "end", "gene.id",
                                   "G4.score.sum", "strand",
                                   "G4.score.max",
                                   "G4.antisense.score.sum",
                                   "G4.antisense.score.max"))

# --- gene annotation -----------------------------------------------
# expected columns: gene.id chr strand start end gene.id2 gene.symbol
annotation <- read.table(annotation_file, header = FALSE, sep = "\t")
colnames(annotation) <- c("gene.id", "chr", "strand", "start", "end",
                          "gene.id2", "gene.symbol")
annotation <- annotation[, c("gene.id", "gene.symbol")]

scores.annotated <- merge(scores, annotation, by = "gene.id", all.x = TRUE)
write.table(scores.annotated,
            paste0(out_prefix, ".scores.annotated.tab"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# --- anchor meta-profile ----------------------------------------------
profile <- read.table(profile_file, header = TRUE, sep = "\t")

p.profile <- ggplot(profile, aes(x = position)) +
  theme_cowplot() +
  background_grid(major = "xy", minor = "none") +
  geom_line(aes(y = sense_score_sum, colour = "sense"), linewidth = 0.4) +
  geom_line(aes(y = antisense_score_sum, colour = "antisense"),
            linewidth = 0.4) +
  scale_colour_manual(name = NULL,
                      values = c(sense = "#42476E", antisense = "#D94D3D")) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "G4Hunter score profile around the region anchor",
       x = "position relative to anchor (bp)",
       y = "summed G4Hunter score") +
  theme(text = element_text(size = 10))

ggsave(paste0(out_prefix, ".metaprofile.pdf"), p.profile,
       width = 6, height = 4)

# --- score distributions -------------------------------------------
plot_dist <- function(df, column, xlab_txt) {
  ggplot(df, aes(x = .data[[column]])) +
    theme_cowplot() +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.5,
                   colour = "black", fill = "white") +
    geom_density(alpha = 0.5, colour = "#FF6666") +
    labs(x = xlab_txt, y = "density")
}

nonzero <- scores.annotated[scores.annotated$G4.score.sum > 0, ]
p.sum <- plot_dist(nonzero, "G4.score.sum", "per-gene G4 score sum")
p.max <- plot_dist(nonzero, "G4.score.max", "per-gene G4 score max")

ggsave(paste0(out_prefix, ".score_sum_distribution.pdf"), p.sum,
       width = 5, height = 4)
ggsave(paste0(out_prefix, ".score_max_distribution.pdf"), p.max,
       width = 5, height = 4)

message("Written: ", out_prefix, ".scores.annotated.tab + 3 PDFs")
