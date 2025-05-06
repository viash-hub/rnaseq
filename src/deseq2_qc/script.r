## VIASH START
par <- list(
  counts = "testData/unit_test_resources/counts.tsv",
  id_col = 1,
  sample_suffix = "",
  outprefix = "deseq2",
  count_col = 2,
  deseq2_output = "deseq2",
  pca_multiqc = "pca.vals_mqc.tsv",
  dists_multiqc = "sample.dists_mqc.tsv",
  vst = FALSE,
  outdir = '.'
)
meta <- list(
  resources_dir = "src/deseq2_qc"
)
## VIASH END

# REQUIREMENTS

## PCA, HEATMAP AND SCATTERPLOTS FOR SAMPLES IN COUNTS FILE
## - SAMPLE NAMES HAVE TO END IN e.g. "_R1" REPRESENTING REPLICATE ID. LAST 3 CHARACTERS OF SAMPLE NAME WILL BE TRIMMED TO OBTAIN GROUP ID FOR DESEQ2 COMPARISONS.

# LOAD LIBRARIES
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(stringr)

if (file.exists(par$outdir) == FALSE) {
  dir.create(par$outdir, recursive = TRUE)
}

# READ IN COUNTS FILE
count_table <- read.delim(file = par$counts, header = TRUE, row.names = NULL)
rownames(count_table) <- count_table[, par$id_col]
count_table <- count_table[, par$count_col:ncol(count_table), drop = FALSE]
colnames(count_table) <- gsub(par$sample_suffix, "", colnames(count_table))
colnames(count_table) <- gsub(pattern = '\\.$', replacement = '', colnames(count_table))

# RUN DESEQ2
samples_vec <- colnames(count_table)
name_components <- strsplit(samples_vec, "_")
n_components <- length(name_components[[1]])
decompose <- n_components != 1 && all(sapply(name_components, length) == n_components)
coldata <- data.frame(samples_vec, sample = samples_vec, row.names = 1)
if (decompose) {
  groupings <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
  n_distinct <- sapply(groupings, function(grp) length(unique(grp)))
  groupings <- groupings[n_distinct != 1 & n_distinct != length(samples_vec)]
  if (ncol(groupings) != 0) {
    names(groupings) <- paste0("Group", 1:ncol(groupings))
    coldata <- cbind(coldata, groupings)
  } else {
    decompose <- FALSE
  }
}

DDSFile <- paste(par$outdir, "/", par$outprefix, ".dds.RData", sep = "")

counts <- count_table[, samples_vec, drop = FALSE]
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~1)
dds <- estimateSizeFactors(dds)

# No point if only one sample, or one gene
if (min(dim(count_table)) <= 1)  {
  save(dds, file = DDSFile)
  saveRDS(dds, file = sub("\\.dds\\.RData$", ".rds", DDSFile))
  warning("Not enough samples or genes in counts file for PCA.", call. = FALSE)
  quit(save = "no", status = 0, runLast = FALSE)
}

if (!par$vst) {
  vst_name <- "rlog"
  rld <- rlog(dds)
} else {
  vst_name <- "vst"
  rld <- varianceStabilizingTransformation(dds)
}

assay(dds, vst_name) <- assay(rld)
save(dds, file = DDSFile)
saveRDS(dds, file = sub("\\.dds\\.RData$", ".rds", DDSFile))

# PLOT QC

##' PCA pre-processeor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on.
##'
##' @param object The DESeq2DataSet object.
##' @param ntop number of top genes to use for principla components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly

plotPCA_vst <- function(object,  ntop = 500, assay = length(assays(object))) {
  rv <- rowVars(assay(object, assay), useNames=FALSE)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object, assay)[select, ]), center = TRUE, scale = FALSE)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  df <- cbind(as.data.frame(colData(object)), pca$x)
  # Order points so extreme samples are more likely to get label
  ord <- order(abs(rank(df$PC1) - median(df$PC1)), abs(rank(df$PC2) - median(df$PC2)))
  df <- df[ord, ]
  attr(df, "percentVar") <- data.frame(PC = seq(along = percentVar), percentVar = 100 * percentVar)
  return(df)
}

PlotFile <- paste(par$outdir, "/", par$outprefix, ".plots.pdf", sep = "")

pdf(file = PlotFile, onefile = TRUE, width = 7, height = 7)

## PCA
ntop <- c(500, Inf)
for (n_top_var in ntop) {
  pca_data <- plotPCA_vst(dds, assay = vst_name, ntop = n_top_var)
  percentVar <- round(attr(pca_data, "percentVar")$percentVar)
  plot_subtitle <- ifelse(n_top_var == Inf, "All genes", paste("Top", n_top_var, "genes"))
  pl <- ggplot(pca_data, aes(PC1, PC2, label = paste0(" ", sample, " "))) +
    geom_point() +
    geom_text(check_overlap = TRUE, vjust = 0.5, hjust="inward") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  print(pl)

  if (decompose) {
    pc_names <- paste0("PC", attr(pca_data, "percentVar")$PC)
    long_pc <- reshape(pca_data, varying=pc_names, direction="long", sep="", timevar="component", idvar="pcrow")
    long_pc <- subset(long_pc, component <= 5)
    long_pc_grp <- reshape(long_pc, varying = names(groupings), direction = "long", sep = "", timevar = "grouper")
    long_pc_grp <- subset(long_pc_grp, grouper <= 5)
    long_pc_grp$component <- paste("PC", long_pc_grp$component)
    long_pc_grp$grouper <- paste0(long_pc_grp$grouper, c("st", "nd", "rd", "th", "th")[long_pc_grp$grouper], " prefix")
    pl <- ggplot(long_pc_grp, aes(x = Group, y = PC)) +
      geom_point() +
      stat_summary(fun = mean, geom = "line", aes(group = 1)) +
      labs(x = NULL, y = NULL, subtitle = plot_subtitle, title = "PCs split by sample-name prefixes") +
      facet_grid(component ~ grouper, scales = "free_x") +
      scale_x_discrete(guide = guide_axis(n.dodge = 3))
    print(pl)
  }
} # at end of loop, we'll be using the user-defined ntop if any, else all genes

## WRITE PC1 vs PC2 VALUES TO FILE
pca_vals <- pca_data[, c("PC1", "PC2")]
colnames(pca_vals) <- paste0(colnames(pca_vals), ": ", percentVar[1:2], '% variance')
pca_vals <- cbind(sample = rownames(pca_vals), pca_vals)
pca_vals_file <- paste(par$outdir, "/", par$outprefix, ".pca_vals.txt", sep = "")
write.table(pca_vals, file = pca_vals_file,
            row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = TRUE)

## SAMPLE CORRELATION HEATMAP
sampleDists <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  main = paste("Euclidean distance between", vst_name, "of samples")
)

## WRITE SAMPLE DISTANCES TO FILE
sample_dist_file <- paste(par$outdir, "/", par$outprefix, ".sample.dists.txt", sep = "")
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),
            file = sample_dist_file, row.names = FALSE,
            col.names = TRUE, sep = "\t", quote = FALSE)
dev.off()

# SAVE SIZE FACTORS
SizeFactorsDir <- paste0(par$outdir, "/size_factors/")
if (file.exists(SizeFactorsDir) == FALSE) {
  dir.create(SizeFactorsDir, recursive = TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir, par$outprefix, ".size_factors.RData", sep = "")

normFactors <- sizeFactors(dds)
save(normFactors, file = NormFactorsFile)

for (name in names(sizeFactors(dds))) {
  sizeFactorFile <- paste(SizeFactorsDir, name, ".txt", sep = "")
  write(as.numeric(sizeFactors(dds)[name]), file = sizeFactorFile)
}

# R SESSION INFO
RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
a <- sessionInfo()
print(a)
sink()

# Prepare files for MultiQC

readLines(paste0(meta$resources_dir, "/deseq2_pca_header.txt")) |>
  stringr::str_replace(pattern = "#id: 'deseq2_pca'",
                      replace = paste0("#id: '", par$label, "_deseq2_pca'")) |>
  writeLines(con = "tmp.txt")

readLines(paste0("tmp.txt")) |>
  stringr::str_replace(pattern = "#section_name: 'DESeq2 PCA plot'",
                      replace = paste0("#section_name: 'DESeq2 PCA plot - '", par$label)) |>
  writeLines(con = "tmp.txt")

system2("cat", args = paste0("tmp.txt ", pca_vals_file), stdout = par$pca_multiqc)

readLines(paste0(meta$resources_dir, "/deseq2_clustering_header.txt")) |>
  stringr::str_replace(pattern = "#id: 'deseq2_clustering'",
                       replace = paste0("#id: '", par$label, "_deseq2_clustering'")) |>
  writeLines(con = "tmp.txt")

readLines(paste0("tmp.txt")) |>
  stringr::str_replace(pattern = "#section_name: 'DESeq2 sample similarity'",
                       replace = paste0("#section_name: 'DESeq2 sample similarity - '", par$label)) |>
  writeLines(con = "tmp.txt")

system2("cat", args = paste0("tmp.txt ", sample_dist_file), stdout = par$sample_dists_multiqc)
