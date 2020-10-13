#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(tximport)
library(DESeq2)
library(rtracklayer)
library(data.table)
library(dplyr)
library(tidyr)

#files and directories
gff <- snakemake@input[["gff"]]
samples_file <- snakemake@input[["samples_file"]]
salmon_quant_dir <- snakemake@input[["salmon_quant_dir"]]
dds_file <- snakemake@output[["dds_file"]]
deseq2_result_file <- snakemake@output[["deseq2_result_file"]]

samples <- read.table(samples_file, header = TRUE)
files <- file.path(salmon_quant_dir, "quant", samples$sample, "quant.sf")

# read gff
gr <- import.gff3(gff, feature.type = c("exons", "CDS", "mRNA", "gene"))

# extract a data.frame of tx to gene
mrnas <- gr[gr$type == "mRNA",]
mrna_dt <- as.data.table(mcols(mrnas))
tx2gene <- data.frame(mrna_dt[, .(TXNAME = Name, GENEID = as.character(gene))])

names(files) <- paste0(samples$sample)

#import quant files
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

#generate & save DESeq2 object
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~group)
saveRDS(dds, dds_file)

# run DESeq2
dds <- DESeq(dds)
res <- results(dds, tidy = TRUE)
res <- res[order(res$padj),]

# write results table
write.table(res, deseq2_result_file, append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)

#join emapper annotations with DEseq2 results
#files and directories
emapper_annotation_file <- snakemake@input[["emapper_annotation_file"]]
deseq2_result_with_annotation_file <- snakemake@output[["deseq2_result_with_annotation_file"]]

#read gff file
feature_types <- c("exon")
my_gff <- import.gff(gff, feature.type = feature_types)
my_gff_dt <- unique(as.data.table(my_gff)[, .(transcript_id, gene, product)])

# read results and add column names
my_results <- cbind(rownames(res), data.frame(res, row.names=NULL))
setnames(my_results, old = c('row'), new = c('gene'))
my_results_dt <- as.data.table(my_results)

#join results with gff
setkey(my_gff_dt, gene)
setkey(my_results_dt, gene)
my_join <- my_gff_dt[my_results_dt, nomatch=0, all = TRUE]
my_join_dt <- as.data.table(my_join)

# read and prep. emapper annotations
emapper_annot <- read.delim(emapper_annotation_file, header=FALSE, comment.char="#")
sep_emapper_annot <- separate(emapper_annot, "V1", c('transcript', 'V1'), sep = ".p")
prepd_annot <- subset(sep_emapper_annot, , -c(V1, V2, V3, V4, V5, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21))
prepd_annot <- prepd_annot[, c(1, 2, 4, 3)]
colnames(prepd_annot) <- c('transcript', 'emapper_gene', 'emapper_description', 'emapper_go')

# merge results with emapper annotations
merged <- merge(my_join_dt, prepd_annot, by.x = "transcript_id", by.y = "transcript", all = TRUE)
setnames(merged, old = c('gene','product'), new = c('gff_gene','gff_product'))
merged <-merged[order(merged$padj),]
merged = merged[!duplicated(merged$gff_gene),]

# write results with annotations table
write.table(merged, deseq2_result_with_annotation_file, append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)

# log
sessionInfo()
