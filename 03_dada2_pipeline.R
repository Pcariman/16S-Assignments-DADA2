# ------------------------------------------------------------------------------
# Script: 03_dada2_pipeline.R
# Description: DADA2 pipeline for quality filtering, ASV inference, taxonomic
#              assignment, and Phyloseq object creation from 16S rRNA data.
# Project: 16S-Assignments-DADA2
# ------------------------------------------------------------------------------

##### 1. Load required packages #####
library(dada2)         # Core denoising and ASV inference
library(phyloseq)      # Integration of OTU table, taxonomy, tree and metadata
library(Biostrings)    # Handling of DNA sequences
library(DECIPHER)      # Sequence alignment
library(phangorn)      # Phylogenetic tree construction

##### 2. Define paths and seed #####
set.seed(12)

# Input/output folders relative to project root
input_dir     <- "raw_data/2_trimmed_cutadapt/trimmomatic/paired"
filtered_dir  <- file.path(input_dir, "filtered")
ref_db        <- "reference/silva_nr99_v138.1_train_set.fa"
metadata_file <- "metadata/metadata.txt"
output_rds    <- "results/RESULTADOS_ASIGNACION_16S.rds"

# Create output folders if they don't exist
dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)

##### 3. List input FASTQ files #####
fnFs <- sort(list.files(input_dir, pattern = "_R1_trimmomatic_paired.fastq", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern = "_R2_trimmomatic_paired.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_trimmomatic_paired"), `[`, 1)

##### 4. (Optional) Quality profile plots #####
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])

##### 5. Filter and trim #####
filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(260, 260),
                     maxN = 0,
                     maxEE = c(2, 2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)

##### 6. Learn error rates #####
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

##### 7. Denoise, merge and build sequence table #####
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

##### 8. Track reads through the pipeline #####
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##### 9. Taxonomic assignment #####
taxa <- assignTaxonomy(seqtab.nochim, refFasta = ref_db,
                       multithread = TRUE,
                       tryRC = FALSE,
                       minBoot = 90,
                       outputBootstraps = TRUE)

##### 10. Phylogenetic tree construction #####
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor = NA)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit <- pml(treeNJ, data = phang.align)
fitGTR <- optim.pml(update(fit, k = 4, inv = 0.2),
                    model = "GTR",
                    optInv = TRUE,
                    optGamma = TRUE,
                    rearrangement = "NNI",
                    control = pml.control(trace = 0))

##### 11. Import metadata and create phyloseq object #####
meta <- read.delim(metadata_file, header = TRUE, row.names = 1)
physeq <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(meta),
  tax_table(taxa),
  phy_tree(fitGTR$tree)
)

##### 12. Save results #####
saveRDS(physeq, file = output_rds)
