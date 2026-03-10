# -------------------------------------------------------------------------
#   Author:       Priscilla A San Juan
#   Date edited:  2024-09-05
#   Topic:        Microbiome pipeline for possum microbiome
#   Bioinformatics (cleaning sequences and preparing phyloseq object)
# -------------------------------------------------------------------------

# Load packages  ----------------------------------------------------------
library(dada2)
library(phyloseq)
library(phangorn)
library(DECIPHER)
library(dplyr)
library(ggplot2)
library(vegan)
library(DESeq2)


# Specify forward & reverse read fastqs -----------------------------------
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# split filenames with '_', pick first character string
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) 
sample.names

# Data quality metrics ----------------------------------------------------
# looking at forward fastq files, Phred scores
plotQualityProfile(fnFs[1:9]) # recommended to trim last 10 nucleotides (@ ~311)
plotQualityProfile(fnRs[1:9]) # reverse reads, trim @ ~271

# Filter and trim ---------------------------------------------------------
# Filtering reads based on sequence quality scores
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # creating filtered folder
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Pulling sample names from filtered fasta files
names(filtFs) <- sample.names 
names(filtRs) <- sample.names

# set your primer sequence
FWD <- "NNNNNNGTGYCAGCMGCCGCGGTAA"
REV <- "NNNNNNGGACTACNVGGGTWTCTAAT"
trimLeft = c(FWD,REV)

# Use known primer sequences to trim from your amplicon sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE,trimLeft = c(25,26))

existing_Ffiles <- filtFs[file.exists(filtFs)]
nonexistent_Ffiles <- filtFs[!file.exists(filtFs)]
existing_Rfiles <- filtRs[file.exists(filtRs)]
nonexistent_Rfiles <- filtRs[!file.exists(filtRs)]

# Learning Error rates ----------------------------------------------------
errF <- learnErrors(existing_Ffiles, multithread=TRUE)
errR <- learnErrors(existing_Rfiles, multithread=TRUE)

# Plotting out the errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample inference from forward reads 
dadaFs <- dada(existing_Ffiles, err=errF, multithread=TRUE)

# Sample inference from reverse reads
dadaRs <- dada(existing_Rfiles, err=errR, multithread=TRUE)

# Inspect the dada objects produced above
dadaFs[[1]]
mergers <- mergePairs(dadaFs, existing_Ffiles, dadaRs, existing_Rfiles, verbose = TRUE)

# View mergers 
head(mergers[[1]])

# Making your ASV table ---------------------------------------------------
seqtab <- makeSequenceTable(mergers)

# Length of ASV's and their frequencies
table(nchar(getSequences(seqtab)))

# Identifying and removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",multithread=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Tracks reads through the pipeline
getN <- function(x) sum(getUniques(x))
track<- cbind(out, sapply(dadaFs,getN),sapply(dadaRs, getN), sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
write.csv(track, "Sequencing Statistics_16S.csv")

# -------------------------------------------------------------------------
# Assign taxonomy for 16S using native bayesian classifier built into dada2, against the silva database
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Attach taxonomic/species labels
taxa2 <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

# Check your taxonomy assignment
taxa.print <- taxa2 
dim(taxa.print)
head(taxa.print)

# Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#  Giving taxonomy table corresponding names as above (ASV_1, ASV_2...)
row.names(taxa.print) <- sub(">", "", asv_headers)
write.table(taxa.print, "ASVs_named_correctly.tsv", sep="\t", quote=F, col.names=NA)
