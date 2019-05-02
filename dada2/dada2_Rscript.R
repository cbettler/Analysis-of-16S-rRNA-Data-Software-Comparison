#Install Dada2

#source("http://bioconductor.org/biocLite.R")
#biocLite(suppressUpdates = FALSE)
#biocLite("ShortRead", suppressUpdates = FALSE)
#biocLite("devtools")

#install.packages("c:/Users/okoro/OneDrive/Desktop/BIOI 500/Project/dada2-1.10", repos = NULL,type = "source", dependencies = c("Depends", "Suggests","Imports"))

#library("devtools")
#devtools::install_github("benjjneb/dada2")

ptm <- proc.time()

library(dada2)
packageVersion("dada2") # Version 1.11.3
path <- "c:/Users/okoro/OneDrive/Desktop/BIOI 500/Project/Miseq16srRNA_Dataset/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#Examine quality profiles of forward and reverse reads

#plotQualityProfile(fnFs[1:6])
#plotQualityProfile(fnRs[1:6])

#Perform filtering and trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
#Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, matchIDs = TRUE, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
#head(out)

##Examine quality profiles of filtered forward and reverse reads
#pdf("/home/paul/metstry1/Second/qplotTrimmedF.pdf")
#plotQualityProfile(filtFs[1:6])
#dev.off()
#pdf("/home/paul/metstry1/Second/qplotTrimmedR.pdf")
#plotQualityProfile(filtRs[1:6])
#dev.off()

#Learn the Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plot learned error rates
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

#Depreplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
#We are now ready to apply the core sequence-variant inference algorithm to the dereplicated data.
#Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Merge Paired reads
#Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#mergePairs(ddF, drpF, ddR, drpR, minOverlap=4)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

#Construct the sequence table (ASV or OTU table)
#The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants
#construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods
seqtab <- makeSequenceTable(mergers)
#dim(seqtab)
# Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))
#Remove Chimeric Sequence
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)
#write out the ASV(OTU) table
#write.table(seqtab.nochim, file = "/home/paul/Project/dada2_results/OTU_table.txt")
#asv <- read.table(file = "/home/paul/Project/dada2_results/OTU_table.txt", header = TRUE) #Just to check the OTU is looking good

#Track reads through the pipeline
#As a final check of our progress, will look at the number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
#head(track)

#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "c:/Users/okoro/OneDrive/Desktop/BIOI 500/Project/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
#Assign Specie
taxa <- addSpecies(taxa, "c:/Users/okoro/OneDrive/Desktop/BIOI 500/Project/silva_species_assignment_v128.fa.gz")
#write out the taxonomy with species
#write.table(taxa, file = "/home/paul/Project/dada2_results/taxonomy.txt")

#Lets inspect the taxonomic assignments
#taxa.print <- taxa # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
#head(taxa.print) #Prints the taxonomy with species


#Mock Validation
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
#DADA2 inferred 20 sample sequences present in the Mock community.
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
#Of those, 20 were exact matches to the expected reference sequences.

library(phyloseq)
library(ggplot2)

#Use Phyloseq to generate taxa abundance
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


#ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_names(samples.out), tax_table(taxa))
#ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

#Check abundance of taxa families
top200 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:200]
ps.top200 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top200 <- prune_taxa(top200, ps.top200)
plot_bar(ps.top200, fill="Genus")

proc.time() - ptm
