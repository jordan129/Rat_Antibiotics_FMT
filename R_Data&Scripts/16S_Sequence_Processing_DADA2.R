# References:
# https://benjjneb.github.io/dada2/tutorial.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

### Load all packages
install.packages("pacman")
library("pacman")
pacman::p_load(ggplot2, dada2, phyloseq, gridExtra, dplyr, DECIPHER, phangorn, knitr)

# dada2 for sequence processing; DECIPHER to performs sequence alignment
# phangorn for phylogenetic tree generation; phyloseq for data visualisation and analysis

setwd("D:/dada2/") # working folder to store outputs
getwd()

path <- "D:/16S_data" # CHANGE to the directory containing the fastq files and silva data files (provided in the R_Data&Scripts folder in Github)
list.files(path)

### STEP 1 ### 
## Read sample names

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

### STEP 2 ###
## Inspect read quality profiles

#Visualize the quality profiles of the forward reads:
plotQualityProfile(fnFs[1:6])

#Visualize the quality profiles of the reverse reads:
plotQualityProfile(fnRs[1:6])

### STEP 3 ###
## Filter annd trim (it takes 5 seconds per sample)

# Make filenames for the filtered fastq files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19, 19), truncLen=c(250,185),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

save.image("DADA2_output_1.Rdata")

#If you want to speed up downstream computation, consider tightening maxEE. 
#If too few reads are passing the filter, consider relaxing maxEE, 
#perhaps especially on the reverse reads (eg. maxEE=c(2,5)), 
#and reducing the truncLen to remove low quality tails. 
#Remember you must maintain overlap after truncation (truncLen) in order to merge them later.

#trimLeft=c(19, 19) is used to remove the primers, if you have done it, delete this
#truncQ=2, #truncate reads after a quality score of 2 or less
#truncLen=130, #truncate after 130 bases
#trimLeft=10, #remove 10 bases off the 5’ end of the sequence
#maxN=0, #Don’t allow any Ns in sequence
#maxEE=2, #A maximum number of expected errors
#rm.phix=TRUE, #Remove lingering PhiX (control DNA used in sequencing) as it is likely there is some.
#mutithread cannot work on windows, so it doesnot matter if you select true or false

### STEP 4 ### 
## Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

save.image("DADA2_output_2.Rdata")

### STEP 5 ###
## Dereplication and inference
# Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence
# Dereplication substantially reduces computation time by eliminating redundant comparisons.

##Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

save.image("DADA2_output_3.Rdata")

### STEP 6 ###
## Merge paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

save.image("DADA2_output_4.Rdata")

### STEP 7 ###
## Construct sequence table

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)))

#You can remove non-target-length sequences with base R manipulations of the sequence table
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(272,366)] # change according to histogram
seqtab <- seqtab2
save.image("DADA2_output_5.Rdata")

### STEP 8 ###
## Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

##Sanity check
#Track reads through the pipeline, This is good to report in your methods/results
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "track_reads.csv")

save.image("DADA2_output_6.Rdata")

### STEP 9 ###
## Assign taxonomy using Silva database

taxa_silva <- assignTaxonomy(seqtab.nochim, "D:/16S_data/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa_silva_species <- addSpecies(taxa_silva, "D:/16S_data/silva_species_assignment_v138.fa.gz")

taxa_silva_one_direction <- assignTaxonomy(seqtab.nochim, "D:/16S_data/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa_silva_species_one_direction <- addSpecies(taxa_silva_one_direction, "D:/16S_data/silva_species_assignment_v138.fa.gz")

save.image("DADA2_output_8.Rdata")

#Inspect the taxonomic assignments:
taxa.print <- taxa_silva # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences 
# are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. 
# Tell dada2 to try the reverse-complement orientation and see if this fixes the assignments:
# taxa <- assignTaxonomy(seqtab.nochim, "/silva_nr_v128_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)

write.csv(seqtab.nochim, "USV_counts_silva.csv")
write.csv(taxa_silva, "taxa_silva.csv")
write.csv(taxa_silva_species, "taxa_silva_species.csv")
write.csv(taxa_silva_species_one_direction, "taxa_silva_species_one_direction.csv")

# here we do both one and two direction. When we get both csv, compare and combine. Some bacteria
# will only be assigned in one or two direction file, so you need to compare and combine to get more.

save.image("DADA2_output_9.Rdata")