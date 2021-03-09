#############################################
####==== Metabarcoding Bioinformatics ===####
####====         06.07.2020          ====####
#############################################

##Some Libraries

library(metabarTOAD)
library(dada2)
library(lulu)


#First we set up some folders
Folders()

#Now lets unzip the files
Unzip()

#Lets check if our files have been UnZipped
list.files("1.rawreads/", pattern=".fastq")

#looks good! 

#DADA2 requires no primers so lets have a look at doing that, we can't use the normal Primer strip argument here because it is for paired data

cutadapt <- "/Users/Luke/Bioinformatics/PATH/cutadapt"

ForwardPrimer <- "NNNNNNGGWACWGGWTGAACWGTWTAYCCYCC"

ReversePrimer <- "TANACYTCNGGRTGNCCRAARAAYCA"

forwards <- list.files("1.rawreads/", pattern="R1_001.fastq")

reverses <- list.files("1.rawreads/", pattern="R2_001.fastq")

for (sample in forwards){
  system2(cutadapt,args=paste0("-g ",ForwardPrimer," --trimmed-only -j 7 -o 3.strippedreads/",sample," 1.rawreads/",sample))
}

for (sample in reverses){
  system2(cutadapt,args=paste0("-g ",ReversePrimer," --trimmed-only -j 7 -o 3.strippedreads/",sample," 1.rawreads/",sample))
}


#Now lets start the DADA2 workflow

path <- "3.strippedreads"

# Forward and reverse fastq filenames have format: ****_R1_001.fastq and ****_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Lets look at the quality for forwards
plotQualityProfile(fnFs[1:2])
#...and reverses
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path("7.DADA2/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("7.DADA2/filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter expression

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,180),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     matchIDs=TRUE,compress=TRUE, multithread=7)


# Now lets learn the error rates

errF <- learnErrors(filtFs, multithread=7)

errR <- learnErrors(filtRs, multithread=7)

#lets visualise the error rates

plotErrors(errF, nominalQ=TRUE)

# Now lets dereplicate the samples, some of our files have dissapeared so we use a complex expression to restore them

derepFs <- derepFastq(list.files("7.DADA2/filtered", pattern="F_filt.fastq.gz", full.names = TRUE), verbose=TRUE)
derepRs <- derepFastq(list.files("7.DADA2/filtered", pattern="R_filt.fastq.gz", full.names = TRUE), verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="F_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="R_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)

# Now we can pply the DADA2 algorithm to the data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# We can now merge the seqs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Now we can construct a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=7, verbose=TRUE)
dim(seqtab.nochim)

#How many reads do we lose?
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- names(dadaFs)
head(track)






