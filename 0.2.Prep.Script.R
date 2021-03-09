#############################################
####==== Metabarcoding Decay Analysis ===####
####====         06.07.2020          ====####
#############################################

###Script 0 - QC & Prep Script###

####====0.0 Packages====####

library(data.table)
library(metabarTOAD)
library(vegan)
library(dplyr)
library(seqinr)

####====0.1 Parameters====####

#Set some variables 
minreads <- 2
items <- NULL

#Set the seed 
set.seed("123456")

#Read in metadata
metadat<-read.csv("metadata.csv") 

####====0.2 Taxonomy ====####

##taxonomy
#Tax <- fread(file = "taxonomy/YiMTB.tax.txt",sep="\t")
#names(Tax) <- c("OTUID","qlen","slen","qcov","qcovhsp","sseqid","bitscore","score","evalue","pctid","qstart","qend","start","send","staxid","sscinames")
#write.table(Tax,"taxonomy/RawTaxonomy.txt",row.names=FALSE,sep="\t", quote = FALSE)
#Assignments <- ParseTaxonomy(blastoutput = "taxonomy/RawTaxonomy.txt",lineages = "taxonomy/lineages-2019-08-14.csv.gz",lwrcovpct=80)
#write.table(Assignments,file="taxonomy/Parsedtax.csv",row.names=FALSE,sep=",", quote = FALSE)

Assignments <- read.csv("taxonomy/Parsedtax.csv")

####====0.3 OTU Table Cleaning  ====####
setwd("rawdata/")

rawdat <-read.csv(file="AllSamples.lulu.dada2.csv")
metadat$ID <- gsub("X","S",make.names(metadat$ID))
colnames(rawdat) <- gsub("X","S",colnames(rawdat))
  
  
  #Separate controls and samples
  samples <- rawdat[colnames(rawdat) %in% metadat$ID[metadat$Type=="sample" | metadat$Type=="control" ]]
  check <- rawdat[colnames(rawdat) %in% metadat$ID[metadat$Type=="control.check"]]
  controls <-  rawdat[colnames(rawdat) %in% metadat$ID[metadat$Type=="control.n"]]
  
  #Filter 1 - minimum number of reads for any ID
  samples[samples< minreads] <- 0
  samples <- samples[rowSums(samples) > 0,]
  
  #Filter 2 - within samples OTU must appear in more than one sample (this works because there are lots of reps per site and sample)
  filtersam <- samples
  filtersam[filtersam>0 ] <- 1
  filtersam <-filtersam[rowSums(filtersam) > 1,]
  samples <- samples[rownames(samples) %in% rownames(filtersam),]
  
  #Filter 3 -Make the maximum umber of reads for each OTU in the contam 100
  controlsCONTAM <- controls[rowSums(controls) > 0,]
  for (contamOTU in 1:length(controlsCONTAM[,1])){
    loopOTU <- row.names(controlsCONTAM[contamOTU,])
    loopMax <- 100
    if (any(is.na(samples[loopOTU,]))){next}
    samples[loopOTU,samples[loopOTU,]<loopMax] <- 0
    print(paste("Cleaning contaminants",contamOTU))
  }
  
  
  #now we add some taxonomy
  ##for the main samples
  assigned <- data.frame(matrix(nrow= dim(samples)[1],ncol=10))
  rownames(assigned) <- rownames(samples)
  colnames(assigned) <- colnames(Assignments)
  temp <- as.data.frame(apply(Assignments,2, as.character),stringsAsFactors=FALSE)
  for (row in 1:dim(samples)[1]){
    if(!rownames(samples)[row] %in% Assignments$OTU){
      next()}
    assigned[rownames(assigned)[row],] <-  temp[na.omit(match(rownames(assigned)[row],as.character(temp$OTU))),]
  }
  
  assigned <- cbind(samples,assigned)
  
  ##for the check samples
  
  assigned.c <- data.frame(matrix(nrow= dim(check)[1],ncol=10))
  rownames(assigned.c) <- rownames(check)
  colnames(assigned.c) <- colnames(Assignments)
  temp <- as.data.frame(apply(Assignments,2, as.character),stringsAsFactors=FALSE)
  for (row in 1:dim(check)[1]){
    if(!rownames(check)[row] %in% Assignments$OTU){
      next()}
    assigned.c[rownames(assigned.c)[row],] <-  temp[na.omit(match(rownames(assigned.c)[row],as.character(temp$OTU))),]
  }
  
  assigned.c <- cbind(check,assigned.c)
  
  #Now we write out the cleaned data

  write.csv(assigned,file="../cleaned/cleaned.dada2.data.csv")
  write.csv(assigned.c,file="../cleaned/checkdata.csv")

  
####====0.4 Cross Contam Test  ====####
#1.Identify OTUs unique to the 4 tanks from the natural site at time 0 compared to the artificial sites
contamDat <- assigned[c("S1.T0","S2.T0","S3.T0","S4.T0","S5.T0","S6.T0","S7.T0","S8.T0")]
#contamDat <- assigned[c("C1","C2","C3","C4","C5","C6","C7","C8")]  
contamDat[contamDat>0] <- 1
contamDat <- cbind(rowSums(contamDat[,1:4]),rowSums(contamDat[,5:8]))   
contamDat <- as.data.frame(contamDat)

uniqueOTus <- rownames(contamDat)[ifelse(contamDat$V1==0 & contamDat$V2==4,TRUE,FALSE)]

#2.Create a dataset containing only these OTUs
contamOTUdat <- assigned[match(uniqueOTus,rownames(assigned)),]

contamOTUdatexp <- contamOTUdat[,grep("S[1-4].T[1-9]",colnames(contamOTUdat))]
  
#3.Count the number of reads for OTUs assigned at low or high confidence in the experimental samples

highOTU <- na.omit(assigned$OTU[assigned$assignmentQual=="High"])
lowOTU <- na.omit(assigned$OTU[assigned$assignmentQual=="Low"])


highcontam <- contamOTUdatexp[na.omit(match(highOTU,rownames(contamOTUdatexp))),]
lowcontam <- contamOTUdatexp[na.omit(match(lowOTU,rownames(contamOTUdatexp))),]

#4.Count the total number of reads in the experimental samples

expSamples <- assigned[,grep("S[1-4].T[1-9]",colnames(assigned))]

readsperSample <- colSums(expSamples)
highreadsperSample <- colSums(highcontam)
lowreadsperSample <- colSums(lowcontam)

#5.Express cross contamination as a percentage

mean(highreadsperSample / readsperSample *100)

#6.Identity of high conf cross-contam OTUs 

taxaContam <- Assignments[match(rownames(highcontam),Assignments$OTU),]

test <- cbind(highcontam,Assignments[match(rownames(highcontam),Assignments$OTU),])

#descriptive statistics for reporting
highreadsperSample

mean(highreadsperSample)
sd(highreadsperSample)

mean(highreadsperSample / readsperSample *100)
sd(highreadsperSample / readsperSample *100)

#last four time points showing no reads at all indicating that very low levels of contam occured at the outset of the experiment

  
####====0.5 Rarefaction and Taxonomic Truncation  ====####

#Rarefy experimental dataset 

RexpSamples<-
  expSamples%>%
  filter(rowSums(.)>0)%>%
  t(.)%>%
  rrarefy(.,min(rowSums(.)))%>%
  t(.)
RexpSamples <- as.data.frame(RexpSamples)

#We will also create a proportional dataset incase the rarefaction causes problems
PexpSamples <- prop.table(as.matrix(expSamples),margin=2)
PexpSamples <- as.data.frame(PexpSamples[rowSums(PexpSamples)!=0,])

#1.Subset data with high confidence assignments
assignedHigh <- assigned[match(rownames(RexpSamples),assigned$OTU),]
assignedHigh <- assignedHigh[!is.na(assignedHigh$assignmentQual), ]
assignedHigh <- assignedHigh[assignedHigh$assignmentQual=="High", ]


#2.Combine ASVs with identical species level assignment

for (name in as.character(names(table(assigned[rownames(assignedHigh),"species"])[table(assigned[rownames(assignedHigh),"species"])>1]))){
  collapseOTUs <- as.character(na.omit(rownames(assignedHigh)[assignedHigh$species==name])) 
  MotherOTU <- names(sort(rowSums(RexpSamples[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  RexpSamples[MotherOTU,] <- RexpSamples[MotherOTU,] + colSums(RexpSamples[collapseOTUs,])
  RexpSamples <- RexpSamples[-match(collapseOTUs,rownames(RexpSamples)),]
}

rareTax <- data.frame(matrix(ncol=10,nrow = dim(RexpSamples)[1]))
colnames(rareTax) <- colnames(Assignments)
rareTax[is.na(rareTax)] <- "None"
rareTax$OTU <-rownames(RexpSamples) 
for (row in 1:length(rareTax$OTU)){
  rareTax[row,2:10] <- Assignments[match(rareTax$OTU[row],Assignments$OTU),2:10]
}

#3. Write out
write.csv(RexpSamples,"../cleaned/Rarefied.clean.csv")
write.csv(rareTax,"../taxonomy/Rarefied.clean.taxonomy.csv")



####====0.6 Descriptive Read Count  ====####

#1.Raw read count & SD
rawreadcount <- read.csv("ReadCounts.csv")

#total reads
sum(rawreadcount$reads.in)

mean(rawreadcount$reads.in[grep("[1-4].T[1-9].*",rawreadcount[,1])])
sd(rawreadcount$reads.in[grep("[1-4].T[1-9].*",rawreadcount[,1])])

#2.Read count per sample after QC

mean(colSums(assigned[,grep("[S1-4].T[1-9]",colnames(assigned))]))
sd(colSums(assigned[,grep("[S1-4].T[1-9]",colnames(assigned))]))

negControls <- c("C12","C14","C15","C16")
sum(controls[,na.omit(match(negControls,colnames(controls)))])
#only one negative control has any reads
sum(controls[,na.omit(match(negControls,colnames(controls)))])/4
sd(c(sum(controls[,na.omit(match(negControls,colnames(controls)))]),0,0,0))

#3.ASV count after QC
dim(RexpSamples)[1]


#4.Taxonomic assignment statistics 

table(rareTax$assignmentQual)

#high
sum(rowSums(RexpSamples)[grep("High",rareTax$assignmentQual)])/sum(rowSums(RexpSamples))
#low
sum(rowSums(RexpSamples)[grep("Low",rareTax$assignmentQual)])/sum(rowSums(RexpSamples))
#medium
sum(rowSums(RexpSamples)[grep("None",rareTax$assignmentQual)])/sum(rowSums(RexpSamples)) + sum(rowSums(RexpSamples)[is.na(rareTax$assignmentQual)])/sum(rowSums(RexpSamples))


####====0.7 Spiked Species Dataset  ====####

#1. How many spiked species do we see represented in the data?

spikedSpp <- read.csv("SpikedMetazoans.csv")

#species in the metabarcoding data
test <- rareTax$species[-5128]
spikedSpp$Binomial[spikedSpp$Binomial %in% test]

spikedSpp$Binomial[spikedSpp$Binomial %in% rareTax$species]
#species not found
notfound <- spikedSpp[!spikedSpp$Binomial %in% rareTax$species,] 

#lets look at genus 
notfound$Genus[notfound$Genus %in% rareTax$genus]

#lets get an OTU the represents each spp
spikedSpp$OTU <- rareTax$OTU[match(spikedSpp$Binomial,rareTax$species)]

#Redo the brittlestart by hand there is another OTU that has a low match

spikedSpp$OTU[spikedSpp$Binomial=="Amphipholis squamata"] <-"OTU_951"

#now lets set a representative OTU for each genus 
for (genus in  notfound$Genus[notfound$Genus %in% rareTax$genus]){
  runningOTUs <-na.omit(rareTax$OTU[rareTax$genus==genus])
  spikedSpp$OTU[spikedSpp$Genus==genus] <- names(sort(rowSums(RexpSamples)[runningOTUs],decreasing = TRUE))[1]
}

#Now we can output a fasta containing all the OTUs to check 

OTUs <- read.fasta("../taxonomy/dada2OTUs.fasta")
OTUseqs <- OTUs[na.omit(spikedSpp$OTU)]
write.fasta(OTUseqs,spikedSpp$Binomial[match(na.omit(spikedSpp$OTU),spikedSpp$OTU)],"../taxonomy/SpikedSppCheck.fa")


#2 Prep assigned only species dataset 

#collapse any genus level assignments 
genusout <-as.data.frame(matrix(nrow=length(notfound$Genus[notfound$Genus %in% rareTax$genus]),ncol=dim(RexpSamples)[2]))
colnames(genusout) <- colnames(RexpSamples)
rownames(genusout) <- notfound$Genus[notfound$Genus %in% rareTax$genus]

for (genus in notfound$Genus[notfound$Genus %in% rareTax$genus]){
print(genus)  
  loopOTUs <- na.omit(rareTax$OTU[rareTax$genus==genus])
  genusout[genus,] <- colSums(RexpSamples[loopOTUs,])
}

#subset spp assignments and combine with genus

spikedSpp$Binomial[spikedSpp$Binomial %in% test]
spikedSpp$OTU[spikedSpp$Binomial %in% test]

sppout <- RexpSamples[spikedSpp$OTU[spikedSpp$Binomial %in% test],]

rownames(sppout) <- spikedSpp$Binomial[spikedSpp$Binomial %in% test]

#write out data

write.csv(rbind(sppout,genusout),"../cleaned/RarefySppObservations.csv")

