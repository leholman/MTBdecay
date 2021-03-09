###Metabarcoding decay analysis 
##03.01.2019


#Load in packages
library("dplyr")
library("vegan")
library("reshape")
library("reshape2")
library("ggplot2")
library("stringr")
library("tidyr")
library("RColorBrewer")
library("readxl")
library("maps")
library("mapdata")
require(ggplot2)
require(ggmap)
require(mapproj)
require(rgeos)
require(maptools)
require(sp)
require(raster)
require(rgdal)
require(dismo)
library("RgoogleMaps")
library("reshape")


#Set some variables 
minreads <- 3
items <- NULL

#Set the seed 
set.seed("123456")


#### Data Cleaning ####
setwd("raw.data/")
metadat<-read.csv("../locations.csv") 
metadat$RealID <- gsub("X","S",make.names(metadat$RealID))
file <- system2('ls',stdout=TRUE)


  rawdat <-read.csv(file=file)
  colnames(rawdat) <- gsub("X","S",colnames(rawdat))
  
  
  #Seperate controls and samples
  samples <- rawdat[colnames(rawdat) %in% metadat$RealID[metadat$Type!="control2"]]
  controls <-  rawdat[colnames(rawdat) %in% metadat$RealID[metadat$Type=="control2"]]
  
  
  #Filter 1 - minimum number of reads for any ID
  samples[samples< minreads] <- 0
  samples <- samples[rowSums(samples) > 0,]
  
  #Filter 2 - within samples OTU must appear in more than one sample (this works becuase there are lots of reps per site and sample)
  filtersam <- samples
  filtersam[filtersam>0 ] <- 1
  filtersam <-filtersam[rowSums(filtersam) > 1,]
  samples <- samples[rownames(samples) %in% rownames(filtersam),]
  
  #Filter 3 -Make the maximum umber of reads for each OTU in the contam the zero value in the main data
  controlsCONTAM <- controls[rowSums(controls) > 0,]
  for (contamOTU in 1:length(controlsCONTAM[,1])){
    loopOTU <- row.names(controlsCONTAM[contamOTU,])
    loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
    if (any(is.na(samples[loopOTU,]))){next}
    samples[loopOTU,samples[loopOTU,]<loopMax] <- 0
    print(paste("Cleaning contaminants",contamOTU))
  }
  
  

  #Now we write out the cleaned data
  
  newname <- paste("../cleaned/","Cleaned.",file,sep="")
  write.csv(samples,file=newname)
  norm <-samples%>%
    filter(rowSums(.)>0)%>%
    t(.)%>%
    rrarefy(.,min(rowSums(.)))%>%
    t(.)
  norm <- as.data.frame(norm)
  rownames(norm) <- rownames(samples)
  write.csv(norm,file="../cleaned/normalised.csv")
  
  #First lets muck around with the data to subset the taxonomy by the OTUs
  taxo <- read.csv(file="../taxonomoy/Taxonomy.csv",header=FALSE)
  newtax <- taxo[as.character(taxo$V1) %in% row.names(samples),]
  write.csv(newtax,file="../taxonomoy/newtax.csv")
  
  
  
  
  ##First lets look to see if that number of OTUs decreases as a function of time 
  cleaneddat <- read.csv("../cleaned/Cleaned.Leraylulu.unoise3.csv",row.names = 1)
  rCOI<-
    cleaneddat%>%
    filter(rowSums(.)>0)%>%
    t(.)%>%
    rrarefy(.,min(rowSums(.)))%>%
    t(.)
  cleanedOTUcount <- as.data.frame(apply(as.data.frame(rCOI), 2,function(x) length(x[x>0])))
  cleanedOTUcount$time <- substr(rownames(cleanedOTUcount),5,5)
  cleanedOTUcount$cond <- c(rep("1",40),rep("2",39),rep("NA",11))
  colnames(cleanedOTUcount)[1] <- "OTUcount"
  
  pdf(file="../plots/boxplots.pdf",width=10,height=5)
  par(mfrow=c(1,2))
  boxplot(cleanedOTUcount$OTUcount[cleanedOTUcount$cond=="1"]~cleanedOTUcount$time[cleanedOTUcount$cond=="1"],xlab="TimePoint",ylab="Number of OTUs",col="grey",main="Tanks 1-4")
  boxplot(cleanedOTUcount$OTUcount[cleanedOTUcount$cond=="2"]~cleanedOTUcount$time[cleanedOTUcount$cond=="2"],xlab="TimePoint",ylab="Number of OTUs",col="grey",main="Tanks 5-8")
  dev.off()
  
  ##alldata
  cleaneddat <- read.csv("../cleaned/Cleaned.Leraylulu.unoise3.csv",row.names = 1)
  rCOI<-
    cleaneddat%>%
    filter(rowSums(.)>0)%>%
    t(.)%>%
    rrarefy(.,min(rowSums(.)))%>%
    t(.)
  
  MDSCOI <- metaMDS(t(rCOI),distance = "bray")
  MDSCOIdat <-as.data.frame(MDSCOI$points)
  plot(MDSCOIdat,type="n")
  text(MDSCOIdat$MDS1,MDSCOIdat$MDS2,labels=colnames(cleaneddat))
  
  ##onlydecaysamples
  samplescleaned <- cleaneddat[colnames(cleaneddat) %in% metadat$RealID[metadat$Type=="sample"]]
  rCOI<-
    samplescleaned%>%
    filter(rowSums(.)>0)%>%
    t(.)%>%
    rrarefy(.,min(rowSums(.)))%>%
    t(.)
  MDSCOI <- metaMDS(t(rCOI),distance = "bray")
  MDSCOIdat <-as.data.frame(MDSCOI$points)
  pdf(file = "../plots/nMDS.pdf")
  par(mfrow=c(1,1))
  plot(MDSCOIdat,type="n")
  text(MDSCOIdat$MDS1,MDSCOIdat$MDS2,labels=colnames(samplescleaned))
  
  
  colfunc <- colorRampPalette(c("darkblue", "lightblue"))
  plot(1:10,col=colfunc(10))
  
  dev.off()
  
 

