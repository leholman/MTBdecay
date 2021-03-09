#############################################
####==== Metabarcoding Decay Analysis ===####
####====         16.07.2020          ====####
#############################################

###Script 1 - Analysis###

####====0.0 Packages====####

library(data.table)
library(metabarTOAD)
library(vegan)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(drc)
library(aomisc)

####====0.1 Data Import====####
rData <- read.csv("cleaned/Rarefied.clean.csv",row.names = 1)
rSpp <- read.csv("cleaned/RarefySppObservations.csv",row.names = 1)
rTaxa <- read.csv("taxonomy/Rarefied.clean.taxonomy.csv",row.names = 1)
metadat <- read.csv("metadata.csv")
metadatSpp <- read.csv("rawdata/SpikedMetazoans.csv")
times <- metadat$Hours[match(gsub("S","",colnames(rData)),metadat$ID)]

####====1.0 Entire Dataset Analysis====####

#OTUs over time
pdf("figures/Figure1/WholComm1.pdf",width=8,height=4)
par(mar=c(4.1, 4.1, 1.5, 1.5), mfrow=c(1,2))

#OTUs over time w/ nonlinear modeling


ASVrich <- colSums(rData != 0)
ASVwmod <- drm(ASVrich ~ times , fct =W1.4 ())
predictions <- predict(ASVwmod,data.frame("times"=0:200),interval = "confidence")
plot(times,ASVrich,type="n",xlab="Time (h)",ylab="ASV Richness")
polygon(x = c(0:200, 200:0),
        y = c(predictions[,3],rev(predictions[,2])),
        col =  adjustcolor("red", alpha.f = 0.10), border = NA)
lines(0:200,predictions[,1],lwd=3,col="red")
points(times,ASVrich,pch=16,col="black")

#Exponential Model
#ExpMod <- lm(log(ASVrich) ~ times)
#predictions <- exp(predict(ExpMod,data.frame("times"=0:200),interval = "confidence"))

#nls exponential 
#ExpMod.2 <- nls(ASVrich ~ SSasymp(times, yf, y0, log_alpha))
#predictions <- exp(predict(ExpMod.2,data.frame("times"=0:200),interval = "confidence"))


#nMDS plot
MDSCOI <- metaMDS(vegdist(t(rData),distance = "jaccard",binary = TRUE))
MDSCOIdat <-as.data.frame(MDSCOI$points)

plot(MDSCOIdat,type="n")
points(MDSCOIdat$MDS1,MDSCOIdat$MDS2,pch=16,col="darkgrey")
timefactor <- as.factor(gsub("S[1-4].","",rownames(MDSCOIdat)))
surfout<- ordisurf(MDSCOI,metadat$Hours[match(gsub("S","",colnames(rData)),metadat$ID)],method="REML",add = TRUE,col="red",lty=2,lwd=2)
text(unlist(aggregate(MDSCOIdat$MDS1~timefactor,FUN=mean)[2]),
     unlist(aggregate(MDSCOIdat$MDS2~timefactor,FUN=mean)[2]),
     labels=levels(timefactor))
dev.off()

#Statistics

#4 param weibull model defined as follows
#b - curve steepness
#c - ASV lower limit
#d - ASV initial richness
#e - hours for inflection 

summary(ASVwmod)
modelFit(ASVwmod)
#test a null by randomizing the order of samples to get an idea of a non-significant result 
ASVwmod.null <- drm(ASVrich[sample(1:36,36)] ~ times , fct =W1.4 ())
test <- summary(ASVwmod.null)


ASVwmod <- nls(ASVrich ~ times )
selfStart(ASVrich ~ times)

#Check time explains nMDS differences
envfit(MDSCOI,times)



####====2.0 Species Analysis====####

#Spp detection
#first by number of reads
#muck the data about
trSpp <- as.data.frame(t(rSpp))
trSpp <- gather(trSpp,key="spp",value="reads")
trSpp$times <- rep(times,dim(rSpp)[1])

Rmax <- max(trSpp$reads)

#plot the data
pdf("figures/trial/SppReads.pdf",width=7,height=5)
par(mar=c(4.1, 4.1, 1.5, 1.5), mfrow=c(1,1))
count <- 1
palette(rev(brewer.pal(12, "Set3")))
plot("n",ylim=c(-1,log(Rmax)),xlim =c(0,200),ylab="log N reads",xlab="Time (hours)")

for (spp in unique(trSpp$spp)){
  print(spp)
  
  #make a temperoy object containing a number below one to sum with values in order to avoid logging of 0s
  temp <- rep(0,length(trSpp$reads[trSpp$spp==spp]))
  temp[trSpp$reads[trSpp$spp==spp]==0] <- 0.95
  
  #use modelling to produce log devay
  times <- trSpp$times[trSpp$spp==spp]
  mod.log<-lm(log(trSpp$reads[trSpp$spp==spp]+temp)~times)
  
  # predict along predictor variable range
  newdat <- data.frame(times=seq(0,200,0.5))
  newdat <-cbind(newdat,predict(mod.log, newdat, type="response",interval = "confidence",level=0.95,se.fit=TRUE))
  
  # plot
  points(trSpp$times[trSpp$spp==spp],log(trSpp$reads[trSpp$spp==spp]+temp),pch=16,col=count)
  
  lines(newdat$times,newdat$fit.fit,col=count)
  
  count <- count+1
}
dev.off()

#GGplot

pdf("figures/trial/SppReadsSep.pdf",width=7,height=5)
par(mar=c(4.1, 4.1, 1.5, 2.0), mfrow=c(1,1))
ggplot(trSpp, aes(times,reads)) + geom_point() + facet_wrap(vars(spp),scales="free_y")
#+ geom_smooth(method="lm", formula= (log(y+0.95) ~ x))
dev.off()


#test the fit of natural log

#simulation

sim.time <- 1:200

sim.Reads <- (1000*(1-0.03)^sim.time) + 200


sim.ReadsNoise <- sim.Reads+rnorm(200,0,100)
sim.ReadsNoise[sim.ReadsNoise<0] <- 1

plot(sim.time,sim.Reads)
plot(sim.time,sim.ReadsNoise)

model <- lm(log(sim.ReadsNoise)~sim.time)
summary(model)
predictions <- (exp(predict(model,data.frame("time"=1:200),interval = "confidence")))
lines(1:200,predictions[,1],lwd=3,col="blue")
lines(sim.Reads,col="red",lwd=3)

model <- drm(sim.ReadsNoise~sim.time, fct = EXD.3())
summary(model)
predictions <- predict(model,data.frame("sim.time"=0:200),interval = "confidence")
lines(0:200,predictions[,1],lwd=3,col="red")


spp <- trSpp$spp
time <- trSpp$times

pdf("figures/trial/SppReadsModel.pdf",width=7,height=5)
par(mfrow=c(3,4))
par(mar=c(2.0, 2.0, 1.0, 1.0))
colours <- 1

for (species in unique(spp)){
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  model <- lm(log(sample.reads+0.5)~sample.time)
  print(species)
  print(summary(model))
  predictions <- (exp(predict(model,data.frame("sample.time"=0:200),interval = "confidence")))-0.5
  plot(sample.time,sample.reads,type="n")
  polygon(x = c(0:200, 200:0),
          y = c(predictions[,3],rev(predictions[,2])),
          col =  adjustcolor(colours, alpha.f = 0.30), border = NA)
  lines(0:200,predictions[,1],lwd=3,col=colours)
  points(sample.time,sample.reads,pch=16,cex=1)
  legend("topright",legend=species,bty="n",cex=0.9)
  colours <- colours + 1
  
}

dev.off()



#Now lets try a non-linear model since the fit of the above is poor
sink("rawdata/non-linearmodelout.txt")
write.csv(data.frame("spp"=unique(spp),"model.fit"=rep(NA,12)),file="rawdata/modelfit.csv")
speciesReads <- c()
for (species in unique(rownames(rSpp))){
  print(species)
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  model <- drm(sample.reads ~ sample.time , fct =EXD.2())
  print(summary(model))
  runningReads <- sum(sample.reads)
  names(runningReads) <- species
  speciesReads <- c(speciesReads,runningReads)
  print(modelFit(model))
}
sink()

modelfit <- read.csv("rawdata/modelfit2.csv",row.names = 1)
modelfit$reads <- speciesReads


pdf("figures/trial/SppReadsModelNonlinear.pdf",width=7,height=5)
par(mfrow=c(3,4))
par(mar=c(2.0, 2.0, 1.0, 1.0))
colours <- 1
for (species in unique(spp)[order(modelfit$reads,decreasing = T)]){
  print(species)
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  plot(sample.time,sample.reads,pch=16)
  
  #if(modelfit$model.fit[modelfit$spp==species]=="yes"){
  model <- drm(sample.reads ~ sample.time , fct =EXD.2())
  predictions <- predict(model,data.frame("sample.time"=0:200),interval = "confidence")
  polygon(x = c(0:200, 200:0),
          y = c(predictions[,3],rev(predictions[,2])),
          col =  adjustcolor(colours, alpha.f = 0.30), border = NA)
  lines(0:200,predictions[,1],lwd=3,col=colours)
  points(sample.time,sample.reads,pch=16,cex=1)
 # }
  legend("topright",legend=paste0(species," reads=",sum(sample.reads)),bty="n",cex=0.9)
  colours <- colours +1
}
dev.off()
sink()

#Let's see if we can predict model fit based on read number?

modelfit.binary <- ifelse(modelfit$model.fit=="yes",1,0)
glm.1 <- glm(modelfit.binary~modelfit$reads,binomial)

summary(glm.1)
plot(glm.1)



#try a weibull model? using W1.2() under fct in model description??


#now by number of tanks with positive detection
trSppBinary <- trSpp
trSppBinary$reads[trSppBinary$reads>0] <- 1
trSppBinary <- aggregate(trSppBinary$reads~trSppBinary$times+trSppBinary$spp,FUN=sum)
colnames(trSppBinary) <- c("times","spp","Nd")

#plot the data
pdf("figures/trial/SppNd.pdf",width=7,height=5)
par(mar=c(4.1, 4.1, 1.5, 2.0), mfrow=c(1,1))
count <- 1
palette(rev(brewer.pal(12, "Set3")))
plot("n",ylim=c(0,4),xlim =c(0,200),ylab="N tanks spp detected",xlab="time (hours)")

for (spp in unique(trSppBinary$spp)){
  print(spp)
  
  # plot
  lines(jitter(trSppBinary$times[trSppBinary$spp==spp],2),trSppBinary$Nd[trSppBinary$spp==spp],pch=16,col=count,type="b")
  
  count <- count+1
}
dev.off()


pdf("figures/trial/SppNdSep.pdf",width=7,height=5)
par(mar=c(4.1, 4.1, 1.5, 2.0), mfrow=c(1,1))
ggplot(trSppBinary, aes(times,Nd)) + geom_point()+ facet_wrap(vars(spp))
dev.off()

#What role does approx weight have with reads?

reads <- rowSums(rSpp)[na.omit(match(metadatSpp$altID,rownames(rSpp)))]
weight <- metadatSpp$WeightApprox[metadatSpp$altID!=""]
mismatchF <- metadatSpp$Fmismatch[metadatSpp$altID!=""]
mismatchR <- metadatSpp$Rmismatch[metadatSpp$altID!=""]


#not normally distributed
shapiro.test(reads)
shapiro.test(weight)

#Borderline significance with a one tailed kendalls tau
cor.test(reads, weight, method=c("kendall"),alternative="greater")
cor.test(reads, weight)
#Cant use these tests because poor data quality
summary(lm(reads~weight))
plot(lm(reads~weight))






