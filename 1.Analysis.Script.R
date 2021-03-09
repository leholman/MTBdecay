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
library(nlstools)

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


#ASVwmod <- nls(ASVrich ~ times )
#selfStart(ASVrich ~ times)

#Time explains nMDS differences
envfit(MDSCOI,times)



####====2.0 Species Analysis====####
palette(rev(brewer.pal(12, "Set3")))

#Spp detection
#first by number of reads
#muck the data about
trSpp <- as.data.frame(t(rSpp))
trSpp <- gather(trSpp,key="spp",value="reads")
trSpp$times <- rep(times,dim(rSpp)[1])
trSpp$tank <-rep(c(rep(1,9),rep(2,9),rep(3,9),rep(4,9)),dim(rSpp)[1]) 

spp <- trSpp$spp
time <- trSpp$times
tank <- trSpp$tank

times <- trSpp$times[trSpp$spp==spp]
Rmax <- max(trSpp$reads)

#GGplot of raw data

pdf("figures/trial/SppReadsSep.pdf",width=7,height=5)
par(mar=c(4.1, 4.1, 1.5, 2.0), mfrow=c(1,1))
ggplot(trSpp, aes(times,reads)) + geom_point() + facet_wrap(vars(spp),scales="free_y")
#+ geom_smooth(method="lm", formula= (log(y+0.95) ~ x))
dev.off()

#Create a dataframe to collect data about models
modelData <- data.frame("Species"=unique(spp),
                        "Reads"=unique(spp),
                        "logL.R2"=rep(NA,length(unique(spp))),
                        "logL.p"=rep(NA,length(unique(spp))),
                        "logL.AIC"=rep(NA,length(unique(spp))),
                        "logL.AIC.2"=rep(NA,length(unique(spp))),
                        "logL.k"=rep(NA,length(unique(spp))),
                        "logL.95upr"=rep(NA,length(unique(spp))),
                        "logL.95lwr"=rep(NA,length(unique(spp))))


#Plot & model output showing poor fit of natural log linear model
# See https://stats.stackexchange.com/questions/48714/prerequisites-for-aic-model-comparison
# to correct AIC when log transforming data

sink("rawdata/loglinearmodelout.txt")
for (species in unique(spp)){
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  model <- lm(log(sample.reads+0.5)~sample.time)
  print(species)
  print(summary(model))
  modelData[modelData$Species==species,2] <- sum(sample.reads)
  modelData[modelData$Species==species,3] <- summary(model)$adj.r.squared
  modelData[modelData$Species==species,4] <- summary(model)$coefficients[2,4]
  modelData[modelData$Species==species,5] <- AIC(model)
  modelData[modelData$Species==species,6] <- AIC(model)+2*sum(log(sample.reads+0.5))
  modelData[modelData$Species==species,7] <- model$coefficients[2]
  modelData[modelData$Species==species,8] <- confint(model, level=0.95)[2,1]
  modelData[modelData$Species==species,9] <- confint(model, level=0.95)[2,2]
  
}
sink()

modelData2 <- modelData[modelData$logL.p<0.05,]
modelData2 <- modelData2[order(as.numeric(modelData2$Reads),decreasing = TRUE),]

pdf("figures/Figure2/Spp.pdf",width=7,height=3.5)
par(mfrow=c(2,4))
par(mar=c(2.0, 2.0, 1.0, 1.0))
colours <- 1

for (species in modelData2$Species){
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  model <- lm(log(sample.reads+0.5)~sample.time)
  predictions <- (exp(predict(model,data.frame("sample.time"=0:200),interval = "confidence")))-0.5
  plot(sample.time,sample.reads,type="n")
  polygon(x = c(0:200, 200:0),
          y = c(predictions[,3],rev(predictions[,2])),
          col =  adjustcolor(colours, alpha.f = 0.30), border = NA)
  lines(0:200,predictions[,1],lwd=3,col=colours)
  points(sample.time,sample.reads,pch=16,cex=1)
  legend("topright",paste("reads=",modelData2$Reads[modelData2$Species==species]),bty="n",cex=0.9)
  colours <- colours + 1
}

dev.off()

##Attempt @ tank effects following reviewers comments 

pdf("figures/Figure2/PerTankEffects.pdf",width=8,height=5)

par(mfrow=c(2,4))
par(mar=c(2.0, 2.0, 1.0, 1.0))
colours <- 1

for (species in modelData2$Species){
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  sample.tank <- tank[spp==species] 
  model <- lm(log(sample.reads+0.5)~sample.time)
  predictions <- (exp(predict(model,data.frame("sample.time"=0:200),interval = "confidence")))-0.5
  plot(sample.time,sample.reads,type="n")
  polygon(x = c(0:200, 200:0),
          y = c(predictions[,3],rev(predictions[,2])),
          col =  adjustcolor(colours, alpha.f = 0.30), border = NA)
  lines(0:200,predictions[,1],lwd=3,col=colours)
  #points(sample.time,sample.reads,pch=sample.tank,cex=1)
  text(sample.time,sample.reads,labels=sample.tank,cex=1)
  legend("topright",paste("reads=",modelData2$Reads[modelData2$Species==species]),bty="n",cex=0.9)
  colours <- colours + 1
}

dev.off()

for (species in unique(spp)){
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  sample.tank <- as.factor(tank[spp==species]) 
  model <- lm(log(sample.reads+0.5)~sample.time+sample.tank)
  print(species)
  print(summary(model))
}



#Non linear exponential decay model

#First we set up the model and save some data
sink("rawdata/non-linearmodelout2(nls).txt")

#We are going to grow the list dynamically, this is bad luck.
speciesReads <- c()
speciesK <- c()
speciesKerrUpr <- c()
speciesKerrLwr <- c()
speciesR2 <- c()
speciesAIC <- c()

for (species in unique(rownames(rSpp))){
  print(species)
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  model <- nls(sample.reads ~ SSasymp(sample.time, yf, y0, log_alpha))
  print(summary(model))
  #Extract read count
  runningReads <- sum(sample.reads)
  names(runningReads) <- species
  speciesReads <- c(speciesReads,runningReads)
  #Extract K
  runningK <- exp(coef(model)[3])
  names(runningK) <- species
  speciesK <- c(speciesK,runningK)
  #Extract K 95% cofInt
  runningKerrUpr <- round(exp(confint2(model)[3,1]),7)
  runningKerrLwr <- round(exp(confint2(model)[3,2]),7)
  names(runningKerrUpr) <- species
  names(runningKerrLwr) <- species
  speciesKerrUpr <- c(speciesKerrUpr,runningKerrUpr)
  speciesKerrLwr <- c(speciesKerrLwr,runningKerrLwr)
  #Extract R2
  RSS.p <- sum(residuals(model)^2) # Residual sum of squares
  TSS <- sum((sample.reads  - mean(sample.reads ))^2) # Total sum of squares
  runningR <- (1 - (RSS.p/TSS))  # R-squared measure
  names(runningR) <- species
  speciesR2 <- c(speciesR2,runningR)
  runningAIC <- AIC(model)
  names(runningAIC) <- species
  speciesAIC <- c(speciesAIC,runningAIC)
  
}
sink()

modelData$reads <- speciesReads
modelData$nls.k <- speciesK
modelData$nls.k.95upr <- speciesKerrLwr
modelData$nls.k.95lwr <- speciesKerrUpr 
modelData$nlsR2 <- speciesR2
modelData$nls.AIC <- speciesAIC
modelData <- modelData[order(modelData$reads,decreasing = T),] 

modelDataOut <- modelData[,1:2]
modelDataOut$Species <- metadatSpp$FinalID[match(modelDataOut$Species,metadatSpp$altID)]
modelDataOut <- cbind(modelDataOut,modelData[,c(3,4,7,8,9,6,14,11,12,13,15)])

write.csv(modelDataOut,"model.outputs.csv")


##PLot of data
pdf("figures/Figure2/Decay.pdf",width = 2,height=5)
par(mfrow=c(1,1),mar=c(2.0, 1.0, 1.0, 2.0))

plot(1:8,sqrt(modelData2$logL.k^2),
     ylim=c(sqrt(min(modelData2$logL.95lwr^2)),sqrt(max(modelData2$logL.95upr^2))),
     xlim=c(0.5,8.5),
     col=1:8,pch=3,cex=1,xaxt='n',yaxt="n")
axis(4)
mtext(text="Decay Model k Estimate (hr-1)",side=4,line=3)
for (sppN in 1:8){
  lines(c(sppN,sppN),c(sqrt(modelData2$logL.95upr[sppN]^2),sqrt(modelData2$logL.95lwr[sppN]^2)),col=sppN,lwd=3)
  }
dev.off()

#Now lets plot the relationship

pdf("figures/Figure2/Figure2.supp.pdf",width=7,height=5)
par(mfrow=c(3,4))
par(mar=c(2.0, 2.0, 1.0, 1.0))
colours <- 1
for (species in unique(spp)[order(speciesReads,decreasing = T)]){
  print(species)
  sample.text <- metadatSpp$FinalID[match(species,metadatSpp$altID)]
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  plot(sample.time,sample.reads,pch=16)
  model <- nls(sample.reads ~ SSasymp(sample.time, yf, y0, log_alpha))
  predictions <- predict(model,data.frame("sample.time"=0:200),interval = "confidence")
  lines(0:200,predictions,lwd=3,col=colours)
  points(sample.time,sample.reads,pch=16,cex=1)
  legend("topright",legend=c(sample.text,paste0(" reads=",sum(sample.reads))),text.font = c(3,1),bty="n",cex=0.9)
  colours <- colours +1
}
dev.off()



#What about looking a number of tanks each species was detected in?

#now by number of tanks with positive detection
trSppBinary <- trSpp
trSppBinary$reads[trSppBinary$reads>0] <- 1
trSppBinary.summed <- aggregate(trSppBinary$reads~trSppBinary$times,FUN=sum)
trSppBinary <- aggregate(trSppBinary$reads~trSppBinary$times+trSppBinary$spp,FUN=sum)
colnames(trSppBinary) <- c("times","spp","Nd")

pdf("figures/SppNd.pdf",width=7,height=5)
par(mfrow=c(3,4))
par(mar=c(2.0, 2.0, 1.0, 1.0))
colours <- 1
for (species in unique(spp)[order(speciesReads,decreasing = T)]){
  print(species)
  sample.text <- metadatSpp$FinalID[match(species,metadatSpp$altID)]
  sample.Nd <- trSppBinary$Nd[trSppBinary$spp==species]
  sample.time <- trSppBinary$times[trSppBinary$spp==species]
  plot(sample.time,sample.Nd,pch=16,ylim=c(0,4.2),col=colours,cex=1.2)
  legend("topright",legend=sample.text,text.font = c(3),bty="n",cex=0.9)
  colours <- colours +1
}
dev.off()


pdf("figures/SppNdsummed.pdf",width=5,height=3.5)
par(mfrow=c(1,1))
par(mar=c(4.0, 4.0, 1.0, 1.0))
plot(trSppBinary.summed$`trSppBinary$times`,trSppBinary.summed$`trSppBinary$reads`/4,pch=16,
     xlab="Time (hrs)",
     ylab="Average N detections")
dev.off()




####====Extra Questions====####

#1.does biomass correlate with reads?

#What role does approx weight have with reads?

reads <- rowSums(rSpp[na.omit(match(metadatSpp$altID,rownames(rSpp))),1:4])
weight <- metadatSpp$WeightApprox[metadatSpp$altID!=""]

plot(log10(reads),weight)

#non normal data
shapiro.test(reads)
shapiro.test(weight)
#normal is log10 transformed
shapiro.test(log10(reads))
shapiro.test(log10(weight))

#Kendalls Tau
cor.test(reads, weight, method=c("kendall"))

#linear model
summary(lm(log10(reads)~weight))

summary(lm(log10(reads)~log10(weight)))

#Do tanks have different decay patterns?




####====Code Basement====####

#Let's see if we can predict model fit based on read number?

modelfit.binary <- ifelse(modelfit$model.fit=="yes",1,0)
glm.1 <- glm(modelfit.binary~modelfit$reads,binomial)

summary(glm.1)
plot(glm.1)

#A simulation to test the fit of natural log when modeling exponential decay

sim.time <- 1:200
k <- 0.02
k.range <- seq(0.01,0.10,by=0.01)

sim.Reads <- (1000*(1-k)^sim.time) + 200
sim.ReadsNoise <- sim.Reads+rnorm(200,0,50)
sim.ReadsNoise[sim.ReadsNoise<0] <- 1


plot(sim.time,sim.ReadsNoise,pch=16,cex=0.8)

model <- lm(log(sim.ReadsNoise)~sim.time)
summary(model)
predictions <- (exp(predict(model,data.frame("time"=1:200),interval = "confidence")))
lines(1:200,predictions[,1],lwd=3,col="blue")


model2 <- nls(sim.ReadsNoise ~ SSasymp(sim.time, yf, y0, log_alpha))
summary(model2)
predictions <- predict(model2,data.frame("sim.time"=0:200),interval = "confidence")
lines(0:200,predictions,lwd=3,col="green")



(RSS.p <- sum(residuals(model2)^2)) # Residual sum of squares
(TSS <- sum((sim.ReadsNoise - mean(sim.ReadsNoise))^2))  # Total sum of squares
1 - (RSS.p/TSS)  # R-squared measure

#We simulate a bunch of ecologically typical k values
par(mfrow=c(2,5))
par(mar=c(2.0, 2.0, 1.0, 1.0))
k.range= seq(0.002,0.02,0.002)

for (k in k.range){
  #simulate decay with new k
  sim.time <- 1:200
  sim.Reads <- (1000*(1-k)^sim.time) + 200
  sim.ReadsNoise <- sim.Reads+rnorm(200,0,50)
  sim.ReadsNoise[sim.ReadsNoise<0] <- 1
  ##SubSample to 50 points
  rSample <- sample(1:200,36)
  sim.ReadsNoise <- sim.ReadsNoise[rSample]
  sim.time <- sim.time[rSample]
  
  #plot decay
  plot(sim.time,sim.ReadsNoise,pch=16,cex=0.8,ylim=c(0,1200))
  #model using log linear
  model <- lm(log(sim.ReadsNoise)~sim.time)
  predictions <- (exp(predict(model,data.frame("sim.time"=1:200),interval = "confidence")))
  lines(1:200,predictions[,1],lwd=2,col="blue",lty=2)
  #model using non linear
  model2 <- nls(sim.ReadsNoise ~ SSasymp(sim.time, yf, y0, log_alpha))
  predictions <- predict(model2,data.frame("sim.time"=0:200),interval = "confidence")
  lines(0:200,predictions,lwd=2,col="red",lty=2)
  text(150,1150,paste("k=",k))
  text(150,1050,paste("k.log=",round(sqrt(coef(model)[2]^2),3))) # not sure how to estimate a decay constant from this model, perhaps the slope?
  text(150,950,paste("k.nls=",round(exp(coef(model2)[3]),3))) 
}





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


##OLD non linear MODEL
sink("rawdata/non-linearmodelout(drm).txt")
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

pdf("figures/trial/SppReadsModelNonlinear.pdf",width=7,height=5)
par(mfrow=c(3,4))
par(mar=c(2.0, 2.0, 1.0, 1.0))
colours <- 1
for (species in unique(spp)[order(speciesReads,decreasing = T)]){
  print(species)
  sample.reads <- trSpp$reads[spp==species]
  sample.time <- time[spp==species]
  plot(sample.time,sample.reads,pch=16)
  model <- drm(sample.reads ~ sample.time , fct =EXD.2())
  predictions <- predict(model,data.frame("sample.time"=0:200),interval = "confidence")
  polygon(x = c(0:200, 200:0),
          y = c(predictions[,3],rev(predictions[,2])),
          col =  adjustcolor(colours, alpha.f = 0.30), border = NA)
  lines(0:200,predictions[,1],lwd=3,col=colours)
  points(sample.time,sample.reads,pch=16,cex=1)
  legend("topright",legend=paste0(species," reads=",sum(sample.reads)),bty="n",cex=0.9)
  colours <- colours +1
}
dev.off()
sink()

##
#one idea for the confidence intervals with non linear models
#this doesn't work at the moment as various packages need to be updaded for the lastets verison of R
require(propagate)

pred_model <- predictNLS(model, newdata=mm)
conf_model <- pred$summary
plot(v~S)

lines(conf_model$Prop.Mean.1 ~ S, lwd=2)
lines(conf_model$"Sim.2.5%" ~ S, lwd=1)
lines(conf_model$"Sim.97.5%" ~ S, lwd=1)




##Testing AICs

model1 <- lm(log(sample.reads+0.5)~sample.time)

model2 <- nls(sample.reads ~ SSasymp(sample.time, yf, y0, log_alpha))

AIC(model1,model2)

AIC(model1)+2*sum(log(sample.reads+0.5))
AIC(model2)


#calculation for lm
logL <- 0.5 * (- length(residuals(model1)) * (log(2 * pi) + 1 - log(length(residuals(model1))) + log(sum(residuals(model1)^2))))
logL.nls <- 0.5 * (- length(residuals(model2)) * (log(2 * pi) + 1 - log(length(residuals(model2))) + log(sum(residuals(model2)^2))))


#lm AIC (ordinary)
2 * (model1$rank + 1) - 2 * logL

#lm AIC (log transformation corrected)
(2 * (model1$rank + 1) - 2 * logL)+2*sum(log(sample.reads+0.5))

#nls AIC
2 * (4) - 2 * logL.nls


