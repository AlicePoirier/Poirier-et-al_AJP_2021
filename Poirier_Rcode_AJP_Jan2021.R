#
### R code for:  Poirier et al. 2021, American Journal of Primatology
### On the trail of primate scent signals: a field analysis of callitrichid 
### scent-gland secretions by portable gas chromatography-mass spectrometry
#


#### 1) Datasets  ----
# Dataset 1: Presence/absence of compounds, removed non-breeders
dataAll<-read.csv(file.choose(), header=T, na.strings="NA", sep=",")
dataAll[,-c(8,10)]<- lapply(dataAll[,-c(8,10)], factor)
head(dataAll)

# Dataset 2: Only samples containing compounds of interest
dataComp<-read.csv(file.choose(), header=T, na.strings = "NA", sep=",")
dataComp[]<- lapply(dataComp, factor)
head(dataComp)




### 2) Effect of storage time on the likelihood of detecting compounds ----
head(dataAll)
mean(dataAll$Storage)
median(dataAll$Storage)
sd(dataAll$Storage)

# Generalized linear model
Storage.glm<-glm(Presence ~ Storage, data=dataAll, family=binomial(link = "logit")) 
summary(Storage.glm) #ns

Storage.glm.R2<- 1-(Storage.lm$deviance/Storage.lm$null.deviance)
Storage.glm.R2

# Check glm residuals
library(DHARMa)
simulationOutput<- simulateResiduals(fittedModel=Storage.glm, n=1000, refit=F, use.u=F)
plot(simulationOutput) # all OK




### 3) Differences in chemical composition between categories of samples ----
### at the levels of species, group, sex, breeding status and sample type

### 3.1. Nonmetric multidimensional scaling (NMDS) ordination  ----
library(dplyr)
library(vegan)
head(dataComp)

A.class<-unique(dataComp[,c(1:7)])
A.class2<-arrange(A.class, Sample)
A.matrix<- xtabs(~Sample+Compound, dataComp)
A.vegdist<-vegdist(A.matrix, method="bray")
A.nmds<- metaMDS(A.vegdist, k=2, autotransform=F, noshare=F, wascores=F, 
                expand=T, trymax=200, trace=F)
A.nmds.scores<- scores(A.nmds) 
A.nmds$stress

A.NMDS.data<- cbind(A.nmds.scores, A.class2)
head(A.NMDS.data)




### 3.2. Analysis of similarity (ANOSIM) tests ----
# Species
anosim(A.vegdist, A.NMDS.data$Species, permutations=999, 
       distance="bray", strata=NULL)  #0.001

# Group within species
anosim(A.vegdist, A.NMDS.data$Group, permutations=999, 
       distance="bray", strata=A.NMDS.data$Species)  #0.001

# Sex within species
anosim(A.vegdist, A.NMDS.data$Sex, permutations=999, 
       distance="bray", strata=A.NMDS.data$Species)  #0.008

# Breding status within species
anosim(A.vegdist, A.NMDS.data$Repro, permutations=999, 
       distance="bray", strata=A.NMDS.data$Species)  #0.053

# Sample type within species
anosim(A.vegdist, A.NMDS.data$Type, permutations=999, 
       distance="bray", strata=A.NMDS.data$Species)  #0.146




### 3.3. Figure 2 ----
library(ggplot2)
library(cowplot)

# Species
fig2a<-
ggplot(A.NMDS.data, mapping=aes(x=NMDS1, y=NMDS2, 
                            colour=Species, linetype=Species))+
  geom_point(size=2, shape=16)+
  scale_colour_manual(values=c("#D55E00", "#56B4E9"),
                      labels=c("Saddleback","Emperor"))+
  xlab("NMDS 1") + ylab("NMDS 2")+ 
  labs(title=expression(paste(bold("a."))))+
  guides(linetype=F)+
  theme_bw(base_size=10)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))+
  theme(legend.position="right")

# Sex
fig2b<-
ggplot(A.NMDS.data, mapping=aes(x=NMDS1, y=NMDS2, shape=Sex, 
                            colour=Species, linetype=Species))+
  geom_point(size=2)+
  scale_colour_manual(values=c("#D55E00", "#56B4E9"))+
  scale_shape_manual(values=c(16,1), labels=c("Female","Male"))+
  xlab("NMDS 1") + ylab("NMDS 2")+ 
  labs(title=expression(paste(bold("b."))))+
  guides(colour=F, linetype=F)+
  theme_bw(base_size=10)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))+
  theme(legend.position="right")

# Breeding status
fig2c<-
ggplot(A.NMDS.data, mapping=aes(x=NMDS1, y=NMDS2, shape=Repro, 
                            colour=Species, linetype=Species))+
  geom_point(size=2)+
  scale_colour_manual(values=c("#D55E00", "#56B4E9"))+
  scale_shape_manual(values=c(16,1), name="Breeding\nstatus",
                     labels=c("Primary","Secondary"))+
  xlab("NMDS 1") + ylab("NMDS 2")+ 
  labs(title=expression(paste(bold("c."))))+
  guides(colour=F, linetype=F)+
  theme_bw(base_size=10)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))+
  theme(legend.position="right")

# Sample type
fig2d<-
ggplot(A.NMDS.data, mapping=aes(x=NMDS1, y=NMDS2, shape=Type, 
                            colour=Species, linetype=Species))+
  geom_point(size=2)+
  scale_colour_manual(values=c("#D55E00", "#56B4E9"))+
  scale_shape_manual(values=c(15,16,17,8), name="Sample\ntype")+
  xlab("NMDS 1") + ylab("NMDS 2")+ 
  labs(title=expression(paste(bold("d."))))+
  guides(colour=F, linetype=F)+
  theme_bw(base_size=10)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))+
  theme(legend.position="right")


# Plot and save figure
dev.new()
Figure2<- plot_grid(fig2a, fig2b, fig2c, fig2d, rel_heights=c(1,1,1,1), ncol=2, align="v")

ggsave(plot=Figure2, filename="AJPfig2_Jan21_600dpi.tiff", 
       device="tiff", scale=1.1, dpi=600, units="cm", width=15.9, height=10)


