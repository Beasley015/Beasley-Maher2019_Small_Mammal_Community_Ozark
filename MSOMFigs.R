########################################################
# Running additional analyses & creating figures from  #
# MSOM output.                                         #
#                                                      #
# Data can be found in the file 'summer_mamm.csv.      #
# Model outputs are 'modelallglades.rds' and           #
# 'modelsampledglades.rds'. Covariates saved as RDS    #
# files 'area1', 'shape1', 'dist1', and 'dist2'.       #
########################################################


#Load packages and set working directory -----------------------------------
library(vegan)
library(betapart)
library(patchwork)
library(plyr)
library(tidyverse)
library(glmulti)
library(gam)

wd <- "c:/users/beasley/dropbox/manuscripts/Mamm IBT"
setwd(wd)

#Set up mammal data -------------------------------------------------------------
#all mammal data
all.mamm<-read.csv("mammibtrawdata/summer_mamm.csv", stringsAsFactors = F)
mamm.data<-all.mamm[,c("Date","Year","Glade","Day","Species")]
mamm.data$Occ <- rep(1, dim(mamm.data)[1]) #Add a column to denote presences

#The next few lines are needed to include sites with 0 captures
#If each site had at least 1 capture event, this bit isn't necessary
mamm.data$Date<-as.character(mamm.data$Date)
mamm.data[617,]<-c("9/6/2017",2,23,1,"Peromyscus attwateri",0) 
mamm.data[618,]<-c("12/7/2017",2,28,1,"Peromyscus attwateri",0)

mamm.data$Glade <- as.numeric(mamm.data$Glade)
mamm.data$Year <- as.integer(mamm.data$Year)
mamm.data$Occ <- as.integer(mamm.data$Occ)

#Create blank vector of number of observed species
unique.species<-as.character(unique(mamm.data$Species))
species.abs <- c("REFU", "SIHI", "PEAT", "PELE", "PEMA", "TAST", "NEFL", "SYFL")

#Create vector of site numbers
unique.sites<-as.character(c(1:32))

#Create presence-absence matrices for each year
presabs <- function(year){
  #Convert mammal data to wide form
  mamm.data %>%
    filter(Year == year) %>%
    group_by(Species, Glade) %>%
    summarise(sum = sum(Occ))  %>%
    spread(key = Glade, value = sum, fill = 0) %>%
    {. ->> z}
  
  #change to matrix and covert all values > 1 to 1
  rows <- z$Species
  z <- as.matrix(z[,-1])
  z[which(z >= 1)] <- 1
  rownames(z) <- rows
  return(z)
}

spoc16 <- presabs(year = 1)
spoc17 <- presabs(year = 2)
spoc17 <- rbind(spoc17, rep(0, ncol(spoc17))) #Didn't catch t.striatus in yr 2
rownames(spoc17)[8] <- "Tamias striatus"

#Read in model and covariates ---------------------------------------------
#Sampled glades only
model <- readRDS(file = "mammibtmodeloutputs/modelsampledglades.rds")
#All glades
# model <- readRDS(file = "mammibtmodeloutputs/modelallglades.rds")

#environmental covariates
#Area
polys1<-readRDS(file = "mammibtmodeloutputs/Area1.rds")

#Isolation: sampled glades only
isos<-readRDS(file = "mammibtmodeloutputs/dist1.rds")
#Isolation: all glades
# isos <- readRDS(file = "mammibtmodeloutputs/dist2.rds")

#shape
shape1<-readRDS(file = "mammibtmodeloutputs/shape1.rds") 

#Summarize community-level covariate effects ---------------------------------

# Extract covariates
#Area
a1<-model$sims.list$a1
a1<-a1[,1:8] # This removes unobserved species

#Isolation
a2<-model$sims.list$a2
a2<-a2[,1:8]

#Shape
a3<-model$sims.list$a3
a3<-a3[,1:8]

#Year
a4<-model$sims.list$a4
a4<-a4[,1:8]

#Look at quantiles for significance
quantile(a1, c(0.025, 0.975))
quantile(a2, c(0.025, 0.975))
quantile(a3, c(0.025, 0.975))
quantile(a4, c(0.025, 0.975))

#Summarize species-level covariate effects
a1.means<-apply(a1,2,mean)
apply(a1,2,quantile,probs = c(0.025, 0.975))

a2.means<-apply(a2,2,mean)
apply(a2,2,quantile,probs = c(0.025, 0.975))

a3.means<-apply(a3,2,mean)
apply(a3,2,quantile,probs = c(0.025, 0.975))

a4.means<-apply(a4,2,mean)
apply(a4,2,quantile,probs = c(0.025, 0.975))

#detection cov
b1<-model$sims.list$b1
b1<-b1[,1:8]
quantile(b1, c(0.025, 0.975))

b1.means<-apply(b1,2,mean)
apply(b1,2,quantile,probs=c(0.025, 0.975))

#Set up area matrix- matrix needs to include mean and 95% confidence interval
species.area<-data.frame(Species = species.abs, Mean = a1.means, 
                         Lo = apply(a1, 2, quantile, 0.025), 
                         hi = apply(a1, 2, quantile, 0.975))

#And isolation
species.isolation<-data.frame(Species = species.abs, Mean = a2.means, 
                              Lo = apply(a2, 2, quantile, 0.025), 
                              hi = apply(a2, 2, quantile, 0.975))

#Shape
species.ratio<-data.frame(Species = species.abs, Mean = a3.means, 
                          Lo = apply(a3, 2, quantile, 0.025), 
                          hi = apply(a3, 2, quantile, 0.975))

#Year
species.year<-data.frame(Species = species.abs, Mean = a4.means, 
                         Lo = apply(a4, 2, quantile, 0.025), 
                         hi = apply(a4, 2, quantile, 0.975))

#detection
spec.date<-data.frame(Species = species.abs, Mean = b1.means, 
                      Lo = apply(b1, 2, quantile, 0.025), 
                      hi = apply(b1, 2, quantile, 0.975))

#PLOT THAT NONSENSE! ----------------------------------------------
#Species-area graph
a<-ggplot()+
  geom_point(aes(x=species.area$Species,y=species.area$Mean), size = 1.8)+
  geom_errorbar(aes(x=species.area$Species, ymin = species.area$Lo, 
                    ymax = species.area$hi), width = 0.3)+
  geom_hline(aes(yintercept=0), linetype = 2)+
  labs(y = "Coefficient", title = "A")+
  scale_y_continuous(limits = c(-6,5))+
  theme_bw()+
  theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12), 
        axis.text.x = element_blank(), title = element_text(size = 18),
        panel.grid = element_blank())

#Species-isolation graph
b<-ggplot(data = species.isolation)+
  geom_point(aes(x=species.isolation$Species,y=species.isolation$Mean), size = 1.8)+
  geom_errorbar(aes(x=species.isolation$Species, ymin = species.isolation$Lo, 
                    ymax = species.isolation$hi), width = 0.3)+
  geom_hline(aes(yintercept=0), linetype = 2)+
  labs(title = "B")+
  scale_y_continuous(limits = c(-6,5))+
  theme_bw()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12), title = element_text(size = 18),
        panel.grid = element_blank())

#Species-shape graph
c<-ggplot()+
  geom_point(aes(x=species.ratio$Species,y=species.ratio$Mean), size=1.8)+
  geom_point(aes(x=c("PEAT","PEMA","SIHI"), y=c(1.5, 1.5, 4)), size = 1.8, shape = 8)+
  geom_errorbar(aes(x=species.ratio$Species, ymin = species.ratio$Lo, 
                    ymax = species.ratio$hi), width = 0.3)+
  geom_hline(aes(yintercept=0), linetype = 2)+ 
  scale_y_continuous(limits = c(-6,5))+
  labs(y = "Coefficient", title = "C")+
  theme_bw()+
  theme(axis.title.x=element_blank(), axis.text = element_text(size = 12), 
        title = element_text(size = 18),panel.grid = element_blank())

#Species-year graph
d<-ggplot(data = species.year)+
  geom_point(aes(x=species.year$Species,y=species.year$Mean), size=1.8)+
  geom_point(aes(x=species.year$Species, y = 1), size = 1.8, shape = 8)+
  geom_errorbar(aes(x=species.year$Species, ymin = species.year$Lo, 
                    ymax = species.year$hi, width = 0.3))+
  geom_hline(aes(yintercept=0), linetype = 2)+ 
  scale_y_continuous(limits = c(-6,5))+
  labs(title = "D")+
  theme_bw()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text = element_text(size = 12), title = element_text(size = 18),
        panel.grid = element_blank())

#Put plots in same frame
allcoef<-(a|b)/(c|d)
allcoef

# ggsave(allcoef, filename = "allcoef.tiff", dpi = 600, width = 178, units = "mm",
#        scale = 1.5)

#date graph
ggplot(data = spec.date)+
  geom_point(aes(x=spec.date$Species,y=spec.date$Mean), size=3)+
  geom_errorbar(aes(x=spec.date$Species, ymin = spec.date$Lo, 
                    ymax = spec.date$hi), size=1.2)+
  geom_hline(aes(yintercept=0), linetype = 2, size=1.2)+
  theme_bw()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text.x = element_text(size = 22))

#Empirical and estimated richness -----------------------------------
#Get empirical richness counts
spec.estimate1<-estimateR(t(spoc16))
richness1<-spec.estimate1[1,]

spec.estimate2<-estimateR(t(spoc17))
richness2<-spec.estimate2[1,]

obs.richness<-c(richness1, richness2)

#See estimates of total richness (N) by site
N <- model$sims.list$N
mean(N);summary(N);plot(table(N))

Nsite<-model$sims.list$Nsite
site.richness <- cbind(apply(Nsite,2,mean), apply(Nsite,2,sd))
est.richness<-site.richness

#Create vector of site #s
sites<-c(1:16)

#Glm's: corrected richness vs area/iso ---------------------------
#Vectors of years
year1<-rep(1, 16)
year2<-rep(2, 16)

#Construct data frame for linear model
gldata<-data.frame(Richness = est.richness[,1], Area = polys1, Iso = isos, 
                   Shape = shape1, Year = c(year1, year2))

#Linear models
l1<-lm(data=gldata, formula = Richness~Area+Iso+Shape+Year)
l2<-lm(data = gldata, formula = Richness~Area+Year)

aiccs<-glmulti(y = l1, level = 1, crit = "aicc")

#Plot of corrected richness vs. area
ggplot(data = gldata, aes(x=Area, y=Richness))+
  geom_point(aes(color=factor(gldata$Year)), size=1.5)+
  geom_smooth(aes(color = factor(gldata$Year)), method = "lm", se = FALSE, size = 1)+
  scale_color_manual(breaks=gldata$Year, values = c("black","limegreen"), name="Year")+
  labs(x = expression('Log'[10]*'Area (m'^2*')'), y = "Species Richness")+
  theme_bw()+
  theme(axis.title = element_text(size = 18), legend.title = element_text(size=18), 
        axis.text = element_text(size = 12), legend.text = element_text(size=12),
        panel.grid = element_blank())

# ggsave(filename = "arearegress.tiff", dpi = 600, width = 85, scale = 1.5, height = 60,
#        units = "mm")

#Look at turnover --------------------------------------------------
#use Z (estimated) occurence matrices derived by the model
Zs<-aperm(model$sims.list$Z, c(2,3,1))
Zs<-alply(Zs, 3)
Zs<-as.list(Zs)

#separate into years 1 and 2
z1<-lapply(Zs, '[',1:16,)
z2<-lapply(Zs, '[',17:32,)

turnover<-list()

#Use function beta.temp to compare sampling years
for(i in 1:length(Zs)){
  turnover[[i]]<-beta.temp(x=z1[[i]], y=z2[[i]], index.family = "jaccard")
}

turn.mat<-lapply(turnover, as.matrix)

jacs<-lapply(turn.mat, function(x)x[,3]) 
#calculates jaccard index between years 1 and 2 for each pair of matrices in the arrays

jacsmat<-do.call(rbind, jacs) #put all jaccard estimates in a single matrix

jaccard<-apply(jacsmat, 2, mean) #Mean jaccard estimate per site 

#run some correlation tests between Jaccard index and covariates
cor.test(x = polys1[1:16], y = jaccard, method = "spearman")
cor.test(x = isos[1:16], y = jaccard, method = "spearman")
cor.test(x = shape1[1:16], y = jaccard, method = "spearman")

#Correlation between isolation and Jaccard is significant. Run a lm
jac.model <- lm(jaccard~isos[1:16]) #This is also significant

#Now run a GAM
jac.gam <- gam(jaccard~isos[1:16]) #Also significant

#Quadratic model
isos.square <- (isos[1:16]^2)
jac.quad <- lm(jaccard~isos[1:16] + isos.square) #Not significant

#correlation between isolation and jaccard was significant, so plot it
ggplot()+
  geom_point(aes(x = isos[1:16], y = jaccard), size = 2)+
  geom_smooth(aes(x = isos[1:16], y = jaccard), method = 'loess', span = 1.3, color = "black")+
  theme_bw()+
  labs(x = "Isolation (scaled)", y = "Jaccard Index")+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 18),
        panel.grid = element_blank())

# ggsave(filename = "distjaccard.tiff", dpi = 600, width = 85, height = 60, 
#        scale = 1.5, units = "mm")

# Is the relationship between the jaccard index and isolation due to nestedness or turnover?
turns<-lapply(turn.mat, function(x)x[,1])
turns2<-do.call(rbind, turns)

nests<-lapply(turn.mat, function(x)x[,2])
nests2<-do.call(rbind, nests)

trnvr<-data.frame(apply(turns2, 2, mean), apply(nests2, 2, mean))
colnames(trnvr)<-c("turnover", "nestedness")

#Define number of "successes"
dim(trnvr)
length(which(trnvr$nestedness > trnvr$turnover))
#how many 'successes' are in the data set

bin<-binom.test(12, 16, alternative = "greater")
