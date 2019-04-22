#######################################################################
# Constructing and running a MSOM using R2OpenBUGS                    #
#                                                                     #
# Model estimates species-specific occupancy and detection, evaluates #
# effects of covariates on the above, and estimates total species     #
# richness N using data augmentation.                                 #
#                                                                     #
# Data are found in the files 'summer_mamm.csv', the shapefiles       #
# 'Glades_Polygon' and 'Other-glades', and 'all_plant.csv             #
#######################################################################

# Load req'd packages
library(reshape)
library(R2OpenBUGS)
library(abind)
library(rgeos)
library(rgdal)
library(raster)

#Set working directory
# wd <- "c:/users/beasley/dropbox/manuscripts/Mamm IBT/MammIBTRawData"
# setwd(wd)

set.seed(15)

# Set up mammal data ------------------------------------------------------

#Read in mammal data set & cut unwanted columns
all.mamm<-read.csv("summer_mamm.csv", stringsAsFactors = F)
mamm.data<-all.mamm[,c("Date","Glade","Day","Species")]
mamm.data$Occ <- rep(1, dim(mamm.data)[1]) #Add a column to denote presences

#The next few lines are needed to include sites with 0 captures
#If each site had at least 1 capture event, this bit isn't necessary
mamm.data$Date<-as.character(mamm.data$Date)
mamm.data[617,]<-c("9/6/2017",23,1,"Peromyscus attwateri",0) 
mamm.data[618,]<-c("12/7/2017",28,1,"Peromyscus attwateri",0)

mamm.data$Glade <- as.numeric(mamm.data$Glade)

#Create blank vector of number of observed species
unique.species<-as.character(unique(mamm.data$Species))
n <- length(unique.species)

#Create blank vecter of number of sites
unique.sites<-as.character(unique(mamm.data$Glade))
J <- length(unique.sites)

#Create a 3-dim array where j is site, k is day, and i is species
#Use 'melt' function to do so 
melt.that.sucker<- melt(mamm.data,id.var=c("Species", "Glade", "Day"), 
                        measure.var="Occ")
mamm<-cast(melt.that.sucker, Glade~Day~Species)
mamm[mamm>1]<-1

#Add all-zero encounter histories to account for undetected
zeroes<-3 #Number of all-zero species to be added
#Create all-zero matrix with same 'site' and 'day' as encounter matrix
mamm.aug <- array(0, dim=c(dim(mamm)[1],dim(mamm)[2],zeroes)) 

#Combine all-zero matrix and encounter matrix
mamm2<-abind(mamm, mamm.aug, along = 3)

#Create vector of length j indicating number of days at site j
K<-c(rep(4, 18), 3, rep(4, 13))
#n, j, and k will be used within the model itself but it's easier to define them beforehand

# Site-level covariates ------------------------------------------------

#Load shapefile for sampled glades
glades.sampled<-readOGR(dsn= wd, layer="Glades_Polygon")

# #Load shapefile for non-sampled glades ("stepping stones")
# glades.other<-readOGR(dsn = wd, layer = "Other-glades")
# #Change projection to match glades.sampled
# glades.other<-spTransform(glades.other, CRS("+proj=utm +zone=15 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#sampled glades only
glades <- glades.sampled
#all glades
# glades <- union(glades.sampled, glades.other)

#Patch area
polys.area<-gArea(glades.sampled, byid = TRUE)
area<-scale(log(polys.area))
area1<-as.vector(c(area,area))

# saveRDS(area1, "area1.RDS")

#Patch Isolation 

#Calculate centroids
centroids<-gCentroid(glades, byid=TRUE)

#split centroids into management units
rrcent<-rbind(centroids[which(glades@data$LandUnit=="RRSP")])
cacent<-rbind(centroids[which(glades@data$LandUnit=="CMCA")])
mtcent<-rbind(centroids[which(glades@data$LandUnit=="MTNF")])
dmcent<-rbind(centroids[which(glades@data$LandUnit=="DMCA")])

cents <- list(rrcent, cacent, mtcent, dmcent)

#gDistance calculates pairwise distances between objects- in this case, points
cent.dist <- lapply(cents, gDistance, byid = T)

#Split polygons into mgt units
rrpoly<-SpatialPolygons(glades@polygons[which(glades@data$LandUnit == "RRSP")])
capoly<-SpatialPolygons(glades@polygons[which(glades@data$LandUnit == "CMCA")])
mtpoly<-SpatialPolygons(glades@polygons[which(glades@data$LandUnit == "MTNF")])
dmpoly<-SpatialPolygons(glades@polygons[which(glades@data$LandUnit == "DMCA")])

polys <- list(rrpoly, capoly, mtpoly, dmpoly)

polys.dist <- lapply(polys, gDistance, byid = T)

#Get standardized distance for each management unit
standard.dist <- function(x, y){
  stdist <- list()
  for(i in 1:length(x)){
    stdist[[i]] <- x[[i]]/(x[[i]]-y[[i]])
  }
  return(stdist)
}

st.dists <- standard.dist(x = cent.dist, y = polys.dist)

#Change NaN values to 0
clean.dists <- lapply(st.dists, function(x) ifelse(is.nan(x),0,x))

#Calculate mean standardized distance for each glade
mean.dists <- lapply(clean.dists, apply, 2, mean)

#Unlist the above
dist<-unlist(mean.dists)

#Extract sampled glades
# dist <- dist[-c(5:10, 15:19, 23, 29:30)]

#Scale it
scaledist <- scale(dist)
dist1 <- as.vector(c(scaledist,scaledist))

#sampled glades only
# saveRDS(dist1, "dist1.RDS")
#all glades
# saveRDS(dist1, "dist2.RDS")

#Patch shape
polys.perimeter<-gBoundary(glades.sampled, byid=TRUE)
perimeter.lengths<-gLength(polys.perimeter, byid=TRUE)
shape<-(0.282*perimeter.lengths)/sqrt(polys.area) #corrected perimeter/area ratio
shapelog<-scale(log(shape))
shape1<-c(shapelog, shapelog)

# saveRDS(shape1, "shape1.RDS")

#Year
year<-c(rep(1, J/2), rep(2, J/2))

#Vegetation cover (not included in model)
#Read in plant data
plant.data.raw<-read.csv("all_plant.csv")
#Remove unwanted columns
plant.data.mod<-plant.data.raw[,-c(3:5)]
plant.data<-plant.data.mod[,-8]
#Remove NA values
plant.data<-na.omit(plant.data)

#Subset into years
plant.1<-subset(plant.data, Year=="1")
plant.2<-subset(plant.data, Year=="2")

#Get mean values for each site
plant1<-aggregate(plant.1[,3:7], list(plant.1[,2]), mean)
plant2<-aggregate(plant.2[,3:7], list(plant.2[,2]), mean)

#Put years back together
allplant<-rbind(plant1,plant2)

#PCA
pca<-prcomp(allplant[,-1], center=TRUE, scale.=TRUE)
#PC1 and PC2 have eigenvalues greater than 1
#PC1 roughly corresponds to cover type, PC2 to total cover

cover.type <- pca$x[,1]
total.cover <- pca$x[,2]

#Significant relationships present between the pcs and other covariates
summary(lm(cover.type ~ area1))
summary(lm(total.cover ~ shape1))

#Detection covariate: Julian date -------------------------------------------
#Date
year1<-subset(all.mamm, Year=="1")
year2<-subset(all.mamm, Year=="2")

#Convert into date format, then convert to Julian dates
date1<-as.Date(year1[,1], format="%d/%m/%Y")
date2<-as.Date(year2[,1], format="%d/%m/%Y")

jdate1<-julian(date1, origin = as.Date("2016-01-01"))
jdate2<-julian(date2, origin = as.Date("2017-01-01"))

#put it in a matrix
date.day1<-matrix(c(jdate1, year1[,5]), ncol = 2)
date.day2<-matrix(c(jdate2, year2[,5]), ncol = 2)

#In year 1, there was a capture event on day 2 for all glades
#So I could make a vector of those dates and construct matrix from there
twos<-unique(date.day1[,1][which(date.day1[,2]==2)])
dates1<-cbind(twos-1,twos,twos+1,twos+2)

#Year 2 had to be done manually because of fewer captures
ones<-c(181, 182, 185, 185, 152, 153, 158, 159, 192, 146, 146, 206, 192, 205, 
        201, 202)
dates2<-cbind(ones, ones+1, ones+2, ones+3)

dates<-rbind(dates1, dates2)
colnames(dates)<-c("Day 1","Day 2","Day 3","Day 4")

date<-scale(dates)

# saveRDS(date, 'date.RDS')

#Write the model for OpenBUGS --------------------------------------------------
cat("
    model{
    
    #Define prior distributions before creating the priors
      omega~dunif(0,1) 
      
      u.mean~dunif(0,1)
      mu.u<-log(u.mean)-log(1-u.mean)

      v.mean~dunif(0,1)
      mu.v<-log(v.mean)-log(1-v.mean)   
      #u and v are species-specific intercepts on occupancy and detection, respectively

      tau.u ~ dgamma(0.1,0.1)
      tau.v ~ dgamma(0.1,0.1)  #stdevs for u and v, will use when defining priors
    
    #Covs

    mua1 ~ dnorm(0, 0.001)
    tau.a1 ~ dgamma(0.1, 0.1)
    mua2 ~ dnorm(0,0.001)
    tau.a2 ~ dgamma(0.1, 0.1)
    mua3 ~ dnorm(0, 0.001)
    tau.a3 ~ dgamma(0.1, 0.1)
    mua4 ~ dnorm(0, 0.001)
    tau.a4 ~ dgamma(0.1, 0.1)

    mub1 ~ dnorm(0, 0.001)
    tau.b1 ~ dgamma(0.1, 0.1)

    for(i in 1:(n+zeroes)){

      #create your priors from distributions above
        w[i]~dbern(omega)  #mean occupancy species i
        u[i]~dnorm(mu.u, tau.u)  #species level variation in occupancy
        v[i]~dnorm(mu.v, tau.v)  #species level variation in detection

        a1[i]~dnorm(mua1, tau.a1)
        a2[i]~dnorm(mua2, tau.a2)
        a3[i]~dnorm(mua3, tau.a3)
        a4[i]~dnorm(mua4, tau.a4)

        b1[i]~dnorm(mub1, tau.b1)

      #Loop within a loop to estimate occupancy of species i at point j
        for (j in 1:J) {
          logit(psi[j,i]) <- u[i]+a1[i]*area[j]+a2[i]*dist[j]+a3[i]*shape[j]+a4[i]*year[j]
          mu.psi[j,i] <- psi[j,i]*w[i]
          Z[j,i] ~ dbern(mu.psi[j,i])
        
          #Loop within a loop within a loop to estimate detection of i at point j during sampling period k
            for(k in 1:K[j]){
              logit(p[j,k,i]) <-  v[i]+b1[i]*date[j,k]
              mu.p[j,k,i] <- p[j,k,i]*Z[j,i] #The addition of Z means that detecting a species depends on its occupancy state
              mamm[j,k,i] ~ dbern(mu.p[j,k,i])
              }
            }
    }  

    #Estimate total richness (N) by adding observed (n) and unobserved (n0) species
      n0<-sum(w[(n+1):(n+zeroes)])
      N<-n+n0

    #Create a loop to determine point level richness estimates
    for(j in 1:J){
    Nsite[j]<- inprod(Z[j,1:(n+zeroes)],w[1:(n+zeroes)])
    }
    }
    ", file="modmitcovs.txt")

#Send model to OpenBUGS -------------------------------------------------
#Load data and specify parameters
data<-list(n=n, zeroes=zeroes, J=J, K=K, mamm=mamm2, area=area1, dist=dist1,
           shape=shape1, date=date, year=year)

# Data when including all glades
# data <- list(n=n, zeroes=zeroes, J=J, K=K, mamm=mamm2, area=area1, dist=dist2,
#              shape=shape1, date=date, year=year)

params<-list('Z','u','v','mu.u','mu.v','tau.u','tau.v','omega','N','Nsite',
             'a1','a2','a3','a4','b1')

#Specify initial values
init.values<-function(){
  omega.guess<-runif(1, n/(n+zeroes), 1)
  mu.psi.guess<-runif(1, 0.25, 1)
  list(omega=omega.guess, w=c(rep(1,n), rbinom(zeroes,size=1,prob=omega.guess)),
          u=rnorm(n+zeroes), v=rnorm(n+zeroes), 
          Z=matrix(rbinom((n+zeroes)*J, size = 1, prob = mu.psi.guess), 
              nrow = J, ncol = (n+zeroes)),
        a1=rnorm(n+zeroes), a2=rnorm(n+zeroes), a3=rnorm(n+zeroes), 
        a4=rnorm(n+zeroes), 
        b1=rnorm(n+zeroes))
  }

#Run the model. Curse as needed.
model<-bugs(data, inits = init.values, parameters.to.save=params,
            model.file = "modmitcovs.txt", debug = TRUE, n.chains = 3,
            n.burnin = 8500, n.iter = 16000, n.thin = 35)

# saveRDS(model, file = "modelsampledglades.rds")
# saveRDS(model, file = "modelallglades.rds")

model <- readRDS(file = "modelsampledglades.rds")
# model <- readRDS(file = "modelallglades.rds")

#Let's do some summarizing -----------------------------------------
#Occupancy and detection probabilities
species.occ <- model$sims.list$u
species.det <- model$sims.list$v

#Show occupancy and detection estimates for only the observed species (1:n)
psi <- plogis(species.occ[,1:n])
p <- plogis(species.det[,1:n])

occ.matrix <- cbind(apply(psi,2,'mean'),apply(psi,2,'sd'))
colnames(occ.matrix) <- c("mean occupancy", "sd occupancy")
rownames(occ.matrix) <- unique.species
det.matrix <- cbind(apply(p,2,'mean'),apply(p,2,'sd'))
colnames(det.matrix) <- c("mean detection", "sd detection")
rownames(det.matrix) <- unique.species

#See estimates of total richness (N)
N <- model$sims.list$N
mean(N)
summary(N)
plot(table(N))

plot(occ.matrix[,1],det.matrix[,1])

#Get point richness estimates
Nsite<-model$sims.list$Nsite
Nsite<-apply(Nsite, 2, 'mean')

#Regression with area and plot it
ibt<-lm(Nsite~area1 + dist1 + shape1 + year)
plot(Nsite~area1)
