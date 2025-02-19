
#The purpose of this script is to simulate SNP data under assortative mating WITHOUT any LD
#This script started with Simulate.AM5.R, but has been modified to be simpler, faster, and not require LD
#Start date: May 1, 2017

#This script is now an update of 
#/Users/matthewkeller/Google Drive/DriveDocuments/RESEARCH/KongVF/Simulation/Multivariate/AM.FUNCTIONS/OLD/SimulateNOLD.AM.FUNCTIONS_May2020-mck4.R
#It simulates data under multivariate AM.

#This script *-mck7.R now emulates the expectations script in that VT & AM only begin AFTER t0. At t0, only VA and VE influence the trait, which has VY=1.





##############################################################################
#LOAD PACKAGES REQUIRED
##############################################################################

#Load packages & functions needed
require(MASS)
require(fields)
require(matrixcalc)
# scripts_path <- "/projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS/Summary/R-dist"

# r_files <- list.files(scripts_path, pattern = "\\.[rR]$", full.names = TRUE)
# lapply(r_files, source)
#Global options for this script - important to run this or script breaks
#op <-  options(stringsAsFactors=FALSE)

##############################################################################
#END LOAD PACKAGES REQUIRED
##############################################################################











##############################################################################
#FUNCTIONS
##############################################################################


##########################
#F1 - SHORT FUNCTIONS
is.even <- function(x) {x%%2 == 0}
is.odd <- function(x) {x%%2 != 0}


next.smallest <- function(x,value,return.index=TRUE){
x2 <- sort(x,decreasing=FALSE)
if(return.index){
    ans <- max(which(x2<value))
} else {
    ans <- x2[max(which(x2<value))]
}
return(ans)}

cor2cov <- function(X,var1,var2) {
    sd.mat <- matrix(c(sqrt(var1),0,0,sqrt(var2)),nrow=2,byrow=T)
    ans <- sd.mat %*% X %*% sd.mat
    return(ans)}  #only works for 2x2 matrix

##########################







#########################################################
#F2 - REMOVE DUPLICATIONS - Function used in the assort.mate() function below

#----------------------------
# PROPOSED WAY TO CREATE PAIRS SUCH THAT EACH INDIVIDUAL OCCURS IN EXACTLY ONE MATE PAIR AND GIVES DESIRED CORRELATION
#OK, so now deal with duplicate elements in closeM & closeF
#Let's rewrite the function so that the order of individuals in rows cannot matter
#In the old remove.dups function, there were order effects (individuals near end of phenotype matrix are more likely not to mate). Note that this order effect shouldn't matter (because the order of the phenotype matrix should be random anyway), but I've made sure it can't matter in the above function just in case someone ever has an ordered phenotype matrix.
#Finally, I've rewritten this so it takes in both the male and female distance matrix, and in doing so, ensures that we have on average as many complete mate pairs as we had before we reordered the matrix
#DISTm=DM; DISTf=DF; max.ord=100
# remove.dups <- function(DISTm,DISTf,max.ord=100){
#     if(sum(dim(DISTm)==dim(DISTf)) != 2) stop("the dimensions of the two distance matrices must be identical")
#     N.obs <- dim(DISTm)[1]
#     rand.index <- sample(ncol(DISTm),ncol(DISTm),replace=FALSE) 
#     DIST.ran.ord.M <- DISTm[,rand.index] #mixes the columns randomly
#     Dist.Ord.Mat.M <- Dist.Ord.Mat.orig.M <-  apply(DIST.ran.ord.M,2,order)
#     DIST.ran.ord.F <- DISTf[,rand.index] #mixes the columns randomly
#     Dist.Ord.Mat.F <- Dist.Ord.Mat.orig.F <-  apply(DIST.ran.ord.F,2,order)
#     actual.rank.F <- new.closest.F <- actual.rank.M <- new.closest.M <- vector(mode="numeric",length=N.obs)
    
#     for (i in 1:N.obs){ #Note that multi-processor foreach loop will not work here
#         #Males first
#         my.min.M <- Dist.Ord.Mat.M[1,i]
#         k.M <- 0
#         while(my.min.M %in% new.closest.M){
#             k.M <- k.M+1
#             Dist.Ord.Mat.M[1:max.ord,i] <- Dist.Ord.Mat.M[2:(max.ord+1),i]
#             my.min.M <- Dist.Ord.Mat.M[1,i]
#             if (k.M > max.ord) my.min.M <- 1e9 + i #ensures it's a unique value that we can select upon at end
#         } #end while loop
#         new.closest.M[i] <- my.min.M
#         actual.rank.M[i] <- k.M+1
        
#         #Females second
#         my.min.F <- Dist.Ord.Mat.F[1,i]
#         k.F <- 0
#         while(my.min.F %in% new.closest.F){
#             k.F <- k.F+1
#             Dist.Ord.Mat.F[1:max.ord,i] <- Dist.Ord.Mat.F[2:(max.ord+1),i]
#             my.min.F <- Dist.Ord.Mat.F[1,i]
#             if (k.F > max.ord) my.min.F <- 1e9 + i #ensures it's a unique value that we can select upon at end
#         } #end while loop
#         new.closest.F[i] <- my.min.F
#         actual.rank.F[i] <- k.F+1
#     } #end for loop
    
#     new.closest.M <- new.closest.M[order(rand.index)] #gets them back to original order
#     actual.rank.M <- actual.rank.M[order(rand.index)] #gets them back to original order
#     new.closest.M[new.closest.M > 1e9] <- NA
#     actual.rank.M[is.na(new.closest.M)] <- NA
#     Dist.Ord.Mat.orig.M <- Dist.Ord.Mat.orig.M[,order(rand.index)]
    
#     new.closest.F <- new.closest.F[order(rand.index)] #gets them back to original order
#     actual.rank.F <- actual.rank.F[order(rand.index)] #gets them back to original order
#     new.closest.F[new.closest.F > 1e9] <- NA
#     actual.rank.F[is.na(new.closest.F)] <- NA
#     Dist.Ord.Mat.orig.F <- Dist.Ord.Mat.orig.F[,order(rand.index)]
    
#     return(list(closestM=new.closest.M,rankM=actual.rank.M,original.closestM=Dist.Ord.Mat.orig.M[1,],closestF=new.closest.F,rankF=actual.rank.F,original.closestF=Dist.Ord.Mat.orig.F[1,]))
# }

# Xuanyu 11/13/24 - modify the function see if it uses less ram
remove.dups <- function(DISTm,max.ord=100){
    N.obs <- dim(DISTm)[1]
    rand.index <- sample(ncol(DISTm),ncol(DISTm),replace=FALSE) 
    #DIST.ran.ord.M <- DISTm[,rand.index] #mixes the columns randomly
    Dist.Ord.Mat.M <-  apply(DISTm[,rand.index],2,order)
    new.closest.M <- vector(mode="numeric",length=N.obs)
    
    for (i in 1:N.obs){ #Note that multi-processor foreach loop will not work here
        #Males first
        my.min.M <- Dist.Ord.Mat.M[1,i]
        k.M <- 0
        while(my.min.M %in% new.closest.M){
            k.M <- k.M+1
            Dist.Ord.Mat.M[1:max.ord,i] <- Dist.Ord.Mat.M[2:(max.ord+1),i]
            my.min.M <- Dist.Ord.Mat.M[1,i]
            if (k.M > max.ord) my.min.M <- 1e9 + i #ensures it's a unique value that we can select upon at end
        } #end while loop
        new.closest.M[i] <- my.min.M
        #actual.rank.M[i] <- k.M+1
    } #end for loop
    
    new.closest.M <- new.closest.M[order(rand.index)] #gets them back to original order
    #actual.rank.M <- actual.rank.M[order(rand.index)] #gets them back to original order
    new.closest.M[new.closest.M > 1e9] <- NA
    #actual.rank.M[is.na(new.closest.M)] <- NA
    #Dist.Ord.Mat.orig.M <- Dist.Ord.Mat.orig.M[,order(rand.index)]

    
    return(new.closest.M)
}



#########################################################







#########################################################
#F3 - ASSORTATIVE MATING FUNCTION

#PHENDATA=PHEN;  MATCOR=PHENO.MATE.CUR; POPSIZE=POP.SIZE; avoid.inbreeding=AVOID.INB
assort.mate <- function(PHENDATA, MATCOR, POPSIZE, avoid.inbreeding=TRUE){

#probability that each individual marries for generation currun; capped at .95
(size.mate <- nrow(PHENDATA))

### MARRIAGEABLE PEOPLE - this section creates an equal number of males and females who will be paired off below
#We need to make the number of marriageable people the same in females & males
(num.males.mate  <- sum(PHENDATA[,"MALE"]))
(num.females.mate <- size.mate-num.males.mate)

### MATING VIA PRIMARY ASSORTMENT - ANCESTORS
#split up the PHENDATA matrix
males.PHENDATA <- PHENDATA[PHENDATA[,"MALE"]==1,]
females.PHENDATA <- PHENDATA[PHENDATA[,"MALE"]==0,]

#make exactly equal numbers of males & females
if(num.males.mate>num.females.mate) {males.PHENDATA <- males.PHENDATA[- sample(1:nrow(males.PHENDATA),1),]}
if(num.males.mate<num.females.mate) {females.PHENDATA <- females.PHENDATA[- sample(1:nrow(females.PHENDATA),1),]}

## Xuanyu 02/15/24: The correlation matrix are not always positive definite. Adapt some code to use SVD to make it positive definite

#print(MATCOR)

#cat(is.positive.definite(MATCOR),"\n")
if (!is.positive.definite(MATCOR)) {
    cat("Matrix is not positive definite. Using SVD to make it positive definite.\n")
    s <- svd(MATCOR)
    V <- s$v
    D <- sqrt(zapsmall(diag(s$d)))
    X <- (V %*% D) %*% matrix(rnorm(4*1e5), 4) 
    MATCOR <- cov(t(X))
}
#print(MATCOR)
#Create simulated m1, m2, f1, f2 data with a specified correlation structure. This is the template correlated data
Xsim <- mvrnorm(nrow(males.PHENDATA),rep(0,4),MATCOR,empirical=TRUE)
XsimM <- Xsim[,1:2]
XsimF <- Xsim[,3:4]
MATCOR ; round(cor(Xsim),4) #CHECK - should be the same

#get distances between all Xm pairs and XsimM pairs & similarly for females
#NOTE: if there are more males than females, remove a number of males at random s.t. the two are equal, and vice-versa.
DM <- rdist(scale(males.PHENDATA[,c("Y1","Y2")]),XsimM) #this means that XsimM goes along y-dimension; Xm along x-dimension
cat("DM\n")

N.obs <- dim(DM)[1]
rand.index <- sample(ncol(DM),ncol(DM),replace=FALSE) 
#DIST.ran.ord.M <- DISTm[,rand.index] #mixes the columns randomly
DM <-  apply(DM[,rand.index],2,order)
new.closest.M <- vector(mode="numeric",length=N.obs)
max.ord = 100
for (i in 1:N.obs){ #Note that multi-processor foreach loop will not work here
    #Males first
    my.min.M <- DM[1,i]
    k.M <- 0
    while(my.min.M %in% new.closest.M){
        k.M <- k.M+1
        DM[1:max.ord,i] <- DM[2:(max.ord+1),i]
        my.min.M <- DM[1,i]
        if (k.M > max.ord) my.min.M <- 1e9 + i #ensures it's a unique value that we can select upon at end
    } 
    new.closest.M[i] <- my.min.M

} #end for loop

new.closest.M <- new.closest.M[order(rand.index)] #gets them back to original order
new.closest.M[new.closest.M > 1e9] <- NA

# delete DM from memory 
rm(DM)
Xm.ord2 <- males.PHENDATA[new.closest.M,]

DF <- rdist(scale(females.PHENDATA[,c("Y1","Y2")]),XsimF) 
cat("DF\n")
N.obs <- dim(DF)[1]
rand.index <- sample(ncol(DF),ncol(DF),replace=FALSE) 
#DIST.ran.ord.M <- DISTm[,rand.index] #mixes the columns randomly
DF <-  apply(DF[,rand.index],2,order)
new.closest.F <- vector(mode="numeric",length=N.obs)
max.ord = 100
for (i in 1:N.obs){ #Note that multi-processor foreach loop will not work here
    #Males first
    my.min.F <- DF[1,i]
    k.F <- 0
    while(my.min.F %in% new.closest.F){
        k.F <- k.F+1
        DF[1:max.ord,i] <- DF[2:(max.ord+1),i]
        my.min.F <- DF[1,i]
        if (k.F > max.ord) my.min.F <- 1e9 + i #ensures it's a unique value that we can select upon at end
    } 
    new.closest.F[i] <- my.min.F

} #end for loop

new.closest.F <- new.closest.F[order(rand.index)] #gets them back to original order
new.closest.F[new.closest.F > 1e9] <- NA
# delete DF from memory
rm(DF)
Xf.ord2 <- females.PHENDATA[new.closest.F,]

#remove.dups 
# TOT <- remove.dups(DM,DF,max.ord=100) 
# Xm.ord2 <- males.PHENDATA[TOT$closestM,]
# Xf.ord2 <- females.PHENDATA[TOT$closestF,]
has.nas <- is.na(Xm.ord2[,"ID"]) | is.na(Xf.ord2[,"ID"])

#Recreate [fe]males.PHENDATA 
males.PHENDATA <- Xm.ord2[!has.nas,]
females.PHENDATA <- Xf.ord2[!has.nas,]
cat("phenotype data got \n")
cov(cbind(males.PHENDATA[,c("Y1","Y2")],females.PHENDATA[,c("Y1","Y2")])) #Check

###AVOID INBREEDING
if (avoid.inbreeding==TRUE){

#for now, we just drop the inbreeding pairs; an alternative is to remate them
no.sib.inbreeding <- males.PHENDATA[,"Father.ID"]!=females.PHENDATA[,"Father.ID"]
no.cousin.inbreeding <- {males.PHENDATA[,"Fathers.Father.ID"]!=females.PHENDATA[,"Fathers.Father.ID"]&
                      males.PHENDATA[,"Fathers.Father.ID"]!=females.PHENDATA[,"Mothers.Father.ID"]&
                      males.PHENDATA[,"Mothers.Father.ID"]!=females.PHENDATA[,"Fathers.Father.ID"]&
                      males.PHENDATA[,"Mothers.Father.ID"]!=females.PHENDATA[,"Mothers.Father.ID"]&
                      males.PHENDATA[,"Fathers.Mother.ID"]!=females.PHENDATA[,"Fathers.Mother.ID"]&
                      males.PHENDATA[,"Fathers.Mother.ID"]!=females.PHENDATA[,"Mothers.Mother.ID"]&
                      males.PHENDATA[,"Mothers.Mother.ID"]!=females.PHENDATA[,"Fathers.Mother.ID"]&
                      males.PHENDATA[,"Mothers.Mother.ID"]!=females.PHENDATA[,"Mothers.Mother.ID"]}
no.inbreeding <- (no.sib.inbreeding & no.cousin.inbreeding)


### CREATE FINAL MARRIED PAIRS WHO AREN'T INBRED
#for now, we just drop inbred pairs & create the male & female files (both in order of married pairs)
     males.PHENDATA <- males.PHENDATA[no.inbreeding,]
     females.PHENDATA <- females.PHENDATA[no.inbreeding,]

} ###End avoid inbreeding


### CREATE COR.SPOUSES & "SPOUSE.ID" VARIABLE
males.Spouse.ID <- females.PHENDATA[,'ID']
females.Spouse.ID <- males.PHENDATA[,'ID']

#Create EXACTLY number of offspring desired (note it's possible that, rarely, number offspring won't be exactly that desired, if the # it is off initially is greater than or less than the number of rows in males.PHENDATA)
num.offspring <- rpois(nrow(males.PHENDATA),lambda=POPSIZE/nrow(males.PHENDATA))
(nof <- sum(num.offspring))
if (nof < POPSIZE){num.offspring <- num.offspring + sample(c(rep(1,POPSIZE-nof),rep(0,nrow(males.PHENDATA)-POPSIZE+nof)),nrow(males.PHENDATA),replace=FALSE)}
if (nof > POPSIZE){
subtractor <- rep(0,nrow(males.PHENDATA))
(which.pos <- which(num.offspring>0))
(which.subt <- sample(which.pos,nof-POPSIZE,replace=FALSE))
subtractor[which.subt] <- -1
num.offspring <- num.offspring+subtractor}
#summary(num.offspring) #CHECK

#Also put the spousal info into PHENDATA
males.PHENDATA <- cbind(males.PHENDATA,males.Spouse.ID,num.offspring)
females.PHENDATA <- cbind(females.PHENDATA,females.Spouse.ID,num.offspring)


#make a list of what is returned from the function
list(males.PHENDATA=males.PHENDATA,females.PHENDATA=females.PHENDATA)
}
#########################################################





#########################################################
#F4 Reproduce Function

#MATES=MATES; XO=XO; XL=XL; PHEN=PHEN; CV.INFO=CV.INFO; CUR.GEN=CUR.GEN; cove.mat=cove.mat; fmat=f.mat; amat=a.mat; dmat=delta.mat; cor.list=am.list; covy.mat=COVY; k2.matrix=k2.matrix
reproduce <- function(MATES=MATES, XO=XO, XL=XL, PHEN=PHEN, CV.INFO=CV.INFO, CUR.GEN=CUR.GEN, cove.mat=cove.mat, fmat=f.mat, amat=a.mat, dmat=delta.mat, cor.list=am.list, covy.mat=COVY, k2.matrix=k2.matrix){

#Create male & female data in mating order for only those who had offspring
males.PHEN <- MATES$males.PHENDATA[MATES$males.PHENDATA[,"num.offspring"]>0,]
females.PHEN <- MATES$females.PHENDATA[MATES$females.PHENDATA[,"num.offspring"]>0,]

#Put male & female genotypes in order of those above
males.index <- match(males.PHEN[,"ID"],PHEN[,"ID"])
males.gentp1.obs <- XO[males.index,]
males.gentp1.lat <- XL[males.index,]
females.index <- match(females.PHEN[,"ID"],PHEN[,"ID"])
females.gentp1.obs <- XO[females.index,]
females.gentp1.lat <- XL[females.index,]

#CHECK - should be 1; it is
#bvm <- males.gentp1 %*% CV.INFO$alpha
#cor(bvm,males.PHEN[,"BV"])
#bvf <- females.gentp1 %*% CV.INFO$alpha
#cor(bvf,females.PHEN[,"BV"])

#Create 1 row per new offspring
males.rep <- rep(1:nrow(males.PHEN),times=males.PHEN[,"num.offspring"])
males.gentp.obs <- males.gentp1.obs[males.rep,]
males.gentp.lat <- males.gentp1.lat[males.rep,]
females.rep <- rep(1:nrow(females.PHEN),times=females.PHEN[,"num.offspring"])
females.gentp.obs <- females.gentp1.obs[females.rep,]
females.gentp.lat <- females.gentp1.lat[females.rep,]


#Create matrices for which of the 2 alleles are passed on at random (male haplotypes) - TRANSMITTED, OBSERVED
male.adder.obs <- matrix(sample(x=c(0,1),size=nrow(males.gentp.obs)*ncol(males.gentp.obs),replace=TRUE),nrow=nrow(males.gentp.obs),ncol=ncol(males.gentp.obs))
XM1.obs <- (males.gentp.obs==1)*1
XM1.adder.obs <- male.adder.obs*XM1.obs #give one or other allele IF you are heterozygous
XM2.obs <- (males.gentp.obs==2)*1
males.haps.obs <- XM1.adder.obs+XM2.obs

#Create matrices for which of the 2 alleles are passed on at random (male haplotypes) - TRANSMITTED, LATENT
male.adder.lat <- matrix(sample(x=c(0,1),size=nrow(males.gentp.lat)*ncol(males.gentp.lat),replace=TRUE),nrow=nrow(males.gentp.lat),ncol=ncol(males.gentp.lat))
XM1.lat <- (males.gentp.lat==1)*1
XM1.adder.lat <- male.adder.lat*XM1.lat #give one or other allele IF you are heterozygous
XM2.lat <- (males.gentp.lat==2)*1
males.haps.lat <- XM1.adder.lat+XM2.lat


#Create matrices for the NON-TRANSMITTED alleles (male haplotypes) - NONTRANSMITTED, OBSERVED
male.NT.adder.obs <- (male.adder.obs*-1)+1
XM1.NT.adder.obs <- male.NT.adder.obs*XM1.obs #give one or other allele IF you are heterozygous
males.NT.haps.obs <- XM1.NT.adder.obs+XM2.obs

#Create matrices for the NON-TRANSMITTED alleles (male haplotypes) - NONTRANSMITTED, LATENT
male.NT.adder.lat <- (male.adder.lat*-1)+1
XM1.NT.adder.lat <- male.NT.adder.lat*XM1.lat #give one or other allele IF you are heterozygous
males.NT.haps.lat <- XM1.NT.adder.lat+XM2.lat


#Create matrices for which of the 2 alleles are passed on at random (female haplotypes) - TRANSMITTED, OBSERVED
female.adder.obs <- matrix(sample(x=c(0,1),size=nrow(females.gentp.obs)*ncol(females.gentp.obs),replace=TRUE),nrow=nrow(females.gentp.obs),ncol=ncol(females.gentp.obs))
XF1.obs <- (females.gentp.obs==1)*1
XF1.adder.obs <- female.adder.obs*XF1.obs
XF2.obs <- (females.gentp.obs==2)*1
females.haps.obs <- XF1.adder.obs+XF2.obs

#Create matrices for which of the 2 alleles are passed on at random (female haplotypes) - TRANSMITTED, LATENT
female.adder.lat <- matrix(sample(x=c(0,1),size=nrow(females.gentp.lat)*ncol(females.gentp.lat),replace=TRUE),nrow=nrow(females.gentp.lat),ncol=ncol(females.gentp.lat))
XF1.lat <- (females.gentp.lat==1)*1
XF1.adder.lat <- female.adder.lat*XF1.lat
XF2.lat <- (females.gentp.lat==2)*1
females.haps.lat <- XF1.adder.lat+XF2.lat


#Create matrices for the NON-TRANSMITTED alleles (female haplotypes) - NONTRANSMITTED, OBSERVED
female.NT.adder.obs <- (female.adder.obs*-1)+1
XF1.NT.adder.obs <- female.NT.adder.obs*XF1.obs #give one or other allele IF you are heterozygous
females.NT.haps.obs <- XF1.NT.adder.obs+XF2.obs

#Create matrices for the NON-TRANSMITTED alleles (female haplotypes) - NONTRANSMITTED, LATENT
female.NT.adder.lat <- (female.adder.lat*-1)+1
XF1.NT.adder.lat <- female.NT.adder.lat*XF1.lat #give one or other allele IF you are heterozygous
females.NT.haps.lat <- XF1.NT.adder.lat+XF2.lat


#Create new genotypes for offspring
XO <- males.haps.obs + females.haps.obs
XL <- males.haps.lat + females.haps.lat

#Create 1 row per new offspring from the phenotype matrix
fathers.PHEN <- males.PHEN[males.rep,]  
mothers.PHEN <- females.PHEN[females.rep,] 

#Create Phenotypes & haplotype PRS's - OBSERVED
AO <- XO %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])  #this is the observed (PGS) breeding values for traits 1 and 2
TPO <- males.haps.obs %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
TMO <- females.haps.obs %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
NTPO <- males.NT.haps.obs %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
NTMO <- females.NT.haps.obs %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
BV.NT.O <- NTPO+NTMO


#Create Phenotypes & haplotype breeding values - LATENT
AL <- XL %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])  #this is the latent (PGS) breeding values for traits 1 and 2
TPL <- males.haps.lat %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
TML <- females.haps.lat %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
NTPL <- males.NT.haps.lat %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
NTML <- females.NT.haps.lat %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])
BV.NT.L <- NTPL+NTML

#Create components of Y (the above multiplied by their path coefficients) and Y at time 0
(CURRENT.N <- nrow(BV.NT.O))
F <- as.matrix(fathers.PHEN[,c("Y1","Y2")]) %*% t(fmat) + as.matrix(mothers.PHEN[,c("Y1","Y2")]) %*% t(fmat)
E <- mvrnorm(CURRENT.N,c(0,0),Sigma=cove.mat,empirical=T)
AOy <- AO %*% t(dmat)
ALy <- AL %*% t(amat)
Fy <- F %*% diag(2)
Ey <- E %*% diag(2)
Y <- (AOy + ALy + Fy + Ey)
#CHECK
#BV2 <- TPO+TMO #exactly same as AO
#cor(cbind(AO,BV2))
#cov(cbind(TPO,NTPO,TMO,NTMO,AO,BV.NT.O))
#var(AOy);var(Fy);var(Ey);var(Y)

#Create relative information for this generation
ID <- sample(1000000:9999999,size=CURRENT.N,replace=FALSE)
Father.ID <- males.PHEN[,"ID"][males.rep]
Mother.ID <- females.PHEN[,"ID"][females.rep]
Fathers.Father.ID <- males.PHEN[,"Father.ID"][males.rep]
Fathers.Mother.ID <- males.PHEN[,"Mother.ID"][males.rep]
Mothers.Father.ID <- females.PHEN[,"Father.ID"][females.rep]
Mothers.Mother.ID <- females.PHEN[,"Mother.ID"][females.rep]
TOT.ID <- cbind(ID,Father.ID,Mother.ID,Fathers.Father.ID,Fathers.Mother.ID,Mothers.Father.ID,Mothers.Mother.ID)
Father.Y <- males.PHEN[males.rep,c("Y1","Y2")]
Father.F <- males.PHEN[males.rep,c("F1","F2")]
Mother.Y <- females.PHEN[females.rep,c("Y1","Y2")]
Mother.F <- females.PHEN[females.rep,c("F1","F2")]


#Create gender for this gen, making an equal number of males & females
if (is.even(CURRENT.N)) {SEXVEC <- c(rep(0,CURRENT.N/2),rep(1,CURRENT.N/2))}
if (is.odd(CURRENT.N)) {SEXVEC <- c(rep(0,floor(CURRENT.N/2)),rep(1,floor(CURRENT.N/2)),0)}
MALE <- sample(SEXVEC,size=CURRENT.N,replace=FALSE)

#Create Phenotype Matrix
PHEN <- cbind(TOT.ID,MALE,
AO,AL,F,E,
AOy,ALy,Fy,Ey,Y,
BV.NT.O,TPO,TMO,NTPO,NTMO,
BV.NT.L,TPL,TML,NTPL,NTML,
Father.Y,Father.F,Mother.Y,Mother.F)

colnames(PHEN) <- c('ID','Father.ID','Mother.ID','Fathers.Father.ID','Fathers.Mother.ID','Mothers.Father.ID','Mothers.Mother.ID','MALE',
'AO1','AO2','AL1','AL2','F1','F2','E1','E2',
'AOy1','AOy2','ALy1','ALy2','Fy1','Fy2','Ey1','Ey2','Y1','Y2',
'BV.NT.O1','BV.NT.O2','TPO1','TPO2','TMO1','TMO2','NTPO1','NTPO2','NTMO1','NTMO2',
'BV.NT.L1','BV.NT.L2','TPL1','TPL2','TML1','TML2','NTPL1','NTPL2','NTML1','NTML2',
'Y1P','Y2P','F1P','F2P','Y1M','Y2M','F1M','F2M')




#return PHEN & genotypes
list(PHEN=PHEN,XO=XO,XL=XL)
}
#########################################################








#########################################################
#F5 - Assortative mating function

#User supplied to function, for checking only (comment out)
#CV.INFO=cv.info; NUM.GENERATIONS=num.gen; POP.SIZE=pop.size; AVOID.INB=avoid.inb; SAVE.EACH.GEN=save.history; SAVE.COVS=save.covariances; SEED=seed; cove.mat=cove.mat; fmat=f.mat; amat=a.mat; dmat=delta.mat; cor.list=am.list; covy.mat=COVY; k2.matrix=k2.matrix

AM.SIMULATE <- function(CV.INFO, NUM.GENERATIONS, POP.SIZE, AVOID.INB, SAVE.EACH.GEN, SAVE.COVS, SEED, cove.mat, fmat, amat, dmat, cor.list, covy.mat, k2.matrix){
    
cat("Starting simulation\n")

###################
#A1 IMPLIED VARIABLES
NUM.CVs <- nrow(CV.INFO)
POP.VECTOR <- rep(POP.SIZE,NUM.GENERATIONS)  #population size over time
###################


###################
#A2 Create Genotypes and Phenotypes for GEN0

#Set seed if needed
if (SEED != 0) {
	set.seed(SEED)
}



#TEMP
#num.cvs <- 10
#min.maf <- .1
#max.maf <- .5
#rg <- .5
#maf.vector <- runif(num.cvs,min.maf,max.maf)  #Can change the distribution of MAFs here
#gentp.var <- maf.vector*(1-maf.vector)*2

#alphas.pre <- mvrnorm(num.cvs,c(0,0),matrix(c(1,rg,rg,1),nrow=2),empirical=TRUE)
#CAREFUL: when you specify empirical=TRUE, sum(alphas.pre^2) is NOT num.cvs. It is num.cvs-1
#alphas <- alphas.pre * cbind(sqrt(1/((num.cvs-1)*gentp.var)),sqrt(1/((num.cvs-1)*gentp.var)))

#cv.info <- data.frame(maf=maf.vector,alpha1=alphas[,1],alpha2=alphas[,2]) #we'll use this for both the observed and latent
#CV.INFO <- cv.info
#NUM.CVs <- nrow(CV.INFO)
#*******


#Create observed and latent (uncorrelated) Genotypes of GEN0. They share MAFs and (below) alphas, but they are nevertheless indep
NUM.GENTPS <- NUM.CVs*POP.SIZE
XO <- matrix(rbinom(NUM.GENTPS,size=2,prob=CV.INFO$maf),nrow=POP.SIZE,ncol=NUM.CVs,byrow=TRUE)
XL <- matrix(rbinom(NUM.GENTPS,size=2,prob=CV.INFO$maf),nrow=POP.SIZE,ncol=NUM.CVs,byrow=TRUE)
#Check
#obs.mafs <- apply(XO,2,mean)/2
#plot(CV.INFO$maf,obs.mafs)
#exp.vr <- CV.INFO$maf*2*(1-CV.INFO$maf)
#obs.vr <- diag(var(XO))
#plot(exp.vr,obs.vr)
#summary(lm(obs.vr~exp.vr))
cat("Check1 \n")
#Create Phenotypes of GEN0; VA will be 1 at T0; other components will be relative to VA @t0
AO <- XO %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])  #this is the observed (PGS) breeding values for traits 1 and 2; they should  have var=1 and cov = rg at t0

AL <- XL %*% as.matrix(CV.INFO[,c("alpha1","alpha2")])  #this is the latent (PGS) breeding values for traits 1 and 2; they should  have var=1 and cov = rg at t0. Note that, for convenience and ensuring no basic difference in genetic architectures between latent and observed, we reuse the same alphas. This does NOT create a corr between them because the genotypes are independent (XO indep of XL)
var(cbind(AO,AL)) #Check - the rg should be there on two off-diagonals. The bottom 2x2 square should be ~0.

#Create equivalent of k*2 and j*2 matrices 
(var.bvo.gen0 <- var(AO)) #this should be exactly k2.matrix as m -> infinity
(var.bvl.gen0 <- var(AL)) #this should be exactly j2.matrix as m -> infinity
k2.matrix

#Create components of Y (the above multiplied by their path coefficients) and Y at time 0
#For F, we need to go ahead and account for any covariance due to pleiotropy and cov(E). Otherwise, the phenotypic variances at time 0 won't equal 1.
F <- mvrnorm(POP.SIZE,c(0,0),covy.mat,empirical=T) %*% t(fmat) + mvrnorm(POP.SIZE,c(0,0),covy.mat,empirical=T) %*% t(fmat)
E <- mvrnorm(POP.SIZE,c(0,0),Sigma=cove.mat,empirical=T)
AOy <- AO %*% t(dmat)
ALy <- AL %*% t(amat)
Fy <- F %*% diag(2)
Ey <- E %*% diag(2)
Y <- (AOy + ALy + Fy + Ey)

#Observed variances & check they're right
(var.e <- var(Ey)); cove.mat
(var.f <- var(Fy)); 2 * (fmat %*% covy.mat %*% t(fmat))
(var.ao <- var(AOy)); delta.mat %*% k2.matrix %*% t(delta.mat) 
(var.al <- var(ALy)); a.mat %*% k2.matrix %*% t(a.mat)  
(var.y <- var(Y)); covy.mat #should be close but not exact bc var.ao and var.al are approximates
cat("Check2 \n")
#Create relative information for GEN0
TOT.ID.VEC <- matrix(sample(10000000:99999999,size=POP.SIZE*7,replace=FALSE),nrow=POP.SIZE,ncol=7)

#Create gender for GEN0, making an equal number of males & females
if (is.even(POP.SIZE)) {SEXVEC <- c(rep(0,POP.SIZE/2),rep(1,POP.SIZE/2))}
if (is.odd(POP.SIZE)) {SEXVEC <- c(rep(0,floor(POP.SIZE/2)),rep(1,floor(POP.SIZE/2)),0)}
MALE <- sample(SEXVEC,size=POP.SIZE,replace=FALSE)
cat("Check3 \n")
#Create Phenotype Matrix
na.mat <- matrix(NA,nrow=length(MALE),ncol=28)
PHEN <- cbind(TOT.ID.VEC,MALE,AO,AL,F,E,AOy,ALy,Fy,Ey,Y,na.mat)
colnames(PHEN) <- c('ID','Father.ID','Mother.ID','Fathers.Father.ID','Fathers.Mother.ID','Mothers.Father.ID','Mothers.Mother.ID','MALE',
                    'AO1','AO2','AL1','AL2','F1','F2','E1','E2',
                    'AOy1','AOy2','ALy1','ALy2','Fy1','Fy2','Ey1','Ey2','Y1','Y2',
                    'BV.NT.O1','BV.NT.O2','TPO1','TPO2','TMO1','TMO2','NTPO1','NTPO2','NTMO1','NTMO2',
                    'BV.NT.L1','BV.NT.L2','TPL1','TPL2','TML1','TML2','NTPL1','NTPL2','NTML1','NTML2',
                    'Y1P','Y2P','F1P','F2P','Y1M','Y2M','F1M','F2M')

#Potentially save GEN0
HISTORY <- list()
if (SAVE.EACH.GEN){
HISTORY$MATES[[1]] <- PHEN
HISTORY$PHEN[[1]] <- PHEN
HISTORY$XO[[1]] <- XO
HISTORY$XL[[1]] <- XL
}

COVS <- list()
if (SAVE.COVS){
COVS[[1]] <- NULL}

cat("check4 \n")
#Save basic information
SUMMARY.RES <- vector(mode="list",length=NUM.GENERATIONS+1)
SUMMARY.RES[[1]] <-      list(GEN=0,
                              NUM.CVs=NUM.CVs,
                              MATE.COR=cor.list[[1]],
                              POPSIZE=POP.SIZE,
                              VAO=var(PHEN[,c('AOy1','AOy2')]), #new0 - 
                              VAL=var(PHEN[,c('ALy1','ALy2')]), #new1  - 
                              VF=var(PHEN[,c('F1','F2')]),
                              VE=var(PHEN[,c('E1','E2')]),
                              VP=var(PHEN[,c('Y1','Y2')]),
                              h2=(var(PHEN[,c('AOy1','AOy2')]) + var(PHEN[,c('ALy1','ALy2')])) /var(PHEN[,c('Y1','Y2')]), #new0 - 
                              h2.obs=(var(PHEN[,c('AOy1','AOy2')])) /var(PHEN[,c('Y1','Y2')]), #new1 - 
                              h2.lat=(var(PHEN[,c('ALy1','ALy2')])) /var(PHEN[,c('Y1','Y2')]), #new1 - 
                              
                              covY=NA,
                              covG=NA,
                              covH=NA,
                              covI=NA,
                              w=NA,
                              q=NA,
                              covF=NA,
                              covE=NA,
                              
                              hapsO.covs= NA,
                              hapsL.covs= NA,
                              
                              omega=NA,
                              gamma=NA,
                              thetaNT=NA,
                              thetaT=NA
)


###################



##########################
#A3 - Loop through the generations, beginning with generation 1
for (CUR.GEN in  1:NUM.GENERATIONS){

#WILDCARDS THIS GEN
POP.SIZE <- POP.VECTOR[CUR.GEN] #scalar - the pop size for the children of this generation

#NOTE: this is a tricky part that I may have wrong. Issue = the phenotypic corr bw Y1 & Y2 will change across generations due to AM/VT. But if cor(Y1,Y2) increases, in order to keep the 4 AM correlations constant, the actual values of mu will have to decrease. This was also the case in the original model and simulation. So I think this (having the 2x2 AM correlations be constant even though the cor(Y1,Y2) increase) is OK, but something to keep an eye on.
(MATE.CUR <- cor.list[[CUR.GEN]]) #matrix - correlation this gen based on original var-covar matrix of Y
PHENO.MATE.CUR <- diag(4)
PHENO.MATE.CUR[1,2] <- PHENO.MATE.CUR[2,1] <- PHENO.MATE.CUR[4,3] <- PHENO.MATE.CUR[3,4] <- cor(PHEN[,c("Y1","Y2")])[1,2] #Here's where we dynamically change cor(Y1,Y2)
PHENO.MATE.CUR[1:2,3:4] <- MATE.CUR
PHENO.MATE.CUR[3:4,1:2] <- t(MATE.CUR)
PHENO.MATE.CUR
##########################

cat("Beginning generation",CUR.GEN,"\n")
###################
#A4 Assortatively mate the people
#MATES FOR EACH GROUP THIS GENERATION
#PHENDATA=PHEN;  MATCOR=PHENO.MATE.CUR; POPSIZE=POP.SIZE; avoid.inbreeding=AVOID.INB
MATES <- assort.mate(PHENDATA=PHEN, MATCOR=PHENO.MATE.CUR, POPSIZE=POP.SIZE, avoid.inbreeding=AVOID.INB)
#cor(cbind(MATES$males.PHENDATA[,c('Y1','Y2')],MATES$females.PHENDATA[,c('Y1','Y2')])) #CHECK
###################
cat("Generation",CUR.GEN,"assortative mating done\n")

###################
#A5 Have mates reproduce and create new genotypes and phenotypes of offspring for this generation
#Reproduce
#MATES=MATES,XO=XO,XL=XL,PHEN=PHEN,CV.INFO=CV.INFO,CUR.GEN=CUR.GEN,cove.mat=cove.mat; fmat=f.mat; amat=a.mat; dmat=delta.mat; cor.list=am.list; covy.mat=COVY; k2.matrix=k2.matrix)
OFFSPRING <- reproduce(MATES=MATES,XO=XO,XL=XL,PHEN=PHEN,CV.INFO=CV.INFO,CUR.GEN=CUR.GEN,cove.mat=cove.mat, fmat=f.mat, amat=a.mat, dmat=delta.mat, cor.list=am.list, covy.mat=COVY, k2.matrix=k2.matrix)

cat("Generation",CUR.GEN,"reproduction done\n")
XO <- OFFSPRING$XO
XL <- OFFSPRING$XL
PHEN <- OFFSPRING$PHEN
###################


###################
#A6 Finish loop

covs <- cov(PHEN[,c('TPO1','TMO1','NTPO1','NTMO1','TPL1','TML1','NTPL1','NTML1','TPO2','TMO2','NTPO2','NTMO2','TPL2','TML2','NTPL2','NTML2','AO1','AO2','AL1','AL2','F1','F2','E1','E2','BV.NT.O1','BV.NT.O2','Y1','Y2','Y1P','Y2P','Y1M','Y2M','F1P','F2P','F1M','F2M')])

hapsO.covs <- cov(PHEN[,c('TPO1','TMO1','NTPO1','NTMO1','TPO2','TMO2','NTPO2','NTMO2')])
(hapsO.covs1 <- hapsO.covs[c(1:4),c(1:4)])
(hapsO.covs2 <- hapsO.covs[c(5:8),c(5:8)])
(hapsO.covs12 <- hapsO.covs[5:8,1:4])
g11pre <- hapsO.covs1 - diag(4)*(k2.matrix/2)[1,1]
(g11 <- g11pre[lower.tri(g11pre,diag=TRUE)])
g22pre <- hapsO.covs2 - diag(4)*(k2.matrix/2)[2,2]
(g22 <- g22pre[lower.tri(g22pre,diag=TRUE)])
g12pre <- hapsO.covs12 - diag(4)*(k2.matrix/2)[1,2]
(g12 <- g12pre[lower.tri(g12pre,diag=TRUE)])
(covG <- matrix(c(mean(g11),mean(g12),mean(g12),mean(g22)),byrow=T,nrow=2))

(hapsL.covs <- cov(PHEN[,c('TPL1','TML1','NTPL1','NTML1','TPL2','TML2','NTPL2','NTML2')]))
(hapsL.covs1 <- hapsL.covs[c(1:4),c(1:4)])
(hapsL.covs2 <- hapsL.covs[c(5:8),c(5:8)])
(hapsL.covs12 <- hapsL.covs[5:8,1:4])
h11pre <- hapsL.covs1 - diag(4)*(k2.matrix/2)[1,1]
(h11 <- h11pre[lower.tri(h11pre,diag=TRUE)])
h22pre <- hapsL.covs2 - diag(4)*(k2.matrix/2)[2,2]
(h22 <- h22pre[lower.tri(h22pre,diag=TRUE)])
h12pre <- hapsL.covs12 - diag(4)*(k2.matrix/2)[1,2]
(h12 <- h12pre[lower.tri(h12pre,diag=TRUE)])
(covH <- matrix(c(mean(h11),mean(h12),mean(h12),mean(h22)),byrow=T,nrow=2))

(hapsOL1.covs <- cov(PHEN[,c('TPO1','TMO1','NTPO1','NTMO1','TPL1','TML1','NTPL1','NTML1')]))
(i11 <- as.vector(hapsOL1.covs[5:8,1:4]))
(hapsOL2.covs <- cov(PHEN[,c('TPO2','TMO2','NTPO2','NTMO2','TPL2','TML2','NTPL2','NTML2')]))
(i22 <- as.vector(hapsOL2.covs[5:8,1:4]))
(hapsOL12.covs <- cov(PHEN[,c('TPO1','TMO1','NTPO1','NTMO1','TPL1','TML1','NTPL1','NTML1','TPO2','TMO2','NTPO2','NTMO2','TPL2','TML2','NTPL2','NTML2')]))
(i12 <- as.vector(hapsOL12.covs[1:4,13:16])) #Note: I believe i12 != i21, so am keeping the separate to check. Note that i is defined as [N]T -> L[N]T, thus starting at O -> L, meaning O defines rows (is on y-axis) and L defines columns (is on x-axis). Therefore, i12 is cov(O1,L2) 
(i21 <- as.vector(hapsOL12.covs[9:12,5:8])) #whereas i21 is cov(O2,L1)
(covI <- matrix(c(mean(i11),mean(i12),mean(i21),mean(i22)),byrow=T,nrow=2))

(covAF <- cov(PHEN[,c("F1","F2","AO1","AO2","AL1","AL2")]))
w <- covAF[1:2,3:4] #careful. w = F -> NT. F needs to be on y-axis and AO on x-axis
q <- covAF[1:2,5:6] #same issue. q = F -> LNT.
covF <- covAF[1:2,1:2] #here, it doesn't matter bc cov(F) is symm

(covE <- cov(PHEN[,c("E1","E2")])) #covE is symm

(covY <- cov(PHEN[,c("Y1P","Y2P","Y1M","Y2M","Y1","Y2")]))

(omega.p.T <- cov(PHEN[,c("Y1P","Y2P","TPO1","TPO2")]))
(omega.m.T <- cov(PHEN[,c("Y1M","Y2M","TMO1","TMO2")]))
(omega.T <- (omega.p.T[1:2,3:4] + omega.m.T[1:2,3:4])*.5)
(omega.p.NT <- cov(PHEN[,c("Y1P","Y2P","NTPO1","NTPO2")]))
(omega.m.NT <- cov(PHEN[,c("Y1M","Y2M","NTMO1","NTMO2")]))
(omega.NT <- (omega.p.NT[1:2,3:4] + omega.m.NT[1:2,3:4])*.5)
#(omega <- (omega.T + t(omega.NT))*.5) #the t(omega.NT) made no sense - not sure why I did that
(omega <- (omega.T + omega.NT)*.5) #the t(omega.NT) made no sense - not sure why I did that

(thetaNTp <- cov(PHEN[,c("Y1","Y2","NTPO1","NTPO2")]))
(thetaNTm <- cov(PHEN[,c("Y1","Y2","NTMO1","NTMO2")]))
(thetaNT <- (thetaNTp[1:2,3:4] + thetaNTm[1:2,3:4]))

(thetaTp <- cov(PHEN[,c("Y1","Y2","TPO1","TPO2")]))
(thetaTm <- cov(PHEN[,c("Y1","Y2","TMO1","TMO2")]))
(thetaT <- (thetaTp[1:2,3:4] + thetaTm[1:2,3:4]))

(gamma.p.T <- cov(PHEN[,c("Y1P","Y2P","TPL1","TPL2")]))
(gamma.m.T <- cov(PHEN[,c("Y1M","Y2M","TML1","TML2")]))
(gamma.T <- (gamma.p.T[1:2,3:4] + gamma.m.T[1:2,3:4])*.5)
(gamma.p.NT <- cov(PHEN[,c("Y1P","Y2P","NTPL1","NTPL2")]))
(gamma.m.NT <- cov(PHEN[,c("Y1M","Y2M","NTML1","NTML2")]))
(gamma.NT <- (gamma.p.NT[1:2,3:4] + gamma.m.NT[1:2,3:4])*.5)
#(gamma <- (gamma.T + t(gamma.NT))*.5) #again, can't figure out why I did this
(gamma <- (gamma.T + gamma.NT)*.5)


#SUMMARY.RES WITH NAMES OF THE COLUMNS NEXT TO THEM
SUMMARY.RES[[CUR.GEN+1]] <- list(GEN=CUR.GEN,
                          NUM.CVs=NUM.CVs,
                          MATE.COR=cor.list[[CUR.GEN+1]],
                          POPSIZE=POP.SIZE,
                          VAO=var(PHEN[,c('AOy1','AOy2')]), #new0 - 
                          VAL=var(PHEN[,c('ALy1','ALy2')]), #new1  - 
                          VF=var(PHEN[,c('F1','F2')]),
                          VE=var(PHEN[,c('E1','E2')]),
                          VP=var(PHEN[,c('Y1','Y2')]),
                          h2=(var(PHEN[,c('AOy1','AOy2')]) + var(PHEN[,c('ALy1','ALy2')])) /var(PHEN[,c('Y1','Y2')]), #new0 - 
                          h2.obs=(var(PHEN[,c('AOy1','AOy2')])) /var(PHEN[,c('Y1','Y2')]), #new1 - 
                          h2.lat=(var(PHEN[,c('ALy1','ALy2')])) /var(PHEN[,c('Y1','Y2')]), #new1 - 
                          
                          covY=covY,
                          covG=covG,
                          covH=covH,
                          covI=covI,
                          w=w,
                          q=q,
                          covF=covF,
                          covE=covE,
                          
                          hapsO.covs= hapsO.covs,
                          hapsL.covs= hapsL.covs,
                          
                          omega=omega,
                          gamma=gamma,
                          thetaNT=thetaNT,
                          thetaT=thetaT
                          )

if (SAVE.EACH.GEN){
HISTORY$MATES[[CUR.GEN+1]] <- MATES
HISTORY$PHEN[[CUR.GEN+1]] <- PHEN
HISTORY$XO[[CUR.GEN+1]] <- XO
HISTORY$XL[[CUR.GEN+1]] <- XL
}

cat("Generation",CUR.GEN,"done\n")

if (SAVE.COVS){
    COVS[[CUR.GEN+1]] <- round(covs,3)}

} #End for loop
###################


###################
#A7 Return variables of interest
return(list(SUMMARY.RES=SUMMARY.RES,XO=XO,XL=XL,PHEN=PHEN,HISTORY=HISTORY,COVARIANCES=COVS))

} #End function
###################
#########################################################




##############################################################################
#END FUNCTIONS
##############################################################################














#Checks
#round(cov(PHEN[,c("TPO1","TPL1","TPO2","TPL2")]),3)

#f2 <- fmat
#f2[1,2] <- 0
#(omega <- matrix(c(.4,.05,.05,.1),nrow=2,byrow=T))
#f2 %*% omega














