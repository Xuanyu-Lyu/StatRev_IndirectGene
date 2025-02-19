
#make reltype = 6+number.generations-currun OR just the specific number
#(PAR$id.vector[number.generations+8]+1):PAR$id.vector[number.generations+9]


#for debugging

#phenotype <- make.phenotype(males.effects.mate, females.effects.mate, patgenes, matgenes, TEMP$females.married, TEMP$reltpe,
#                            TEMP$children.vector, corrections=list(), founder.pop=TRUE, PARAM=PAR, A.effects.key=a.effects.key,
#                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix)

#CHANGED - comment this after
#females.married <- TEMP$females.married
#reltpe <- TEMP$reltpe
#children.vector <- TEMP$children.vector
#corrections=list() #OR corrections=correct
#corrections=correct
#founder.pop=TRUE  #OR FALSE! Look to see
#founder.pop=FALSE
#PARAM=PAR
#A.effects.key=a.effects.key
#D.effects.key=d.effects.key
#VARIANCE=VAR
#BETA.matrix=beta.matrix
#center=TRUE

#save(males.effects.mate,females.effects.mate,patgenes,matgenes,TEMP$size.mate,PAR$num.par,a.effects.key,aa.effects.key,
#                           d.effects.key,reltpe,PAR$id.vector,PAR$genenames.row,PAR$parnames,a.mean.correction,
#                           a.sd.correction,aa.mean.correction,aa.sd.correction,d.mean.correction,d.sd.correction,
#                           f.mean.correction,f.sd.correction,s.mean.correction,
#                           s.sd.correction,e.mean.correction,e.sd.correction,age.mean.correction,
#                           age.sd.correction,slope.mean.correction,slope.sd.correction,vt.model,sibling.age.cor,
#                           females.married,children.vector,mvrnormal,effects.mate,beta.matrix,PAR$am.multiplier,file="toMPFunc.RData")

#load("toMPFunc.RData")



make.phenotype <- function(males.effects.mate, females.effects.mate, patgenes, matgenes, females.married, reltpe, children.vector,
                           corrections, twin.age=1, founder.pop=FALSE, center=TRUE,PARAM, A.effects.key, D.effects.key, VARIANCE,
                           BETA.matrix){




############################################################################
#Create GENERATION CUR PHENOTYPES

#make size.gen
size.gen <- number.children <- dim(patgenes)[2]
number.genes <- dim(patgenes)[1]

#create ID's for each person & relatives two generations back
Subj.ID <- sample((PARAM$id.vector[reltpe+1]+1):PARAM$id.vector[reltpe+2],size=size.gen,replace=FALSE)

#Name the rows & columns of genes.mate
genenames.col <- paste(rep(paste("subj",reltpe,".",sep=""),size.gen),Subj.ID,sep="")

#Stack the genes - one complete genotype per column
genes.stack <- rbind(patgenes,matgenes)
dimnames(genes.stack) <- list(PARAM$genenames.row,genenames.col)

#Make new effects.cur matrix
effects.cur <- matrix(0,nrow=PARAM$num.par,ncol=size.gen,dimnames=list(PARAM$parnames,NULL)) 
effects.cur["female",] <- (runif(size.gen)<.5)*1
if (reltpe==1) effects.cur["female",] <- rep(sample(c(0,1),size=size.gen/2,replace=TRUE),each=2)   #MZ twins are the same sex

Father.ID <- males.effects.mate["Subj.ID",]
Mother.ID <- females.effects.mate["Subj.ID",]
Fathers.Father.ID <- males.effects.mate["Father.ID",]
Fathers.Mother.ID <- males.effects.mate["Mother.ID",]
Mothers.Father.ID <- females.effects.mate["Father.ID",]
Mothers.Mother.ID <- females.effects.mate["Mother.ID",]
Spouse.ID <- rep(0,size.gen)
Relative.type <- rep(reltpe,size.gen) #7=twin/spouse parents; 8 on up are grandparents, greatgrandparents, etc...

effects.cur <- rbind(Subj.ID, Father.ID, Mother.ID, Fathers.Father.ID,Fathers.Mother.ID,
                     Mothers.Father.ID,Mothers.Mother.ID,Spouse.ID,Relative.type,genes.stack,effects.cur)
dimnames(effects.cur)[[2]] <- genenames.col
################################



#General note: the variable scores (A, F, U, etc) in effects.cur and (fe)males.effects.mate are NOT the contributions to the phenotype - they typically have a mean=0 and var=1. They are instead just the scores on the latent variables, and their effects on the phenotype have to be multiplied by beta.


### ADDITIVE GENETIC EFFECTS 
#Note: we need to correct the mean and SD of genetic effects with the *same* correction factor from generation 0, making
#all variance components relative to the mean=0 and var=1 in generation 1; if we standardized each generation, this would
#be wrong because the mean=0 and var=1 always, by definition, and thereby be disallowed to change
a.indx <- c(genes.stack)
a.matrix.prelim <- matrix(A.effects.key[a.indx,1],ncol=ncol(genes.stack))
a.prelim <- colSums(a.matrix.prelim) #NOTE - the var(a.prelim) is ~ the number of SNPs, so we must adjust it (next 2 lines) to have var=1
#var(a.prelim)
if (founder.pop){corrections$a.mean <- mean(a.prelim);corrections$a.sd <- 1/sd(a.prelim)}
effects.cur["A",]<- (a.prelim-corrections$a.mean)*corrections$a.sd
################################




### ADDITIVE BY ADDITIVE EPISTATIC GENETIC EFFECTS 
aa.pp <-  genes.stack[1:number.genes,]*rbind(genes.stack[2:number.genes,],genes.stack[1,])
aa.mm <-  genes.stack[(number.genes+1):(number.genes*2),]*rbind(genes.stack[(number.genes+2):(number.genes*2),],genes.stack[(number.genes+1),])
aa.pm <- genes.stack[1:number.genes,]*rbind(genes.stack[(number.genes+2):(number.genes*2),],genes.stack[(number.genes+1),])
aa.mp <- rbind(genes.stack[2:number.genes,],genes.stack[1,])*genes.stack[(number.genes+1):(number.genes*2),]
aa.fin <- rbind(aa.pp,aa.mm,aa.pm,aa.mp)

aa.indx <- c(aa.fin %% 1e4)+1   #the "+1" is necessary to avoid zeros in the index
aa.matrix.prelim <- matrix(A.effects.key[aa.indx],nrow=number.genes*4)
aa.prelim <- colSums(aa.matrix.prelim)
if (founder.pop){corrections$aa.mean <- mean(aa.prelim);corrections$aa.sd <- 1/sd(aa.prelim)}
effects.cur["AA",] <- (aa.prelim-corrections$aa.mean)*corrections$aa.sd 
################################




### DOMINANCE GENETIC EFFECTS 
d.indx <- c((patgenes*matgenes) %% 1e4)+1  
d.matrix.prelim <- matrix(D.effects.key[d.indx],ncol=ncol(genes.stack))
d.prelim <- colSums(d.matrix.prelim)
if (founder.pop){corrections$d.mean <- mean(d.prelim);corrections$d.sd <- 1/sd(d.prelim)}
effects.cur["D",] <- (d.prelim-corrections$d.mean)*corrections$d.sd 
################################




### FAMILIAL ENVIRONMENTAL EFFECTS 
#Vertical Transmission;  transmission of parenting phenotypes to offspring;
#Note - unlike any of the other variance components, V(F) < 1 (Nov 2023) because it is the actual variance of F contributing to the trait and, alone among the effects, has a beta=1. 
#OUTDATED (OLD): We force the var(F) to equal its expectation (F) in gen0. Changes in V(F) thereafter are relative to this original variance
#NOTE: I am not correcting the V(F) in gen0 right now - it is whatever it is
f.prelim <- males.effects.mate["parenting.phenotype",]*VARIANCE$pvt + females.effects.mate["parenting.phenotype",]*VARIANCE$mvt
var(f.prelim) # - should be close to what the user specified as "F" in GE-75.R
#if (founder.pop){corrections$f.mean <- mean(f.prelim); corrections$f.sd <- sqrt(VARIANCE$F/var(f.prelim))} #Here, I changed it so that we are not scaling VF by the mean and sd of VF in gen0
corrections$f.mean <- 0; corrections$f.sd <- 1 #the corrections$f.mean is no longer used anywhere as of March 21, 2024
#effects.cur["F",] <- (f.prelim-corrections$f.mean)*corrections$f.sd  #This commented out on March 21,2024
effects.cur["F",] <- (f.prelim-mean(f.prelim))*corrections$f.sd  #this ensure that VF always has mean 0 each generation; CHANGED on March 21, 2024
var(effects.cur["F",])
################################




### SIBLING ENVIRONMENT
#NOTE: "S" is created from the males.effects.mate["sib.env",] so that siblings always share the same S; could have just as easily used females.effects.mate["sib",]
#(correction factor not really needed b/c these are rnorm() each gen, but for consistency, they're used here too)
s.prelim <- sqrt(1-VARIANCE$A.S.cor^2)*males.effects.mate["sib.env",] + VARIANCE$A.S.cor*effects.cur["A",]
#Note that the above is just 1*males.effects.mate["sib.env",] if A.S.cor==0
if (founder.pop){corrections$s.mean <- mean(s.prelim);corrections$s.sd <- 1/sd(s.prelim)}
effects.cur["S",] <- (s.prelim-corrections$s.mean)*corrections$s.sd 
################################




### CHILDREN'S SIBLING ENVIRONMENT 
#(correction factor not really needed b/c these are rnorm() each gen, but for consistency, they're used here too)
#this is a factor that allows the children of this individual to share a sib env. It does not affect the parenting phenotype
#I think this is here so that the kids of *this parent* have the same S effect (Nov 2023)
sib.prelim <- rnorm(size.gen)
if (founder.pop){corrections$sib.mean <- mean(sib.prelim); corrections$sib.sd <- 1/sd(sib.prelim)}
effects.cur["sib.env",] <- (sib.prelim-corrections$sib.mean)*corrections$sib.sd
################################




### UNIQUE ENVIRONMENTAL EFFECTS
#(correction factor not really needed b/c these are rnorm() each gen, but for consistency, they're used here too)
u.prelim <- sqrt(1-VARIANCE$A.U.cor^2)*rnorm(size.gen) + VARIANCE$A.U.cor*effects.cur["A",]
if (founder.pop){corrections$u.mean <- mean(u.prelim);corrections$u.sd <- 1/sd(u.prelim)}
effects.cur["U",] <- (u.prelim-corrections$u.mean)*corrections$u.sd
################################




### MZ & TW ENVIRONMENTAL EFFECTS
#(correction factor not really needed b/c these are rnorm() each gen, but for consistency, they're used here too)
if (reltpe==1){
  mz.prelim <- rep(rnorm(length(children.vector)),times=children.vector)
  tw.prelim <- rep(rnorm(length(children.vector)),times=children.vector)}
if (reltpe==2){
  mz.prelim <- rnorm(size.gen)
  tw.prelim <- rep(rnorm(length(children.vector)),times=children.vector)}
if (reltpe>2){
  mz.prelim <- rnorm(size.gen)
  tw.prelim <- rnorm(size.gen)}
if (founder.pop){corrections$mz.mean <- mean(mz.prelim);corrections$mz.sd <- 1/sd(mz.prelim)
               corrections$tw.mean <- mean(tw.prelim);corrections$tw.sd <- 1/sd(tw.prelim)}
effects.cur["MZ",] <- (mz.prelim-corrections$mz.mean)*corrections$mz.sd
effects.cur["TW",] <- (tw.prelim-corrections$tw.mean)*corrections$tw.sd

#the empirical value of MZ or TW (unstandardized); we'll use this later for summarizing (below, Finding Variance Components)
if (reltpe>2){empirical.tw <- 0}
if (reltpe==2){empirical.tw <- var(effects.cur["TW",])} 
if (reltpe==1){empirical.tw <- var(effects.cur["MZ",])} 
################################




### SEX EFFECTS
sex.prelim <- effects.cur["female",]
effects.cur["SEX",] <- (sex.prelim-.5)/.5        #always use population mean (.5) and sd of sex (.5) - i.e., sqrt(var(x)) = sqrt(.5^2) = .5
if (center==FALSE) effects.cur["SEX",] <- effects.cur["SEX",]+.5    #this uncenters sex (range -.5 to 1.5), making possible scalar interactions with no main effects
################################




### AGE EFFECTS 
#make siblings ages s.t. they are clustered within families according to PARAM$sibling.age.cor
suppressWarnings(varcovar <- matrix(rep(c(1,rep(PARAM$sibling.age.cor,10)),9),nrow=10))
norm.age <- mvrnormal(females.married,mu=rep(0,10),Sigma=varcovar,empirical=FALSE) #normally dist age, s.t. sibs have intraclass correlation
cum.prob.age.spacing <- pnorm(norm.age)          #convert normally dist. ages to cumulative probability; i.e., 0 to 1 uniform
age.new <- qunif(cum.prob.age.spacing,PARAM$range.pop.age[1],PARAM$range.pop.age[2])  #find unif value at these probabilities

#for ancestors, parents, & children of twins
if (reltpe>3){
row.selector <- rep(1:females.married,time=children.vector)
age2 <- age.new[row.selector,]                   #matrix with #rows=number.children, repeated within families
age3 <- as.vector(t(age2))                       #vectorize age
adder1 <- rep(1:10,ceiling(number.children/10))[1:number.children] 
adder2 <- 0:(number.children-1)*10
age.fin <- age3[(adder1+adder2)]                 #choose ages
effects.cur["cur.age",] <- age.fin
effects.cur["age.at.mar",] <-  runif(size.gen,18,36)-(effects.cur["female",]*3) #makes females on avg 3 years younger @ marriage
if (founder.pop){corrections$age.mean <- mean(effects.cur["cur.age",]);corrections$age.sd <- 1/sd(effects.cur["cur.age",])}
effects.cur["AGE",] <- (effects.cur["cur.age",]-corrections$age.mean)*corrections$age.sd
effects.cur["AGE.M",] <- (effects.cur["age.at.mar",]-corrections$age.mean)*corrections$age.sd}

#for sibs
if (reltpe==3){
sibplustwin <- rep(1,females.married)+children.vector
num.tot <- sum(sibplustwin)
row.selector <- rep(1:females.married,time=sibplustwin)
age2 <- age.new[row.selector,]                   #matrix with #rows=number of sibs + twins, repeated within families
age3 <- as.vector(t(age2))                       #vectorize age
adder1 <- rep(1:10,ceiling(num.tot/10))[1:num.tot] 
adder2 <- 0:(num.tot-1)*10
age.fin <- age3[(adder1+adder2)]                 #choose ages for sibs + twins (to be carried forward)
sib.chooser <- rep(which(children.vector != 0),times=children.vector[children.vector != 0])+0:(number.children-1)
sib.age <- age.fin[sib.chooser]
twin.age <- age.fin[-sib.chooser] #twin age gets carried forward, to be input into make.phenotype for MZ and DZ twins
effects.cur["age.at.mar",] <- runif(size.gen,18,36)-(effects.cur["female",]*3) #makes females on avg 3 years younger @ marriage
effects.cur["cur.age",] <- sib.age
effects.cur["AGE",] <- (effects.cur["cur.age",]-corrections$age.mean)*corrections$age.sd
effects.cur["AGE.M",] <- (effects.cur["age.at.mar",]-corrections$age.mean)*corrections$age.sd}

#for dz & mz
if (reltpe<3){
tw.age <- rep(twin.age[which(children.vector == 2)],each=2)
effects.cur["age.at.mar",] <- runif(size.gen,18,36)-(effects.cur["female",]*3) #makes females on avg 3 years younger @ marriage
effects.cur["cur.age",] <-  tw.age
effects.cur["AGE",] <- (effects.cur["cur.age",]-corrections$age.mean)*corrections$age.sd
effects.cur["AGE.M",] <- (effects.cur["age.at.mar",]-corrections$age.mean)*corrections$age.sd}
################################





### SEX by ADDITIVE GENETIC EFFECTS 
sex.slope.prelim <- colSums(matrix(A.effects.key[a.indx,2],ncol=ncol(genes.stack)))   #individual gene effects for A,sex slopes

#correction factors to be used hereafter
if (founder.pop){corrections$sex.slope.mean <- mean(sex.slope.prelim);corrections$sex.slope.sd <- 1/sd(sex.slope.prelim)}
effects.cur["A.slopes.sex",] <- (sex.slope.prelim-corrections$sex.slope.mean)*corrections$sex.slope.sd
effects.cur["A.by.SEX",] <- effects.cur["A.slopes.sex",]*effects.cur["SEX",]
################################





### AGE by ADDITIVE GENETIC EFFECTS 
age.slope.prelim <- colSums(matrix(A.effects.key[a.indx,3],ncol=ncol(genes.stack)))   #individual gene effects for A,age slopes

#correction factors to be used hereafter
if (founder.pop){corrections$age.slope.mean <- mean(age.slope.prelim);corrections$age.slope.sd <- 1/sd(age.slope.prelim)}
effects.cur["A.slopes.age",] <- (age.slope.prelim-corrections$age.slope.mean)*corrections$age.slope.sd
effects.cur["A.by.AGE",] <- effects.cur["A.slopes.age",]*effects.cur["AGE",]
effects.cur["A.by.AGE.M",] <- effects.cur["A.slopes.age",]*effects.cur["AGE.M",]

#NOTE: AGE.M and A.by.AGE.M are the age and A*age interaction effects at the age of marriage. These will be used to calculate the mating phenotype and the parenting phenotype (ostensibly closer in time to the marriatal age than the age when the person's data is collected). (Nov 2023)

################################





### S by ADDITIVE GENETIC EFFECTS 
S.slope.prelim <- colSums(matrix(A.effects.key[a.indx,4],ncol=ncol(genes.stack)))   #individual gene effects for A,S slopes

#correction factors to be used hereafter
if (founder.pop){corrections$S.slope.mean <- mean(S.slope.prelim);corrections$S.slope.sd <- 1/sd(S.slope.prelim)}
effects.cur["A.slopes.S",] <- (S.slope.prelim-corrections$S.slope.mean)*corrections$S.slope.sd
effects.cur["A.by.S",] <- effects.cur["A.slopes.S",]*effects.cur["S",]
################################





### U by ADDITIVE GENETIC EFFECTS 
U.slope.prelim <- colSums(matrix(A.effects.key[a.indx,5],ncol=ncol(genes.stack)))   #individual gene effects for A,U slopes

#correction factors to be used hereafter
if (founder.pop){corrections$U.slope.mean <- mean(U.slope.prelim);corrections$U.slope.sd <- 1/sd(U.slope.prelim)}
effects.cur["A.slopes.U",] <- (U.slope.prelim-corrections$U.slope.mean)*corrections$U.slope.sd
effects.cur["A.by.U",] <- effects.cur["A.slopes.U",]*u.prelim
################################





### CREATE 3 PHENOTYPES (actual phenotype, mating phenotype, and parenting phenotype)

#NOTE: This is a reversion (as of March 21, 2024) of the old way I used to create phenotypes. It does not scale the phenotypes relative to what they were in gen0. I reverted to this because the scaling was forcing the variance of the parental (or mating) phenotype to be what was specified in "F" in the GE-75.R script, instead of allowing it to drift to lower or higher values, as it naturally should.
effects.cur["cur.phenotype",] <-    BETA.matrix %*%  effects.cur[PARAM$varnames,] #Old way
var(effects.cur["cur.phenotype",]) #should be the sum of all variances specified by user that make up the phenotype

effects.cur["mating.phenotype",]  <- (BETA.matrix*PARAM$am.multiplier) %*% (effects.cur[PARAM$varnames.M,]) #Old way
var(effects.cur["mating.phenotype",]) #should be the sum of all variances specified by user that make up the mating phenotype

effects.cur["parenting.phenotype",] <- (BETA.matrix*PARAM$vt.model) %*% (effects.cur[PARAM$varnames.M,])
var(effects.cur["parenting.phenotype",]) #should be the sum of all variances specified by user that make up the mating phenotype



#NOTE: all the below is now commented out. To go back to the newer way, uncomment these and comment the ones above.

#Note - The way you standardize assumes all effects are uncorrelated. That is fine for effects that become correlated naturally due to AM or VT because we are concerned with the base population. However, the 4 correlations between intercepts and slopes in GE-75.R will NOT be accounted for and they should, so the way I'm doing it now (as of Nov 2023) will lead to variances != what you think they should be at the start of the simulation when there are these correlated effects. Essentially, this simulation assumes ALL effects are uncorrelated in the base population.

#Note 2 - this also means that the phenotypes will NOT exactly equal the sum of all the effects - the phenotypes will be correlated at +1.0 with them, but they'll differ slightly, and their mean and var will differ slightly


#Note - We want the phenotype to start with mean=0 & var=its expected variance=par0 in gen0; we use this as the baseline thereafter
#pheno.prelim <- as.vector(BETA.matrix %*% (effects.cur[PARAM$varnames,]))
#var(pheno.prelim) #should be close to sum of all variances specified by user that make up the parenting phenotype
#(phen0 <- sum((BETA.matrix^2)*c(rep(1,3),VARIANCE$F,rep(1,10)))) #expected variance of phenotype assuming all effects are uncorrelated
#if (founder.pop){corrections$pheno.mean <- mean(pheno.prelim); corrections$pheno.sd <- sqrt(phen0/var(pheno.prelim))}
#effects.cur["cur.phenotype",] <- (pheno.prelim-corrections$pheno.mean)*corrections$pheno.sd  
#var(effects.cur["cur.phenotype",]) #should be the sum of all variances specified by user that make up the phenotype


#Note - We want the mating phenotype to start with mean=0 & var=its expected variance=par0 in gen0; we use this as the baseline thereafter. We use varnames.M because these are the values of Age & AxAge at the time the two spouses married
#(am.paths <- (BETA.matrix*PARAM$am.multiplier))
#mating.prelim <- as.vector(am.paths %*% (effects.cur[PARAM$varnames.M,]))
#var(mating.prelim) #should be close to sum of all variances specified by user that make up the parenting phenotype
#(mating0 <- sum((am.paths^2)*c(rep(1,3),VARIANCE$F,rep(1,10)))) #expected variance of mating phenotype assuming all effects are uncorrelated
#if (founder.pop){corrections$mating.mean <- mean(mating.prelim); corrections$mating.sd <- sqrt(mating0/var(mating.prelim))}
#effects.cur["mating.phenotype",] <- (mating.prelim-corrections$mating.mean)*corrections$mating.sd  
#var(effects.cur["mating.phenotype",]) #should be the sum of all variances specified by user that make up the mating phenotype


#NOTE - the problem with VF not halving when pvt = mvt = .5 is probably here; does this keep the parenting phenotype the same each gen?

#Note - We want parenting phenotype to start with mean=0 & var=its expected variance=par0 in gen0; we use this mean & var of the parenting phenotype in gen0 as the baseline thereafter. We use varnames.M because these are the values of Age & AxAge at the time the two spouses married, which would be closer to when they have offspring and pass down their parental traits
#(par.paths <- (BETA.matrix*PARAM$vt.model))
#parenting.prelim <- as.vector(par.paths %*% (effects.cur[PARAM$varnames.M,]))
#var(parenting.prelim) #should be close to sum of all variances specified by user that make up the parenting phenotype
#(par0 <- sum((par.paths^2)*c(rep(1,3),VARIANCE$F,rep(1,10)))) #expected variance of parenting phenotype assuming all effects are uncorrelated
#if (founder.pop){corrections$parenting.mean <- mean(parenting.prelim); corrections$parenting.sd <- sqrt(par0/var(parenting.prelim))}
#effects.cur["parenting.phenotype",] <- (parenting.prelim-corrections$parenting.mean)*corrections$parenting.sd  
#var(effects.cur["parenting.phenotype",]) #should be the sum of all variances specified by user that make up the mating phenotype

################################



#make a list of what is returned from the function
list(effects.mate=effects.cur,size.mate=ncol(effects.cur),empirical.tw=empirical.tw,twin.age=twin.age,corrections=corrections)

#END FUNCTION
}

############################################################################





