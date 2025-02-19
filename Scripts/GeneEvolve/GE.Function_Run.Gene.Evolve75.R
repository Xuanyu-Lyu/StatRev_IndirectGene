



#Updated March 20, 2024 by mck
#For debugging (should be commented out unless debugging):
#Run GeneEvolve
#ALL <- run.gene.evolve(PAR=PAR1,VAR=VAR1,iter=iteration)
#PAR=PAR1
#VAR=VAR1
#iter=1
#source("http://www.matthewckeller.com/R.Class/KellerScript3.R") 





run.gene.evolve <- function(PAR,VAR,iter){ 



############################################################################
#2.1 CORRECT ERRORS

#give errors or warnings #CHANGE - create new errors to do with VT, esp vt.model
if (PAR$number.genes<2) {stop("You must have at least two genes - even if there are no genetic effects")}
if (PAR$number.genes>50) {stop("Please specify 50 genes or fewer")}

if (PAR$range.dataset.age[1]>=max(PAR$range.twin.age)) {
  stop("Woops: the minimum ages in your dataset (range.dataset.age[1]) is >= the maximum ages of twins you want - You're not going to have any twins!")}
if (PAR$range.dataset.age[2]<=min(PAR$range.twin.age)) {
  stop("Woops: the maximum ages in your dataset (range.dataset.age[2]) is <= the minimum ages of twins you want - You're not going to have any twins!")}

if (length(VAR$latent.AM) != PAR$number.generations) {stop("You must have as many entries in the VAR$latent.AM vector as the number of runs
                                           For constant AM, try 'rep(AM,number.generations)' where AM is your chosen AM coefficient")}
{if (PAR$start.pop.size<500) {warning("The starting population size is less than 500. The script may crash with such small sample sizes
                        under certain circumstances (e.g., there are too few husbands or wives).
                        Don't blame me if it crashes on you!")
                        PAR$startpop.warn <- "Script may crash when n < 500"}
 else PAR$startpop.warn <- ""}

if (odd(length(PAR$range.twin.age))){
  stop("range.twin.age must be a VECTOR of length 2,4,6,8,... (etc).")}

if (max(PAR$range.pop.age) > PAR$range.dataset.age[2]) {
  stop("Your range of ages in your population exceeds your max age (range.dataset.age[2]). Please change one or the other.")}

if (min(PAR$range.pop.age) < PAR$range.dataset.age[1]) {
  stop("Your range of ages in your population exceeds your max.age. Please change one or the other.")}

############################################################################

#22222222222222222222222222222222222222222222222222222222222222222222222222222222222222







#3 CREATING PARAMETERS THAT DO NOT CHANGE AND STARTING THE ITERATION FOR EACH RUN

############################################################################
#make beta coefficients
beta <- list(); TEMP<- list()
beta$A <- sqrt(VAR$A)
beta$AA <- sqrt(VAR$AA)
beta$D <- sqrt(VAR$D)
beta$F <- 1 #CHANGED from sqrt(VAR$F) - this now mimics how VT is treated in the Cascade model
beta$S <- sqrt(VAR$S)
beta$U <- sqrt(VAR$U)
beta$MZ <- sqrt(VAR$MZ)
beta$TW <- sqrt(VAR$TW)
beta$SEX <- sqrt(abs(VAR$SEX))*sign(VAR$SEX)
beta$AGE <- sqrt(abs(VAR$AGE))*sign(VAR$AGE)
beta$A.by.SEX <- sqrt(abs(VAR$A.by.SEX))*sign(VAR$A.by.SEX)                
beta$A.by.AGE <- sqrt(abs(VAR$A.by.AGE))*sign(VAR$A.by.AGE)                
beta$A.by.S <- sqrt(abs(VAR$A.by.S))*sign(VAR$A.by.S)                
beta$A.by.U <- sqrt(abs(VAR$A.by.U))*sign(VAR$A.by.U)                

#create beta matrix for creating phenotype at random age
(beta.matrix <-matrix(c(beta$A,beta$AA,beta$D,beta$F,beta$S,beta$U,beta$MZ,beta$TW,beta$SEX,beta$AGE,beta$A.by.SEX,beta$A.by.AGE,beta$A.by.S,beta$A.by.U),nrow=1))
(var.matrix <- abs(matrix(c(VAR$A,VAR$AA,VAR$D,VAR$F,VAR$S,VAR$U,VAR$MZ,VAR$TW,VAR$SEX,VAR$AGE,VAR$A.by.SEX,VAR$A.by.AGE,VAR$A.by.S,VAR$A.by.U),nrow=1)))
(VAR$TOTAL <- sum(var.matrix))
colnames(var.matrix) <- colnames(beta.matrix) <- PAR$varnames <- c("A","AA","D","F","S","U","MZ","TW","SEX","AGE","A.by.SEX","A.by.AGE","A.by.S","A.by.U")
PAR$number.parameters <- dim(beta.matrix)[2]

#create matrices for creating phenotype at marriage
PAR$varnames.M <- c("A","AA","D","F","S","U","MZ","TW","SEX","AGE.M","A.by.SEX","A.by.AGE.M","A.by.S","A.by.U")

#create multiplier vector that will create the phenotype on which spouses choose each other
#if (PAR$am.model=="I") {PAR$am.multiplier <- rep(1,PAR$number.parameters)}   #phenotypic homogamy
#if (PAR$am.model=="II") {PAR$am.multiplier <- c(0,0,0,1,1,1,1,1,0,0,0,0,0,0)}  #social homogamy - all environmental factors
#if (PAR$am.model=="III") {PAR$am.multiplier <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0)} #convergence - only through unique environmental factors
#if (PAR$am.model=="IV") {PAR$am.multiplier <- c(1,1,1,0,0,0,0,0,0,0,0,0,0,0)} #genetic homogamy - only through genetic factors
PAR$am.multiplier <- PAR$am.model #here, we just make the multipler the same as am.model as a quick fix to allow for any type of AM


#make expected population sizes
PAR$popsize <- vector(length=PAR$number.generations+1)
PAR$popsize[1] <- PAR$start.pop.size
for(i in 1:PAR$number.generations){PAR$popsize[i+1] <- PAR$popsize[i]*PAR$pop.growth[i]}; remove(i)
PAR$popsize <- PAR$popsize[-1]

##create id.vector from which subject numbers for all ancestors & relatives will be drawn
# insures that no-one will ever have the same ID number
PAR$id.vector <- c(1,(1:(PAR$number.generations+9))*(max(PAR$popsize)*10))
         
#make other parameters
PAR$parnames <- c("marriageable","female","age.at.mar","cur.age","A","AA","D","F","S","U","MZ","TW","SEX","AGE","A.by.SEX","A.by.AGE",
              "A.by.S","A.by.U","A.slopes.sex","A.slopes.age","A.slopes.S","A.slopes.U","cur.phenotype","mating.phenotype","parenting.phenotype",
               "AGE.M","A.by.AGE.M","sib.env")  #names for parameter matrix
PAR$num.par <- length(PAR$parnames)                                                  #should be equal to all pars from above
PAR$frontnames <- c("Subj.ID","Father.ID","Mother.ID","Fathers.Father.ID","Fathers.Mother.ID","Mothers.Father.ID","Mothers.Mother.ID","Spouse.ID","Relative.type")

#create parameters matrix
TEMP$row.par <- c(PAR$varnames,"AS.cor","AU.cor","sib.cor","R.slope.sex","R.slope.age","R.slope.S","R.slope.U","num.genes","expected.pop","expec.cor.sps") #CHANGED - REMOVED am.mod
TEMP$col.par <- paste(rep("reltp",PAR$number.generations-1),(PAR$number.generations+6):7,sep="")
parameters <- matrix(0,nrow=length(TEMP$row.par),ncol=length(TEMP$col.par))
rownames(parameters) <- TEMP$row.par
colnames(parameters) <- TEMP$col.par

#if (PAR$am.model=="I"){parameters["expec.cor.sps",] <- VAR$latent.AM} #CHANGED - REMOVED am.mod
#if (PAR$am.model=="II"){parameters["expec.cor.sps",] <- VAR$latent.AM*(VAR$F+VAR$S+VAR$U+VAR$TW+VAR$MZ)/VAR$TOTAL}#CHANGED - REMOVED am.mod
#if (PAR$am.model=="III"){parameters["expec.cor.sps",] <- VAR$latent.AM*(VAR$U)/VAR$TOTAL}#CHANGED - REMOVED am.mod
#if (PAR$am.model=="IV"){parameters["expec.cor.sps",] <- VAR$latent.AM*(VAR$A+VAR$D+VAR$AA)/VAR$TOTAL}#CHANGED - REMOVED am.mod
am.beta.matrix <-matrix(c(A,AA,D,F,S,U,MZ,TW,SEX,AGE,A.by.SEX,A.by.AGE,A.by.S,A.by.U),nrow=1) #CHANGED 4/18 - does it work?
parameters["expec.cor.sps",] <- VAR$latent.AM*sum(am.beta.matrix*am.model)/VAR$TOTAL #CHANGED 4/18 - does it work?


parameters["expected.pop",] <- PAR$popsize
parameters[1:(nrow(parameters)-2),] <- matrix(rep(c(var.matrix,VAR$A.S.cor,VAR$A.U.cor,PAR$sibling.age.cor,VAR$R.Alevel.Aslope.sex,VAR$R.Alevel.Aslope.age,VAR$R.Alevel.Aslope.S,VAR$R.Alevel.Aslope.U,PAR$number.genes),each=PAR$number.generations),nrow=nrow(parameters)-2,ncol=ncol(parameters),byrow=TRUE)

#create list ("keep") of objects which will not be erased after each iteration
PAR$keep <- c(ls(),"keep","run.number")
############################################################################

#33333333333333333333333333333333333333333333333333333333333333333333333333333333333333



 




#4 CREATING THE ANCESTRAL (GENERATION 0) POPULATION

############################################################################
#4.1 Make Genetic Effects
#make genes;  see the function make.gene.matrix; change defaults by changing that function
gene.list <- make.gene.matrix(gene.number=PAR$number.genes,gm=PAR$gene.model,psize=PAR$start.pop.size,allele.num=PAR$number.alleles)

#Name the rows of genes.cur
PAR$genenames.row <- c(paste("Pat.Loc",1:PAR$number.genes,sep=""),paste("Mat.Loc",1:PAR$number.genes,sep=""))

#Stack the genes - one complete genotype per column
patgenes <- gene.list$data.pat
matgenes <- gene.list$data.mat

#Create genetic conversion tables ("effects.key"); this converts the allele name/label to its two effects 
#Each allele has four effects - its level (main effect, e.g., A) and its 3 slopes (changes - i.e. interactions - due to age, S, or U)
#these effects are correlated at r=R.Alevel.Aslope.X where X is age, S, or U;
#we multiply the correlations to ensure that the sigma matrix is positive definite
TEMP$sigmat <- matrix(
     c(1,VAR$R.Alevel.Aslope.sex,VAR$R.Alevel.Aslope.age,VAR$R.Alevel.Aslope.S,VAR$R.Alevel.Aslope.U,
     VAR$R.Alevel.Aslope.sex,1,VAR$R.Alevel.Aslope.sex*VAR$R.Alevel.Aslope.age, VAR$R.Alevel.Aslope.sex*VAR$R.Alevel.Aslope.S, VAR$R.Alevel.Aslope.sex*VAR$R.Alevel.Aslope.U,
     VAR$R.Alevel.Aslope.age, VAR$R.Alevel.Aslope.sex*VAR$R.Alevel.Aslope.age, 1, VAR$R.Alevel.Aslope.age*VAR$R.Alevel.Aslope.S, VAR$R.Alevel.Aslope.age*VAR$R.Alevel.Aslope.U,
     VAR$R.Alevel.Aslope.S, VAR$R.Alevel.Aslope.sex*VAR$R.Alevel.Aslope.S, VAR$R.Alevel.Aslope.age*VAR$R.Alevel.Aslope.S, 1, VAR$R.Alevel.Aslope.S*VAR$R.Alevel.Aslope.U,
     VAR$R.Alevel.Aslope.U, VAR$R.Alevel.Aslope.sex*VAR$R.Alevel.Aslope.U, VAR$R.Alevel.Aslope.age*VAR$R.Alevel.Aslope.U, VAR$R.Alevel.Aslope.S*VAR$R.Alevel.Aslope.U, 1),
                      nrow=5,byrow=TRUE)
a.effects.key <- mvrnormal(10^4, mu=c(0,0,0,0,0),Sigma=TEMP$sigmat, empirical=TRUE)
d.effects.key <- rnorm(10^4)
aa.effects.key <- rnorm(10^4)

if (PAR$save.objects=="no") {remove(gene.list) } #@@@@@CLEANUP
############################################################################






############################################################################
#4.2 CREATE GENERATION CUR PHENOTYPES - GENERATION -1 AND THEN GENERATION 0

#make reltype for ancestors -1, males.effects.mate & females.effects.mate, etc for founder pops
TEMP$reltpe <- PAR$number.generations +7
TEMP$size.mate <- TEMP$females.married <- PAR$start.pop.size
if (odd(TEMP$size.mate)==TRUE){ TEMP$size.mate <- TEMP$size.mate +  sample(c(-1,1),1)} #make number.married always even

#make males & females.effects.mate for gen -1
#This creates RANDOM rnorm() distributed variables for all of the effects for the -1 gen
males.effects.mate <-rbind(matrix(0,nrow=length(PAR$frontnames)+(PAR$number.genes*2),ncol=TEMP$size.mate),matrix(rnorm(length(PAR$parnames)*TEMP$size.mate),nrow=length(PAR$parnames),ncol=TEMP$size.mate))
females.effects.mate <-rbind(matrix(0,nrow=length(PAR$frontnames)+(PAR$number.genes*2),ncol=TEMP$size.mate),matrix(rnorm(length(PAR$parnames)*TEMP$size.mate),nrow=length(PAR$parnames),ncol=TEMP$size.mate))
rownames(females.effects.mate) <- rownames(males.effects.mate) <-c(PAR$frontnames,PAR$genenames.row,PAR$parnames)
TEMP$children.vector <- rep(2,TEMP$size.mate)   #for gen -1, each parent has exactly one kid #CHANGED - is it really 1? or is it 2?

#Change F to be function of cur.phenotype for gen -1
males.effects.mate["F",] <- rnorm(ncol(males.effects.mate),mean=0,sd=sqrt(VAR$F)) #changes F to have correct variance as specified by VF
females.effects.mate["F",] <- rnorm(ncol(females.effects.mate),mean=0,sd=sqrt(VAR$F)) #changes F to have correct variance as specified by VF

#Change parenting.phenotype to be function of cur.phenotype for gen -1
VAR$par.paths <- beta.matrix*PAR$vt.model #creates paths from parental phenotype to latent parenting trait
par0 <- sum((VAR$par.paths^2)*c(rep(1,3),VAR$F,rep(1,10)))
males.effects.mate["parenting.phenotype",] <- males.effects.mate["cur.phenotype",] * sqrt(par0/VAR$TOTAL)
females.effects.mate["parenting.phenotype",] <- females.effects.mate["cur.phenotype",] * sqrt(par0/VAR$TOTAL)

#Change mating.phenotype to be function of cur.phenotype for gen -1
VAR$am.paths <- beta.matrix*PAR$am.multiplier #creates paths from parental phenotype to latent mating trait
mat0 <- sum((VAR$am.paths^2)*c(rep(1,3),VAR$F,rep(1,10)))
males.effects.mate["mating.phenotype",] <- males.effects.mate["cur.phenotype",] * sqrt(mat0/VAR$TOTAL)
females.effects.mate["mating.phenotype",] <- females.effects.mate["cur.phenotype",] * sqrt(mat0/VAR$TOTAL)


#figure out other variables of gen -1
TEMP$Subj.ID <- sample((PAR$id.vector[TEMP$reltpe+2]+1):PAR$id.vector[TEMP$reltpe+3],size=TEMP$size.mate*2,replace=FALSE)
males.effects.mate["Subj.ID",] <- TEMP$Subj.ID[1:TEMP$size.mate]
females.effects.mate["Subj.ID",] <- TEMP$Subj.ID[(TEMP$size.mate+1):(TEMP$size.mate*2)]

if (PAR$save.objects=="yes"){effects.neg1 <- cbind(males.effects.mate,females.effects.mate)}

#make gen 0 phenotypes
phenotype <- make.phenotype(males.effects.mate, females.effects.mate, patgenes, matgenes, TEMP$females.married, TEMP$reltpe,
                            TEMP$children.vector, corrections=list(), founder.pop=TRUE, PARAM=PAR, A.effects.key=a.effects.key,
                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix)

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- phenotype$size.mate
correct <- phenotype$corrections #THIS is where we set the correction factor (mean & sd of each factor) for the rest of the simulation

if (PAR$save.objects == "no") {remove(phenotype)}
############################################################################








############################################################################
#4.3 CREATE MATRIX TO TRACK CHANGES OVER TIME

#make the proto-track.changes matrix
TEMP$track.col <- c(paste(rep("reltp",PAR$number.generations+6),(PAR$number.generations+7):1,sep=""),"reltp.1-7",
                    paste(rep("data.t", length(PAR$range.twin.age)/2),1:(length(PAR$range.twin.age)/2),sep=""))
TEMP$names <- c(rownames(effects.mate)[(14+PAR$number.genes*2):(nrow(effects.mate)-3)])
TEMP$ne <- length(TEMP$names)
TEMP$mean.names <- paste("mean.",TEMP$names,sep="")
TEMP$var.names <- paste("var.",TEMP$names,sep="")
TEMP$cov.names <- c(paste("Cov.",TEMP$names[1],",",TEMP$names[2:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[2],",",TEMP$names[3:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[3],",",TEMP$names[4:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[4],",",TEMP$names[5:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[5],",",TEMP$names[6:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[6],",",TEMP$names[7:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[7],",",TEMP$names[8:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[8],",",TEMP$names[9:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[9],",",TEMP$names[10:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[10],",",TEMP$names[11:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[11],",",TEMP$names[12:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[12],",",TEMP$names[13:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[13],",",TEMP$names[14:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[14],",",TEMP$names[15:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[15],",",TEMP$names[16:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[16],",",TEMP$names[17:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[17],",",TEMP$names[18:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[18],",",TEMP$names[19:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[19],",",TEMP$names[20:TEMP$ne],sep=""),
              paste("Cov.",TEMP$names[20],",",TEMP$names[21],sep=""))

TEMP$track.row <-c(TEMP$mean.names,TEMP$var.names,TEMP$cov.names,"cor.spouses","popsize")

track.changes <- matrix(-999,nrow=length(TEMP$track.row),ncol=length(TEMP$track.col))
dimnames(track.changes) = list(TEMP$track.row,TEMP$track.col)
TEMP$track.beta <- c(beta.matrix,rep(1,7))

#create the track.changes matrix
track.changes <- change.tracker(run=0,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=sum(var.matrix*PAR$am.multiplier)*VAR$latent.AM[1])


#create a matrix to hold observed intraclass coefficients
if (PAR$sib.cor=="yes"){
intra <- matrix(-999,nrow=PAR$num.par-4,ncol=PAR$number.generations+6)
TEMP$track.col <- paste(rep("reltp",PAR$number.generations+5),(PAR$number.generations+6):1,sep="")
dimnames(intra) <- list(PAR$parnames[5:length(PAR$parnames)],TEMP$track.col)}


#write these files out to the working directory if PAR$save.gen.data = "yes"
if (PAR$save.gen.data=="yes"){
  write.table(t(effects.mate),file='Ancestral.Pop.0', append=FALSE,row.names=FALSE)}


#save this file if PAR$save.objects=="yes"
if (PAR$save.objects=="yes"){effects.0 <- effects.mate}
############################################################################

#44444444444444444444444444444444444444444444444444444444444444444444444444444444444444










#5 CREATING THE ANCESTRAL POPULATIONS THROUGH ITERATION (GENERATIONS 1 THROUGH PAR$NUMBER.GENERATIONS) 

    #########################################################
    #########################################################
    #5.05 START ANCESTRAL (INNER) LOOP
    for(currun in 1:PAR$number.generations){
    #########################################################
    #########################################################


       
     
############################################################################
#5.1 ASSORTATIVE MATING - ANCESTORS
#one thing this does not account for is AM on multiple traits

#make effects.mate matrix by lopping off the previous generation from effects.mate (at the beginning of effects.mate)
TEMP$header.cols <- nrow(effects.mate)

#assortative mating
amcur <- assort.mate(effects.mate,VAR$latent.AM,currun,stop.inbreeding=TRUE,mu.constant=PAR$mu.constant)
TEMP$females.married <- amcur$females.married
TEMP$males.married <- amcur$males.married
males.effects.mate <- amcur$males.effects.mate
females.effects.mate <- amcur$females.effects.mate
TEMP$size.mate <- amcur$size.mate
effects.mate <- amcur$effects.mate
TEMP$cor.spouses <- amcur$cor.spouses

#write these files out to the working directory if PAR$save.gen.data == "yes"
if (PAR$save.gen.data=="yes"){
write.table(t(males.effects.mate),file=paste("Ancestral.BreedingMales.", currun, sep=""),append=FALSE,row.names=FALSE)
write.table(t(females.effects.mate),file=paste("Ancestral.BreedingFemales.", currun, sep=""),append=FALSE,row.names=FALSE)}

if (PAR$save.objects=="no") {remove(amcur)}

if (PAR$save.objects=="yes"){
assign(paste("males.effects.mate.",currun,sep=""),males.effects.mate)
assign(paste("females.effects.mate.",currun,sep=""),females.effects.mate)}
############################################################################




############################################################################
#5.2 CREATING CHILDREN - ANCESTORS

#number of children needed to have population growth as expected this generation
TEMP$nec.child.per.marriage <- PAR$popsize[currun]/TEMP$females.married

#vector of the actual # children for each marriage; a non-iid poisson process
TEMP$children.vector <- rpois(TEMP$females.married,lambda=TEMP$nec.child.per.marriage)
TEMP$children.vector[TEMP$children.vector>9] <- 10          #can never be more than 10 children per female (for computational reasons)

#reproduction - over-writing males.effects.mate to be those of the children rather than parents above
children.effects <- reproduce(TEMP$children.vector,TEMP$females.married,males.effects.mate,females.effects.mate,PAR$number.genes)
males.effects.mate <- children.effects$males.effects.mate
females.effects.mate <- children.effects$females.effects.mate
patgenes <- children.effects$patgenes
matgenes <- children.effects$matgenes

if (PAR$save.objects=="no") {remove(children.effects)}
############################################################################






############################################################################
#5.3 CREATE GENERATION CUR PHENOTYPES - ANCESTORS

#make reltype for ancestors
TEMP$reltpe <- PAR$number.generations-currun +7 

phenotype <- make.phenotype(males.effects.mate,females.effects.mate,patgenes,matgenes,TEMP$females.married,TEMP$reltpe,TEMP$children.vector,
                            corrections=correct, founder.pop=FALSE, PARAM=PAR, A.effects.key=a.effects.key,
                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix)

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- phenotype$size.mate

if (PAR$save.objects=="no") {remove(phenotype)}
############################################################################






############################################################################
#5.4 UPDATE MATRICES THAT TRACK CHANGES OVER TIME - ANCESTORS

#update sibling intraclass correlation matrix; for currun=1 it is half of effects.mate because gen 0 sibs aren't related properly
if (PAR$sib.cor=="yes") intra[,currun] <- multi.intra(effects.mate, PAR$parnames[5:length(PAR$parnames)]) 

#update mean/variance/covariance matrix
track.changes <- change.tracker(run=currun,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=TEMP$cor.spouses)

#write these files out to the working directory if PAR$save.gen.data = "yes"
if (PAR$save.gen.data=="yes"){
write.table(t(effects.mate),file=paste("Ancestral.Pop.", currun, sep=""),append=FALSE,row.names=FALSE)}

if (PAR$save.objects=="yes"){
assign(paste("effects.",currun,sep=""),effects.mate)}
############################################################################




#########################################################
#########################################################
#5.55 END ANCESTORS LOOP
}

#how long has the iterating the ancestor populations taken?
PAR$t2 <- Sys.time()
#########################################################
#########################################################

#55555555555555555555555555555555555555555555555555555555555555555555555555555555555555










#6 CREATE SPOUSES OF TWINS & TWIN PARENTS

############################################################################
#Note that all parameters from here on will be those from the last run; i.e., we'll use currun = number.run

#here we just split the population into two - note that (except for possibly a single family)
#entire families go into one population (effects.spouseparents) or the other (effects.twinparents)
TEMP$number.spouse.parents <- floor(TEMP$size.mate/2)
effects.spouseparents <- effects.mate[,1:TEMP$number.spouse.parents]

#here is the second half of the population, which will be the parents of twins & twin sibs
effects.twinparents <- effects.mate[,(TEMP$number.spouse.parents+1):ncol(effects.mate)]  #(see #6.1.2 above)
if (PAR$save.objects=="no") {remove(effects.mate);    } #@@@@@CLEANUP
############################################################################

 


############################################################################
#6.1 ASSORTATIVE MATING - PARENTS OF SPOUSES

#assortative mating
amcur <- assort.mate(effects.spouseparents,VAR$latent.AM,currun,stop.inbreeding=TRUE,mu.constant=PAR$mu.constant)
TEMP$females.married <- amcur$females.married
TEMP$males.married <- amcur$males.married
males.effects.mate <- amcur$males.effects.mate
females.effects.mate <- amcur$females.effects.mate
TEMP$size.mate <- amcur$size.mate
effects.mate <- amcur$effects.mate
TEMP$cor.spouses <- amcur$cor.spouses

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.spouseparents,file=paste("Current.SpouseParents.", currun, sep=""),ncolumns=TEMP$header.cols,append=FALSE)
write(males.effects.mate,file=paste("Current.BreedingMales.SpousesParents.", currun, sep=""),ncolumns=TEMP$header.cols,append=FALSE)
write(females.effects.mate,file=paste("Current.BreedingFemales.SpousesParents.", currun, sep=""),ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that they exist in .RData if PAR$save.objects=="yes"
if (PAR$save.objects=="yes"){
effects.fathers.spouses <- males.effects.mate
effects.mothers.spouses <- females.effects.mate}

if (PAR$save.objects=="no") {remove(effects.spouseparents)   } #@@@@@CLEANUP
############################################################################






############################################################################
#6.2 CREATING CHILDREN - SPOUSES

#number of children needed to have population growth as expected this generation
TEMP$needed.children <- PAR$popsize[currun]/2

#Note: This is the "old way" of doing it;
#INDUCES SELECTION FOR 'YOUTHFUL' LOOKING MALES IN CONTEXT OF AM
#years females have remaining (45-age) to have offspring; after age 45, 0 years remaining
#reproductive.years <- (45-females.effects.mate["age.at.mar",]>0)*(45-females.effects.mate["age.at.mar",]) 
#make lambda (the mean of the poisson dist) a function of reproductive.years
#rate.per.year <- min(1,TEMP$needed.children/sum(reproductive.years))  #the min() assures that females can't average more than 1 child/year
#rate.per.female <- rate.per.year*reproductive.years
#it will induce selection (see line XXX), but only for 2 generations -
#a trade-off for realism with regard to parental age given your age
#years females have remaining (45-age) to have offspring; after age 45, 0 years remaining
TEMP$reproductive.years <- (45-females.effects.mate["age.at.mar",]>0)*(45-females.effects.mate["age.at.mar",]) 

#make lambda (the mean of the poisson dist) a function of TEMP$reproductive.years
TEMP$rate.per.year <- min(1,TEMP$needed.children/sum(TEMP$reproductive.years))  #the min() assures that females can't average more than 1 child/year
TEMP$rate.per.female <- TEMP$rate.per.year*TEMP$reproductive.years

#vector of the actual # children for each marriage; a non-iid poisson process
TEMP$children.vector <- rpois(TEMP$females.married,lambda=TEMP$rate.per.female)
TEMP$children.vector[TEMP$children.vector>9] <- 10          #can never be more than 10 children per female (for computational reasons)

#reproduction
children.effects <- reproduce(TEMP$children.vector,TEMP$females.married,males.effects.mate,females.effects.mate,PAR$number.genes)
males.effects.mate <- children.effects$males.effects.mate
females.effects.mate <- children.effects$females.effects.mate
patgenes <- children.effects$patgenes
matgenes <- children.effects$matgenes

if (PAR$save.objects=="no") {remove(children.effects)}
############################################################################






############################################################################
#6.3 CREATE GENERATION CUR PHENOTYPES - SPOUSES

#make reltype for spouses
TEMP$reltpe <- 6

phenotype <- make.phenotype(males.effects.mate,females.effects.mate,patgenes,matgenes,TEMP$females.married,TEMP$reltpe,TEMP$children.vector,
                            corrections=correct, founder.pop=FALSE, PARAM=PAR, A.effects.key=a.effects.key,
                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix)

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- phenotype$size.mate

if (PAR$save.objects=="no") {remove(phenotype)}
############################################################################







############################################################################
#6.4 UPDATE MATRICES THAT TRACK CHANGES OVER TIME - SPOUSES

#update sibling intraclass correlation matrix; for currun=1 it is half of effects.mate because gen 0 sibs aren't related properly
if (PAR$sib.cor=="yes") {intra[,currun+1] <- multi.intra(effects.mate, PAR$parnames[5:length(PAR$parnames)])  }

#update mean/variance/covariance matrix
track.changes <- change.tracker(run=currun+1,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=TEMP$cor.spouses)


#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file="Current.Spouses",ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="max"
if (PAR$save.objects=="yes"){effects.spouses <- effects.mate}

#begin making the matrix "total" that eventually holds all the relative types
total <- effects.mate
############################################################################

#66666666666666666666666666666666666666666666666666666666666666666666666666666666666666










#7 CREATE PARENTS OF TWIN DATA

############################################################################
#7.1 ASSORTATIVE MATING - PARENTS OF TWINS & TWIN SIBS

#assortative mating
amcur <- assort.mate(effects.twinparents,VAR$latent.AM,currun,stop.inbreeding=TRUE,mu.constant=PAR$mu.constant)
TEMP$females.married <- amcur$females.married
TEMP$males.married <- amcur$males.married
males.effects.mate <- amcur$males.effects.mate
females.effects.mate <- amcur$females.effects.mate
TEMP$size.mate <- amcur$size.mate
#effects.mate <- amcur$effects.mate     #not necessary anymore
TEMP$cor.spouses <- amcur$cor.spouses

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.twinparents,file='Current.TwinParents',ncolumns=TEMP$header.cols,append=FALSE)
write(males.effects.mate,file='Current.BreedingMales.TwinParents',ncolumns=TEMP$header.cols,append=FALSE)
write(females.effects.mate,file='Current.BreedingFemales.TwinParents',ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="max"
if (PAR$save.objects=="yes"){
effects.fathers.twins <- males.effects.mate
effects.mothers.twins <- females.effects.mate}

#Create total matrix which eventually includes all relative types
total <- cbind(total,males.effects.mate,females.effects.mate)
if (PAR$save.objects=="no") {remove(effects.twinparents)   } #@@@@@CLEANUP
############################################################################

#77777777777777777777777777777777777777777777777777777777777777777777777777777777777777









#8 CREATE SIBLINGS OF TWINS DATA

############################################################################
#8.2 NUMBER OF CHILDREN PER MARRIAGE & POPULATION GROWTH - SIBS ONLY

#number of children needed to have population growth as expected this generation
TEMP$needed.children <- (PAR$popsize[currun]/2)-TEMP$females.married #subtract females.married to account for twins - extra kids females will have

#Note: This is the "old way" of doing it; it will induce selection, but only for 2 generations - a trade-off for
TEMP$reproductive.years <- (45-females.effects.mate["age.at.mar",]>0)*(45-females.effects.mate["age.at.mar",]) 

#make lambda (the mean of the poisson dist) a function of TEMP$reproductive.years
TEMP$rate.per.year <- min(1,TEMP$needed.children/sum(TEMP$reproductive.years))  #the min() assures that females can't average more than 1 child/year
TEMP$rate.per.female <- TEMP$rate.per.year*TEMP$reproductive.years

#vector of the actual # sibs, mz's, and dz's for each marriage
TEMP$children.vector <- rpois(TEMP$females.married,lambda=TEMP$rate.per.female)
TEMP$children.vector[TEMP$children.vector>9] <- 10          #can never be more than 10 children per female (for computational reasons)
TEMP$mz.vector <- sample(c(0,2),size=TEMP$females.married,replace=TRUE,prob=c(1-PAR$percent.mz,PAR$percent.mz))
TEMP$dz.vector <- (TEMP$mz.vector-2)*-1

#reproduction
children.effects <- reproduce(TEMP$children.vector,TEMP$females.married,males.effects.mate,females.effects.mate,PAR$number.genes)
paternal.effect <- children.effects$males.effects.mate         #call it a diff name so that you don't erase males.effects.mate
maternal.effect <- children.effects$females.effects.mate       #ditto
patgenes <- children.effects$patgenes
matgenes <- children.effects$matgenes

if (PAR$save.objects=="no") {remove(children.effects)}
############################################################################





############################################################################
#8.3 CREATE GENERATION CUR PHENOTYPES - SIBS
TEMP$reltpe <- 3

phenotype <- make.phenotype(paternal.effect,maternal.effect,patgenes,matgenes,TEMP$females.married,TEMP$reltpe,TEMP$children.vector,
                            corrections=correct, founder.pop=FALSE, PARAM=PAR, A.effects.key=a.effects.key,
                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix) #don't need effects.mate anymore

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- TEMP$number.sib <- phenotype$size.mate
TEMP$twinage <- phenotype$twin.age

if (PAR$save.objects=="no") {remove(phenotype)}
############################################################################





############################################################################
#8.4 UPDATE MATRICES THAT TRACK CHANGES OVER TIME - SIBS

#update sibling intraclass correlation matrix; for currun=1 it is half of effects.mate because gen 0 sibs aren't related properly
if (PAR$sib.cor=="yes") {intra[,currun+4] <- multi.intra(effects.mate, PAR$parnames[5:length(PAR$parnames)])  }

#update mean/variance/covariance matrix
track.changes <- change.tracker(run=currun+4,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=NA)

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file="Current.Sibs",ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="max"
if (PAR$save.objects=="yes"){effects.sibs <- effects.mate}

total <- cbind(total,effects.mate)
############################################################################

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888










#9 CREATE MZ TWIN DATA

############################################################################
#9.2 CHILDREN - MZ TWINS

#reproduction
children.effects <- reproduce(TEMP$mz.vector,TEMP$females.married,males.effects.mate,females.effects.mate,
                              PAR$number.genes,type.mz=TRUE)
paternal.effect <- children.effects$males.effects.mate
maternal.effect <- children.effects$females.effects.mate
patgenes <- children.effects$patgenes
matgenes <- children.effects$matgenes

if (PAR$save.objects=="no") {remove(children.effects)}
############################################################################





############################################################################
#9.3 CREATE GENERATION CUR PHENOTYPES - MZ TWINS
TEMP$reltpe <- 1

phenotype <- make.phenotype(paternal.effect,maternal.effect,patgenes,matgenes,TEMP$females.married,TEMP$reltpe,TEMP$mz.vector,
                            corrections=correct,twin.age=TEMP$twinage, founder.pop=FALSE, PARAM=PAR, A.effects.key=a.effects.key,
                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix) #don't need effects.mate anymore

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- TEMP$number.mz <- phenotype$size.mate
TEMP$empirical.mz <- phenotype$empirical.tw

if (PAR$save.objects=="no") {remove(phenotype)}
############################################################################









############################################################################
#9.4 UPDATE MATRICES THAT TRACK CHANGES OVER TIME - MZ TWINS

#update sibling intraclass correlation matrix; for currun=1 it is half of effects.mate because gen 0 sibs aren't related properly
if (PAR$sib.cor=="yes") {intra[,currun+6] <- multi.intra(effects.mate, PAR$parnames[5:length(PAR$parnames)])  }

#update mean/variance/covariance matrix
track.changes <- change.tracker(run=currun+6,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=NA)

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file="Current.MZs",ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="max"
if (PAR$save.objects=="yes"){effects.mz <- effects.mate}

total <- cbind(total,effects.mate)
############################################################################

#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999











#10 CREATE DZ TWIN DATA


############################################################################
#10.2 CHILDREN - DZ TWINS

#reproduction
children.effects <- reproduce(TEMP$dz.vector,TEMP$females.married,males.effects.mate,females.effects.mate,
                              PAR$number.genes)
paternal.effect <- children.effects$males.effects.mate
maternal.effect <- children.effects$females.effects.mate
patgenes <- children.effects$patgenes
matgenes <- children.effects$matgenes

if (PAR$save.objects=="no") {remove(children.effects)}
############################################################################





############################################################################
#10.3 CREATE GENERATION CUR PHENOTYPES - DZS TWINS
TEMP$reltpe <- 2

phenotype <- make.phenotype(paternal.effect,maternal.effect,patgenes,matgenes,TEMP$females.married,TEMP$reltpe,TEMP$dz.vector,
                            corrections=correct,twin.age=TEMP$twinage, founder.pop=FALSE, PARAM=PAR, A.effects.key=a.effects.key,
                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix) #don't need effects.mate anymore

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- TEMP$number.dz <- phenotype$size.mate
TEMP$empirical.tw <- phenotype$empirical.tw

if (PAR$save.objects=="no") {remove(phenotype)}
############################################################################




############################################################################
#10.4 UPDATE MATRICES THAT TRACK CHANGES OVER TIME - DZ TWINS

#update sibling intraclass correlation matrix; for currun=1 it is half of effects.mate because gen 0 sibs aren't related properly
if (PAR$sib.cor=="yes") {intra[,currun+5] <- multi.intra(effects.mate, PAR$parnames[5:length(PAR$parnames)])  }

#update mean/variance/covariance matrix
track.changes <- change.tracker(run=currun+5,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=NA)

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file='Current.DZs',ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="yes"
if (PAR$save.objects=="yes"){effects.dz <- effects.mate}

#total = dataset of twin parents,twins,sibs, & spouses from the twin's generation (all @ phenotypic values from when married)
total <- cbind(total,effects.mate)

if (PAR$save.objects=="no") {remove(effects.mate,males.effects.mate,females.effects.mate)   } #@@@@@CLEANUP
############################################################################

#10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 











#11 CREATE CHILDREN OF MALE TWINS 

############################################################################
#11.1 ASSORTATIVE MATING - MALE TWINS, FEMALE SPOUSES
#we'll work on female twins below
#we separate them in this way (rather than all males on one side and all females
#on the other) so that we don't have twins marrying twins, which is rare in the real world
sps <- total[,(total["Relative.type",]==6 & total["female",]==1)]  #female spouses
twsps <- total[,(total["Relative.type",]<3 & total["female",]==0)] #male MZ & DZ twins
effects.mate <- cbind(sps,twsps)

#assortative mating
amcur <- assort.mate(effects.mate,VAR$latent.AM,currun,stop.inbreeding=TRUE,mu.constant=PAR$mu.constant)
TEMP$females.married <- amcur$females.married
TEMP$males.married <- amcur$males.married
males.effects.mate <- amcur$males.effects.mate
females.effects.mate <- amcur$females.effects.mate
TEMP$size.mate <- amcur$size.mate
effects.mate <- amcur$effects.mate
TEMP$cor.spouses <- amcur$cor.spouses

#give male twins and female spouses their spouse.id's
TEMP$male.index <- match(males.effects.mate["Subj.ID",],total["Subj.ID",])
total["Spouse.ID",TEMP$male.index] <- males.effects.mate["Spouse.ID",]
TEMP$female.index <- match(females.effects.mate["Subj.ID",],total["Subj.ID",])
total["Spouse.ID",TEMP$female.index] <- females.effects.mate["Spouse.ID",]


#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file='Current.MaleTwin.FemSpouse',ncolumns=TEMP$header.cols,append=FALSE)
write(males.effects.mate,file='Current.BreedingMaleTwins',ncolumns=TEMP$header.cols,append=FALSE)
write(females.effects.mate,file='Current.BreedingFemaleSpouses',ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="yes"
if (PAR$save.objects=="yes"){
effects.fathers.mtc <- males.effects.mate
effects.mothers.mtc <- females.effects.mate}
############################################################################






############################################################################
#11.2 NUMBER OF CHILDREN PER MARRIAGE & POPULATION GROWTH - MALE TWINS CHILDREN (MTC)

#this makes TEMP$needed.childrens - number of children of twins to have population ~ constant after generation
TEMP$needed.children <- .5*PAR$popsize[currun]*((TEMP$number.mz+TEMP$number.dz)/(TEMP$number.mz+TEMP$number.dz+TEMP$number.sib))

#Note: This is the "old way" of doing it; it will induce selection, but only for 2 generations - a trade-off for
TEMP$reproductive.years <- (45-females.effects.mate["age.at.mar",]>0)*(45-females.effects.mate["age.at.mar",]) 

#make lambda (the mean of the poisson dist) a function of TEMP$reproductive.years
TEMP$rate.per.year <- min(1,TEMP$needed.children/sum(TEMP$reproductive.years))  #the min() assures that females can't average more than 1 child/year
TEMP$rate.per.female <- TEMP$rate.per.year*TEMP$reproductive.years

#vector of the actual # children for each marriage; a non-iid poisson process
TEMP$children.vector <- rpois(TEMP$females.married,lambda=TEMP$rate.per.female)
TEMP$children.vector[TEMP$children.vector>9] <- 10          #can never be more than 10 children per female (for computational reasons)

#reproduction
children.effects <- reproduce(TEMP$children.vector,TEMP$females.married,males.effects.mate,females.effects.mate,PAR$number.genes)
males.effects.mate <- children.effects$males.effects.mate         
females.effects.mate <- children.effects$females.effects.mate      
patgenes <- children.effects$patgenes
matgenes <- children.effects$matgenes

if (PAR$save.objects=="no") {remove(children.effects)}
############################################################################





############################################################################
#11.3 CREATE GENERATION CUR PHENOTYPES - MALE TWINS CHILDREN (MTC)
TEMP$reltpe <- 4

phenotype <- make.phenotype(males.effects.mate,females.effects.mate,patgenes,matgenes,TEMP$females.married,TEMP$reltpe,TEMP$children.vector,
                            corrections=correct, founder.pop=FALSE, PARAM=PAR, A.effects.key=a.effects.key,
                            D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix) #don't need effects.mate anymore

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- phenotype$size.mate

if (PAR$save.objects=="no") {remove(phenotype)}
############################################################################





############################################################################
#11.4 UPDATE MATRICES THAT TRACK CHANGES OVER TIME - MTC

#update sibling intraclass correlation matrix; for currun=1 it is half of effects.mate because gen 0 sibs aren't related properly
if (PAR$sib.cor=="yes") {intra[,currun+3] <- multi.intra(effects.mate, PAR$parnames[5:length(PAR$parnames)])  }

#update mean/variance/covariance matrix
track.changes <- change.tracker(run=currun+3,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=NA)

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file='Current.MaleTwinChildren',ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="yes"
if (PAR$save.objects=="yes"){effects.mtc <- effects.mate}

total <- cbind(total,effects.mate)
############################################################################

#11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11











#12 CREATE CHILDREN OF FEMALE TWINS 

############################################################################
#12.1 ASSORTATIVE MATING - FEMALE TWINS, MALE SPOUSES
sps <- total[,(total["Relative.type",]==6 & total["female",]==0)]  #male spouses
twsps <- total[,(total["Relative.type",]<3 & total["female",]==1)] #female MZ & DZ twins
effects.mate <- cbind(sps,twsps)

if (PAR$save.objects=="no") {remove(sps,twsps)   } #@@@@@CLEANUP

#assortative mating
amcur <- assort.mate(effects.mate,VAR$latent.AM,currun,stop.inbreeding=TRUE,mu.constant=PAR$mu.constant)
TEMP$females.married <- amcur$females.married
TEMP$males.married <- amcur$males.married
males.effects.mate <- amcur$males.effects.mate
females.effects.mate <- amcur$females.effects.mate
TEMP$size.mate <- amcur$size.mate
effects.mate <- amcur$effects.mate
TEMP$cor.spouses <- amcur$cor.spouses

#give male twins and female spouses their spouse.id's
TEMP$male.index <- match(males.effects.mate["Subj.ID",],total["Subj.ID",])
total["Spouse.ID",TEMP$male.index] <- males.effects.mate["Spouse.ID",]
TEMP$female.index <- match(females.effects.mate["Subj.ID",],total["Subj.ID",])
total["Spouse.ID",TEMP$female.index] <- females.effects.mate["Spouse.ID",]

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file='Current.FemaleTwin.MaleSpouse',ncolumns=TEMP$header.cols,append=FALSE)
write(males.effects.mate,file='Current.BreedingFemaleTwins',ncolumns=TEMP$header.cols,append=FALSE)
write(females.effects.mate,file='Current.BreedingMaleSpouses',ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="yes"
if (PAR$save.objects=="yes"){
effects.fathers.ftc <- males.effects.mate
effects.mothers.ftc <- females.effects.mate}
############################################################################






############################################################################
#12.2 NUMBER OF CHILDREN PER MARRIAGE & POPULATION GROWTH - FEMALE TWINS CHILDREN (FTC)

#Note: This is the "old way" of doing it; it will induce selection, but only for 2 generations - a trade-off for
TEMP$reproductive.years <- (45-females.effects.mate["age.at.mar",]>0)*(45-females.effects.mate["age.at.mar",]) 

#make lambda (the mean of the poisson dist) a function of TEMP$reproductive.years
TEMP$rate.per.year <- min(1,TEMP$needed.children/sum(TEMP$reproductive.years))  #the min() assures that females can't average more than 1 child/year
TEMP$rate.per.female <- TEMP$rate.per.year*TEMP$reproductive.years

#vector of the actual # children for each marriage; a non-iid poisson process
TEMP$children.vector <- rpois(TEMP$females.married,lambda=TEMP$rate.per.female)
TEMP$children.vector[TEMP$children.vector>9] <- 10          #can never be more than 10 children per female (for computational reasons)

#reproduction
children.effects <- reproduce(TEMP$children.vector,TEMP$females.married,males.effects.mate,females.effects.mate,PAR$number.genes)
males.effects.mate <- children.effects$males.effects.mate        
females.effects.mate <- children.effects$females.effects.mate      
patgenes <- children.effects$patgenes
matgenes <- children.effects$matgenes

if (PAR$save.objects=="no") {remove(children.effects)}
############################################################################





############################################################################
#12.3 CREATE GENERATION CUR PHENOTYPES - FEMALE TWINS CHILDREN (FTC)
TEMP$reltpe <- 5

phenotype <- make.phenotype(males.effects.mate,females.effects.mate,patgenes,matgenes,TEMP$females.married,TEMP$reltpe,TEMP$children.vector,corrections=correct, founder.pop=FALSE, PARAM=PAR, A.effects.key=a.effects.key,D.effects.key=d.effects.key, VARIANCE=VAR, BETA.matrix=beta.matrix) #don't need effects.mate anymore

effects.mate <- phenotype$effects.mate
TEMP$size.mate <- phenotype$size.mate

if (PAR$save.objects=="no") {remove(phenotype)}
############################################################################





############################################################################
#12.4 UPDATE MATRICES THAT TRACK CHANGES OVER TIME - FTC

#update sibling intraclass correlation matrix; for currun=1 it is half of effects.mate because gen 0 sibs aren't related properly
if (PAR$sib.cor=="yes") {intra[,currun+2] <- multi.intra(effects.mate, PAR$parnames[5:length(PAR$parnames)])  }

#update mean/variance/covariance matrix
track.changes <- change.tracker(run=currun+2,x=track.changes,y=effects.mate,mean.names=TEMP$mean.names,effects.names=TEMP$names,
                                beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=TEMP$size.mate,
                                cor.spouses=NA)

#write these files out to the working directory
if (PAR$save.gen.data=="yes"){
write(effects.mate,file='Current.FemaleTwinChildren',ncolumns=TEMP$header.cols,append=FALSE)}

#also, rename these files so that the exist in .RData if PAR$save.objects=="yes"
if (PAR$save.objects=="yes"){effects.ftc <- effects.mate}

total <- cbind(total,effects.mate)

#how long has the iterating taken?
PAR$t3 <- Sys.time()
############################################################################

#12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12










############################################################################
#13 CREATING PEDIGREE DATASETS

#remove unnecessary variables
if (PAR$save.objects=="no") {remove(effects.mate,a.effects.key,aa.effects.key,d.effects.key,maternal.effect,paternal.effect,amcur) } #@@@@@CLEANUP

#Create Pedigree Datasets, one for each measurement timepoint
PedigreeData <- create.pedigree.data(total)



#13.2 START LOOP FOR CREATING REPEATED-MEASURES PEDIGREE DATA SETS 
################
#start loop
for(repmeas in 1:(length(PAR$range.twin.age)/2)){
################

#age everyone, change their phenotype appropriately, and remove those who have died or are outside age range
PData.new <- make.aged.data(PedigreeData,twin.age.range=PAR$range.twin.age,data.age.range=PAR$range.dataset.age,
                            varnames=PAR$varnames,rep=repmeas,corrections=correct,BETA.matrix=beta.matrix)

#use make.data function to turn PData.new into five files in format useful to MX
wide.data <- make.wide.data(PData.new)


#get info on number of missing etc for wide data
TEMP$info.widedata <- widedata.info(wide.data)

if (PAR$real.missing=="yes"){
  twintype <- wide.data$twintype
  TEMP$grp <- c(6,1,1,2,2,3,3,3,3,4,4,5,5,5,5,5,5,5,5)
  TEMP$ozva.miss <- c(0,.002,.724,.846,.729,.942)
  miss.mask <- make.missing.mask(x=wide.data[,2:20],missing=TEMP$ozva.miss,col.grouping=TEMP$grp)
  miss.mask <- cbind(miss.mask,miss.mask[,-1])
  wide.data <- cbind(twintype,wide.data[,2:ncol(wide.data)]*miss.mask)}


#get updated info on number of missing etc for wide data
TEMP$info.widedata.updated <- widedata.info(wide.data)


#create datasets that mimic missing patterns from ozva60k
mz.male <- as.matrix(wide.data[wide.data$twintype=="mzm",2:ncol(wide.data)])
mz.female <- as.matrix(wide.data[wide.data$twintype=="mzf",2:ncol(wide.data)])
dz.male <- as.matrix(wide.data[wide.data$twintype=="dzm",2:ncol(wide.data)])
dz.os <- as.matrix(wide.data[wide.data$twintype=="dzos",2:ncol(wide.data)])
dz.female <- as.matrix(wide.data[wide.data$twintype=="dzf",2:ncol(wide.data)])


#update track.changes for time repmeas
TEMP$index <- match(PData.new[,"Spouse.ID"],PData.new[,"Subj.ID"])
TEMP$cor.spouses <- cor(PData.new[,c("cur.phenotype")],PData.new[TEMP$index,c("cur.phenotype")],use="complete.obs")
track.changes <- change.tracker(run=PAR$number.generations+7+repmeas,x=track.changes,y=t(as.matrix(PData.new[,1:nrow(total)])),mean.names=TEMP$mean.names,effects.names=TEMP$names,beta1=TEMP$track.beta,var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=nrow(PData.new),cor.spouses=TEMP$cor.spouses)


#write out this data as well as track.changes for each measurement for each run number
if (repmeas==1) {

#Track.Changes - for rel types 1-7
track.changes <- change.tracker(run=currun+7,x=track.changes,y=total,mean.names=TEMP$mean.names,effects.names=TEMP$names,beta1=TEMP$track.beta,
                                var.names=TEMP$var.names,cov.names=TEMP$cov.names,popsize=dim(total)[2],cor.spouses=NA)

TEMP$index <- match(PedigreeData[,"Spouse.ID"],PedigreeData[,"Subj.ID"])
track.changes["cor.spouses","reltp.1-7"] <- cor(PedigreeData[,c("cur.phenotype")],PedigreeData[TEMP$index,c("cur.phenotype")],use="complete.obs")
track.changes["cor.spouses","reltp6"] <- cor(wide.data[,"Fa"],wide.data[,"Mo"],use="complete.obs") #track changes

TEMP$tw.sp <- rbind(mz.male[,c("Tw1","sp.tw1")],mz.male[,c("Tw2","sp.tw2")],
                    mz.female[,c("Tw1","sp.tw1")],mz.female[,c("Tw2","sp.tw2")],
                    dz.male[,c("Tw1","sp.tw1")],dz.male[,c("Tw2","sp.tw2")],
                    dz.female[,c("Tw1","sp.tw1")],dz.female[,c("Tw2","sp.tw2")],
                    dz.os[,c("Tw1","sp.tw1")],dz.os[,c("Tw2","sp.tw2")])
track.changes["cor.spouses","reltp2"] <- cor(TEMP$tw.sp,use="complete.obs")[1,2] #track changes
  
write.table(as.data.frame(round(track.changes,3)),"Track.Changes",col.names=TRUE,row.names=TRUE)
if (PAR$sib.cor=="yes"){write.table(as.data.frame(round(intra,3)),"Sib.Correlations",col.names=TRUE,row.names=TRUE)}
}


#create threshold data if necessary
if (is.numeric(PAR$thresholds)){
  mn.all <- mean(unlist(wide.data[,3:20]),na.rm=TRUE)
  sd.all <- sd(unlist(wide.data[,3:20]),na.rm=TRUE)
  thresh <- c(-999,PAR$thresholds,999)

    mz.male[,2:19] <- cut((mz.male[,2:19]-mn.all)/sd.all,breaks=thresh,labels=FALSE) -1
    mz.female[,2:19] <- cut((mz.female[,2:19]-mn.all)/sd.all,breaks=thresh,labels=FALSE) -1
    dz.male[,2:19] <- cut((dz.male[,2:19]-mn.all)/sd.all,breaks=thresh,labels=FALSE) -1
    dz.female[,2:19] <- cut((dz.female[,2:19]-mn.all)/sd.all,breaks=thresh,labels=FALSE) -1
    dz.os[,2:19] <- cut((dz.os[,2:19]-mn.all)/sd.all,breaks=thresh,labels=FALSE) -1

    remove(mn.all,sd.all,thresh)
    }


  
if (repmeas==1) {
#write out for repmeas = 1
write.table(PData.new, file="PedigreeData", col.names=TRUE,row.names=FALSE,na="-999")
write.table(round(mz.male,4), file="MZM",col.names=FALSE,row.names=FALSE,na="-999")
write.table(round(mz.female,4), file="MZF",col.names=FALSE,row.names=FALSE,na="-999")
write.table(round(dz.male,4), file="DZM",col.names=FALSE,row.names=FALSE,na="-999")
write.table(round(dz.female,4), file="DZF",col.names=FALSE,row.names=FALSE,na="-999")
write.table(round(dz.os,4), file="DZOS",col.names=FALSE,row.names=FALSE,na="-999")

}
###################

#write out for repmeas > 1
if (repmeas>1) {
write.table(PData.new, file=paste("PedigreeData.", repmeas,sep=""),col.names=TRUE,row.names=FALSE,na="-999")
write.table(mz.male, file=paste("MZM.",repmeas,sep=""), col.names=TRUE,row.names=FALSE,na="-999")
write.table(mz.female, file=paste("MZF.",repmeas,sep=""), col.names=TRUE,row.names=FALSE,na="-999")
write.table(dz.male, file=paste("DZM.",repmeas,sep=""), col.names=TRUE,row.names=FALSE,na="-999")
write.table(dz.female, file=paste("DZF.",repmeas,sep=""), col.names=TRUE,row.names=FALSE,na="-999")
write.table(dz.os, file=paste("DZOS.",repmeas,sep=""), col.names=TRUE,row.names=FALSE,na="-999")

}


#END DATA CREATION LOOP
################
}
################

remove(repmeas)

############################################################################



#how long has creating datasets taken?
PAR$t4 <- Sys.time()

############################################################################

#13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13












############################################################################
#14 FIND RELATIVE CORRELATIONS & COVARIANCES, AND FIND & WRITE SUMMARY MATRICES


############################################################################
#14.1 FIND RELATIVE CORRELATIONS & COVARIANCES
rel.correlations <- find.rel.correlations(wide.data)

rel.covs <- find.rel.covariances(wide.data)
rel.covs <- c(iter,rel.covs)
names(rel.covs)[1] <- "iteration"
############################################################################




############################################################################
#14.2 FIND CTD VARIANCE COMPONENTS ESTIMATES - to use in graphing
CTD <- list()
{if (rel.correlations[1]/rel.correlations[2] <= 2){
CTD$est.A <- 2*(rel.correlations[1]-rel.correlations[2])
CTD$est.C <- (2*rel.correlations[2]) - rel.correlations[1]
CTD$est.E <- 1-rel.correlations[1]
CTD$est.D <- 0}

else if (rel.correlations[1]/rel.correlations[2] > 2){
CTD$est.A <- (4*rel.correlations[2])-rel.correlations[1]
CTD$est.C <- 0
CTD$est.E <- 1-rel.correlations[1]
CTD$est.D <- (2*rel.correlations[1])-(4*rel.correlations[2])}}
############################################################################




############################################################################
#14.3 CREATE MATRIX TO TRACK CHANGES OF VAR FOR MULTIPLE RUN NUMBERS

#create long PE.Var
TEMP$mult.names <- c(PAR$varnames,"cur.phenotype","mating.phenotype")
TEMP$mult.varnames <- c(paste("mean.",TEMP$mult.names,sep=""),paste("var.",TEMP$mult.names,sep=""))
TEMP$mult.covnames <- c(paste("Cov.",TEMP$mult.names[1],",",TEMP$mult.names[2:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[2],",",TEMP$mult.names[3:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[3],",",TEMP$mult.names[4:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[4],",",TEMP$mult.names[5:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[5],",",TEMP$mult.names[6:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[6],",",TEMP$mult.names[7:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[7],",",TEMP$mult.names[8:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[8],",",TEMP$mult.names[9:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[9],",",TEMP$mult.names[10:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[10],",",TEMP$mult.names[11:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[11],",",TEMP$mult.names[12:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[12],",",TEMP$mult.names[13:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[13],",",TEMP$mult.names[14:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[14],",",TEMP$mult.names[15:length(TEMP$mult.names)],sep=""),
              paste("Cov.",TEMP$mult.names[15],",",TEMP$mult.names[16:length(TEMP$mult.names)],sep=""))
TEMP$col.chooser <- c(1:(PAR$number.generations),(ncol(track.changes)-(length(PAR$range.twin.age)/2)-7),(ncol(track.changes)-(length(PAR$range.twin.age)/2)):ncol(track.changes))
PE.Var <- track.changes[c(TEMP$mult.varnames,TEMP$mult.covnames),TEMP$col.chooser]

#create condensed PE.Var2 - condensed info from only the last time point of empirical data
TEMP$mult.covnames2 <- c(paste("Cov.",PAR$varnames[1],",",PAR$varnames[2:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[2],",",PAR$varnames[3:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[3],",",PAR$varnames[4:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[4],",",PAR$varnames[5:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[5],",",PAR$varnames[6:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[6],",",PAR$varnames[7:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[7],",",PAR$varnames[8:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[8],",",PAR$varnames[9:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[9],",",PAR$varnames[10:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[10],",",PAR$varnames[11:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[11],",",PAR$varnames[12:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[12],",",PAR$varnames[13:length(PAR$varnames)],sep=""),
              paste("Cov.",PAR$varnames[13],",",PAR$varnames[14:length(PAR$varnames)],sep=""))
TEMP$covs <- sum(track.changes[TEMP$mult.covnames2,ncol(track.changes)])*2

PE.Var2 <- matrix(c(track.changes[TEMP$mult.varnames,ncol(track.changes)],TEMP$covs),ncol=1)
rownames(PE.Var2) <- c(TEMP$mult.varnames,"Covs")

#create PE.Var3 - to use in graphing
TEMP$cols <- c(1:(PAR$number.generations+1),ncol(track.changes))
TEMP$E.actual <- colSums(track.changes[c("var.U","var.MZ","var.TW"),TEMP$cols])
PE.Var3 <- rbind(track.changes[c("var.A","var.AA","var.D","var.F","var.S"),TEMP$cols],TEMP$E.actual,
                track.changes[c("var.SEX","var.AGE","var.A.by.SEX","var.A.by.AGE","var.A.by.S","var.A.by.U","Cov.A,F","var.cur.phenotype","cor.spouses"),TEMP$cols])
rownames(PE.Var3) <- c("V(A)","V(AA)","V(D)","V(F)","V(S)","V(E)","V(Sex)","V(Age)","V(AxSex)","V(AxAge)","V(AxS)","V(AxU)","Cov(A,F)","V(P)","r(sps)")


#create long it.summ - last time point + rel covariances 
it.summ <- c(rel.covs,PE.Var[,ncol(PE.Var)])



#Create summary matrices - useful for multiple runs
write.table(as.data.frame(t(PE.Var)),file="GE.Var",col.names=TRUE,row.names=TRUE)

{if (file.exists("Iteration.Summary")) {write.table(t(it.summ),file="Iteration.Summary",col.names=FALSE,row.names=FALSE,append=TRUE)}
         else write.table(t(it.summ),file="Iteration.Summary",col.names=TRUE,row.names=FALSE)} 


############################################################################





############################################################################
#14.4 TIMING

#how long has finding the variance components taken (incl. running MX)?
PAR$t5 <- Sys.time()
TEMP$d <- unlist(strsplit(as.character(Sys.time())," "))
TEMP$date <- c(TEMP$d[1],unlist(strsplit(TEMP$d[2],":"))[1:2])
TEMP$name.date <- paste("GE_",TEMP$date[1],"_",TEMP$date[2],".",TEMP$date[3],sep="")



############################################################################

#14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14


return(list(CTD=CTD,PAR=PAR,PData.new=PData.new,PE.Var=PE.Var,PE.Var2=PE.Var2,PE.Var3=PE.Var3,PedigreeData=PedigreeData,
            TEMP=TEMP,VAR=VAR,beta.matrix=beta.matrix,correct=correct,dz.female=dz.female,dz.male=dz.male,dz.os=dz.os,         
            females.effects.mate=females.effects.mate,males.effects.mate=males.effects.mate,matgenes=matgenes,
            mz.female=mz.female,mz.male=mz.male,parameters=parameters,patgenes=patgenes,
            rel.correlations=rel.correlations,total=total,track.changes=track.changes,
            var.matrix=var.matrix,wide.data=wide.data,it.summ=it.summ))


#CHANGED the below - I think sib.cor doesn't do ANYTHING - the two lists are identical


#jj1 <- list(CTD=CTD,PAR=PAR,PData.new=PData.new,PE.Var=PE.Var,PE.Var2=PE.Var2,PE.Var3=PE.Var3,PedigreeData=PedigreeData,
#            TEMP=TEMP,VAR=VAR,beta.matrix=beta.matrix,correct=correct,dz.female=dz.female,dz.male=dz.male,dz.os=dz.os,         
#            females.effects.mate=females.effects.mate,intra=intra,males.effects.mate=males.effects.mate,matgenes=matgenes,
#            mz.female=mz.female,mz.male=mz.male,parameters=parameters,patgenes=patgenes,
#            rel.correlations=rel.correlations,total=total,track.changes=track.changes,
#            var.matrix=var.matrix,wide.data=wide.data,it.summ=it.summ)

#jj2 <- list(CTD=CTD,PAR=PAR,PData.new=PData.new,PE.Var=PE.Var,PE.Var2=PE.Var2,PE.Var3=PE.Var3,PedigreeData=PedigreeData,
#     TEMP=TEMP,VAR=VAR,beta.matrix=beta.matrix,correct=correct,dz.female=dz.female,dz.male=dz.male,dz.os=dz.os,         
#     females.effects.mate=females.effects.mate,males.effects.mate=males.effects.mate,matgenes=matgenes,
#     mz.female=mz.female,mz.male=mz.male,parameters=parameters,patgenes=patgenes,
#     rel.correlations=rel.correlations,total=total,track.changes=track.changes,
#     var.matrix=var.matrix,wide.data=wide.data,it.summ=it.summ)


#{if (PAR$sib.cor=="yes"){
#return(list(CTD=CTD,PAR=PAR,PData.new=PData.new,PE.Var=PE.Var,PE.Var2=PE.Var2,PE.Var3=PE.Var3,PedigreeData=PedigreeData,
#TEMP=TEMP,VAR=VAR,beta.matrix=beta.matrix,correct=correct,dz.female=dz.female,dz.male=dz.male,dz.os=dz.os,         
#females.effects.mate=females.effects.mate,intra=intra,males.effects.mate=males.effects.mate,matgenes=matgenes,
#mz.female=mz.female,mz.male=mz.male,parameters=parameters,patgenes=patgenes,
#rel.correlations=rel.correlations,total=total,track.changes=track.changes,
#var.matrix=var.matrix,wide.data=wide.data,it.summ=it.summ))}

# else {
#return(list(CTD=CTD,PAR=PAR,PData.new=PData.new,PE.Var=PE.Var,PE.Var2=PE.Var2,PE.Var3=PE.Var3,PedigreeData=PedigreeData,
#TEMP=TEMP,VAR=VAR,beta.matrix=beta.matrix,correct=correct,dz.female=dz.female,dz.male=dz.male,dz.os=dz.os,         
#females.effects.mate=females.effects.mate,males.effects.mate=males.effects.mate,matgenes=matgenes,
#mz.female=mz.female,mz.male=mz.male,parameters=parameters,patgenes=patgenes,
#rel.correlations=rel.correlations,total=total,track.changes=track.changes,
#var.matrix=var.matrix,wide.data=wide.data,it.summ=it.summ))}}

}

