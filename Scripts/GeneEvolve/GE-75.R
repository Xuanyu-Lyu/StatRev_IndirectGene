##################################################################################################
# GeneEvolve version 0.75                                                                        #
# SIMULATION PROGRAM FOR EVOLVING AND CREATING EXTENDED PEDIGREE & TWIN DATA                      #
# by Matthew C Keller, updated April 3, 2012                                                     #
##################################################################################################



#This is newest script from March 20, 2024
#Updated April 18 to allow for any type of latent variable AM (and any type of latent variable parenting, from previous version)

#TO RUN: 1) Create two folder: script.directory (where the scripts will be placed) and output.directory (where the GE output will go)
#        2) Open the GeneEvolve.zip file in your "script.directory" - all scripts must remain in that directory.
#        3) Make changes to the parameters below (#1 of this script) by placing your desired inputs to the *right* of the arrows (<- x  # i.e., user replaces "x").
#        4) Advanced users might consider changing parameters under section 1.5
#        5) Never remove any hash marks (#). Never change anything outside of sections 1 & 1.5
#        6) Run entire script, either in BATCH mode (e.g., R CMD BATCH GeneEvolve69.R on UNIX run line) or by highlighting the entire script and running it in an R session.
#        7) At end, you will have 5 data files in your output.directory (extended families of different twin types) plus a PDF summary of each GeneEvolve run and several matrices





remove(list=ls()) #CHANGED MAR 2024

#1 BASIC USER INPUT PARAMETERS (31)      
#################################################################################################

#BASICS
script.directory <- "/Users/matthewkeller/GoogleCloud/DriveDocuments/RESEARCH/GeneEvolve/GeneEvolve-ETFD/Working"     # directory where all GeneEvolve R scripts exist - NO SPACES in here!
output.directory <- "/Users/matthewkeller/GoogleCloud/DriveDocuments/RESEARCH/GeneEvolve/GeneEvolve-ETFD/GE.Output"     # working directory where all GeneEvolve R output will be written - NO SPACES in here!

make.graphics <- "yes"                # do you want to create a full graphical report at end?

#DEMOGRAPHIC DETAILS:
number.generations <- 10      # number of generations to evolve before pedigree data created
start.pop.size <- 10000       # breeding population size at start of simulation; should be > 500
pop.growth <- rep(1,number.generations)  #vector of length number.generations that tells how big
    # each generation is relative to the one before it; for no growth: rep(1,number.generations); 
    # for constant 5% growth: rep(1.05,number.generations)               

#MODEL DETAILS:
thresholds <- FALSE      # EITHER a vector of thresholds corresponding to location on standard normal OR "FALSE", meaning continuous
am.model <- c(.5,.5,.5,1,1,1,1,1,1,1,0,0,0,0)     #vector of length 14 representing paths from parental P to mating phenotype
                          # In order: A, AA, D, F, S, U, MZ, TW, SEX, AGE, A.by.SEX, A.by.AGE, A.by.S, A.by.U). Common options include:
                          # rep(1,14) is the standard primary phenotypic AM
                          # c(0,0,0,1,1,1,1,1,1,1,0,0,0,0) is social homogamy
                          # c(0,0,0,1,0,0,0,0,0,0,0,0,0,0) is pure cultural homogamy
vt.model <- rep(1,14)     #vector of length 14 representing paths from parental P to parenting phenotype
                           # In order: A, AA, D, F, S, U, MZ, TW, SEX, AGE, A.by.SEX, A.by.AGE, A.by.S, A.by.U). Common options include:
                           # rep(1,14) is the standard parental phenotype to offspring F
                           # c(0,0,0,1,1,1,1,1,1,1,0,0,0,0) is environmental influences of parental phenotype to offspring F
                           # c(0,0,0,1,0,0,0,0,0,0,0,0,0,0) is pure cultural transmission - parental F to offspring F - NOTE: this typically leads to unstable, non-equilibrium values of VF (either increasing to infinity or decreasing to 0)

#DATASET PARAMETERS
percent.mz <- .44              # % of twins who are MZ (in pedigree creation)
range.dataset.age <- c(14,90)  # The pedigree datasets includes twins plus twin offspring,sibs, 
                               # spouses, and parents. This is the range of ages in that dataset
range.twin.age <- c(18,80)     # Range of twin ages in final dataset at each timepoint. Is a 
          # VECTOR twice as long as # of (repeated) measurements
          # Ex. 1, for data of adult twins of all ages measured once: <- c(18,100)
          # Ex. 2, for data of twins measured at exactly 14: <- c(14,14) 
          # Ex. 3, for data of twins measured at 14-18 and again (i.e., REPEATED MEASURES) at 
          #        18-25: <- c(14,18,18,25) 


#VARIANCES: These are variance terms of 1st gen (they may change thereafter). 
# It is nice, but unnecessary, for them to sum to 1

#Genetic Factors:
A  <- .3       # var acct for by A (additive genetic variance)
AA <- .0       # var acct for by AxA (additive-by-additive) epistasis. 
D  <- .1      # var acct for by D (dominance genetic variance)

#Environmental Factors:
U <- .30        # var acct for by U (unique variance). 
               # Note: in non-twins, E (eg, the E from ACE models)=U+MZ+T
MZ <- .0       # var acct for by special MZ twin env (this goes into E for non-mz individuals)
TW <- .0       # var acct for by special twin env (this goes into E for non-twin individuals)
S <- .1       # var acct for by S (sibling environment). 

#Vertical Transmission Variance & Path Coefficients
F <- .20       # this is the var acct for by F in the first generation; it will change thereafter as a function of AM. 
mat2pat <- 1   # "maternal relative to the paternal VT effect" - this is how much greater the maternal path coefficient from maternal phenotype to F (mvt) is than the paternal path coefficient to F (pvt). E.g., if set to 2, the maternal path coefficient is twice the paternal one (mvt=2pvt). Set to 1 for maternal effect = paternal effect (mvt=pvt).
#Note that VF without AM is expected = tau^2(mvt^2 + pvt^2), where tau^2 is the variance of the parenting phenotype. GeneEvolve internally figures out what the actual maternal and paternal path coefficients need to be (given F, mat2pat, and vt.model) to give V(F)=F for the first generation, and in such a way that V(F) is stable over time absent AM (with AM, V(F) tends to increase).
#If users want to input their own specific values of mvt and pvt rather than have GeneEvolve do this for them, simply change mvt and pvt in the final lines under "ADVANCED USER INPUT PARAMETERS" then "OPTIONAL VT INPUT" below. If this is done, the value of F will be VF in the parents of gen0 and mat2pat above will be ignored.


#Covariates & Moderators
AGE <- .0      # var acct for by age; This term can be negative; 
               # sqrt(abs(AGE))= beta(age) = change in phenotype over 1 SD of age (15-45)
SEX <- .0      # var acct for by sex (female=1,male=0); This term can be negative 
               # sqrt(abs(SEX))= beta(sex) = change in phenotype over 1 SD of sex
A.by.AGE <- .0 # var acct for by interaction bw A & Age; 
               # Also, the var of slopes of A across age range specified in range.pop.age
A.by.SEX <- .0 # var acct for by interaction bw A & Sex; 
               # Also, the var of slopes of A across gender. 
A.by.S <- .0   # var acct for by interaction bw A & S (sib env).
A.by.U <- .0   # var acct for by interaction bw A & U (unique env).

#COVARIANCES: (Note: when there are covariances b/w components, total var != sum(variances))
latent.AM <- rep(.2,number.generations) #VECTOR of length = number.generations; 
               # either the correlation (if mu.constant=FALSE) or the mu (if mu.constant=TRUE) b/w spouses' latent mating phenotype each gen., as determined by am.model.
A.S.cor  <- .0 # corr b/w A & S - caused by ACTIVE g-e covariance
A.U.cor <- .0  # corr b/w A & U - caused by ACTIVE g-e covariance
R.Alevel.Aslope.age <- .0 #corr bw intercept & slope. Should be b/w -1 & 1.
R.Alevel.Aslope.sex <- .0 #corr bw intercept & slope. Should be b/w -1 & 1.
R.Alevel.Aslope.S <- .0   #corr bw intercept & slope. Should be b/w -1 & 1.
R.Alevel.Aslope.U <- .0   #corr bw intercept & slope. Should be b/w -1 & 1. 
 # Note on R.Alevel.Aslope.X: these are the corr's b/w the intercepts and slopes of A;
 # e.g., do people with the highest A values also have the highest changes in A across X?
 # if this is 1 (-1), A.by.X becomes the scalar interaction coefficient 
 # - i.e., the PURCEL model where h2 increases (decreases) as X increases
 # if this is 0, A.by.X becomes nonscalar 
 # - i.e., h2 remains ~ constant but diff genes turn on at diff values of X
 # WARNING: for scalar interactions, including only interaction var with no main effect var of A
 # will not produce what is desired!
#################################################################################################










#1.5 ADVANCED USER INPUT PARAMETERS (MULTIPLE RUNS, RANDOM PARAMETERS, RUNNING MX, ETC) (17)
############################################################################

number.runs <- 1    # multiple GE runs- useful for finding bias, var/covar of parameter estimates
rand.parameters <- "no" # do you want GE to use random parameters? 
        # Choices: "no","all","A","D","F","S","AM.mod". "all" randomly chooses A,D,F,S, & AM.mod
save.gen.data <- "no"   # do you want to write out the datasets for each generation? #CHANGED
        # Answer "no" unless you want a detailed ancestral record
continuation <- "no"    # is this a continuation of a job you are now restarting? 
        # Assures the old Parameter.Comparison file is not overwritten
save.rdata <-  "no"    # do you want to save the .RData file?
save.objects <- "no"    # do you want to save all the R objects created?  #CHANGED
        # Generally, you should answer "no" unless you're debugging. Options: "yes" or "no".

number.genes <- 10 # number of genes affecting phenotype; Must be > 2 but <= 50. 5+ is fine 
sibling.age.cor <- .5 # corr b/w sibs' ages within a breeding cohort. 
                      # Higher corr. imply less spacing in years b/w sibs.
range.pop.age <- range.twin.age[1:2]  # range of age in the population over which the beta terms
                      # above explain the AGE and A.by.AGE variation; I recommend placing this 
                      # the same as the first two entries in range.twin.age above

real.missing <- "no"  #do you want missingness patterns to emulate the ozva60k?
run.corr <- "no"      #gets 88 relative covariances using Mx
sib.cor <- "yes"      #do you want GE to compute the correlations b/w sibs for each generation?

gene.model <- "rare"  # this must be equal to "common" or "rare"; most people should choose rare!
                      # If "common", total variance will not usually be sum of var components due
                      # to cor of var components (e.g., A & D)
number.alleles <- 2   # in the "common" model ONLY, how many alleles per locus? 

#Vertical transmission - here GeneEvolve figures out what m & n need to be in order for F to be correct given vt.model and mat2pat
#Do not comment these next 6 lines out
beta.matrix <-matrix(c(sqrt(A),sqrt(AA),sqrt(D),1,sqrt(S),sqrt(U),sqrt(MZ),sqrt(TW),sqrt(SEX),sqrt(AGE),sqrt(A.by.SEX),sqrt(A.by.AGE),sqrt(A.by.S),sqrt(A.by.U)),nrow=1)
par.paths <- beta.matrix*vt.model
par0 <- sum((par.paths^2)*c(rep(1,3),F,rep(1,10))) #the variance of the latent parenting trait in gen0
K <- mat2pat^2 + 1
pvt <- sqrt(F/(K*par0)) #paternal path coefficient - this will be saved in VAR
mvt <- pvt*mat2pat #maternal path coefficient - this will be saved in VAR

#OPTIONAL VT INPUT - to set pvt and mvt to user-specified values. If set to pvt <- pvt (or mvt <- mvt), this has pvt and mvt set by the F value in section 1 (i.e., it doesn't do anything). If you set pvt and mvt to specific values, this will over-ride the value of both "F" (the value set will be the VF in the generation before gen0) and mat2pat input in section 1
pvt <- pvt #.2 #comment this line out if you want to input your own mvt value and override the F value set above
mvt <- mvt #comment this line out if you want to input your own mvt value and override the F value set above
#pvt <- #.5 #uncomment this and set "pvt" to any value if you want to override the F value set above. 
#mvt <- #.5 #uncomment this and set "mvt" to any value if you want to override the F value set above. 


#The below should generally NOT be changed (it should be FALSE). mu is just a trick to be able to provide proper expectations of variances and covariances, and so will decrease naturally over time if the r(sps) is constant and VY is increasing. Forcing mu to be constant isn't something we want to try to simulate (and if you do, it will often create impossible scenarios, such as r(sps) > 1 whenever VY*mu>1)
mu.constant <- FALSE       #the default in GeneEvolve is for the spousal CORRELATION (see "latent.AM" below) to stay constant over generations. However, if you want to make mu constant, change this to "TRUE". If you do this, GeneEvolve will interpret the values of "latent.AM" below to be the mu's rather than the correlations. This is generally NOT recommended.

############################################################################









#################################################################################################
#Start (NO NEED TO ALTER THE SCRIPT BELOW THIS LINE)
#################################################################################################


#2 LOAD FUNCTIONS USED IN SCRIPT & CREATE LISTS OF PARAMETER VALUES ABOVE TO ENTER INTO FUNCTIONS

############################################################################
#2.0 LOAD FUNCTIONS
#change working directory
setwd(script.directory)

#load all PE functions
source("GE.Functions_Short.R")
source("GE.Function_AM.R")
source("GE.Function_Reproduce.R")
source("GE.Function_MakeDataFile.R")
source("GE.Function_MakePhenotype.R")
source("GE.Function_Run.Gene.Evolve75.R")
source("GE.Function_Graph.R")

#change working directory to where output will go
setwd(output.directory)

 
#make PAR and VAR lists; remove unnecessary variables
vars <- c("A", "A.S.cor","A.U.cor","A.by.AGE","A.by.S","A.by.U","AA","AGE","D","F","mvt","pvt","MZ","S","TW",
          "U","SEX","A.by.SEX","R.Alevel.Aslope.S","R.Alevel.Aslope.U","R.Alevel.Aslope.age",
          "R.Alevel.Aslope.sex","am.model","latent.AM","vt.model","make.graphics",
          "number.generations","number.genes","percent.mz","pop.growth","range.dataset.age",
          "range.pop.age","range.twin.age","save.gen.data","save.objects","sib.cor",
          "sibling.age.cor","start.pop.size","output.directory","gene.model","number.alleles","real.missing","save.rdata",
          "continuation","run.corr","beta.matrix","par.paths","par0","K")

PAR1 <- list(output.directory=output.directory, save.gen.data=save.gen.data,
             save.objects=save.objects,sib.cor=sib.cor, 
             make.graphics=make.graphics, number.generations=number.generations,
             start.pop.size=start.pop.size, pop.growth=pop.growth, number.genes=number.genes,
             am.model=am.model, vt.model=vt.model,percent.mz=percent.mz, range.dataset.age=range.dataset.age,
             sibling.age.cor=sibling.age.cor,range.twin.age=range.twin.age,
             range.pop.age=range.pop.age,gene.model=gene.model,
             number.alleles=number.alleles,real.missing=real.missing,
             save.rdata=save.rdata,number.runs=number.runs,continuation=continuation,
             run.corr=run.corr,rand.parameters=rand.parameters,thresholds=thresholds,mu.constant=mu.constant)
          
VAR1 <- list(A=A, AA=AA, D=D,F=F,mvt=mvt, pvt=pvt, S=S, U=U, MZ=MZ, TW=TW, SEX=SEX, AGE=AGE, A.by.SEX=A.by.SEX,
             A.by.AGE=A.by.AGE, A.by.S=A.by.S, A.by.U=A.by.U,latent.AM=latent.AM, A.S.cor=A.S.cor,
             A.U.cor=A.U.cor, R.Alevel.Aslope.sex=R.Alevel.Aslope.sex,
             R.Alevel.Aslope.age=R.Alevel.Aslope.age,R.Alevel.Aslope.S=R.Alevel.Aslope.S,
             R.Alevel.Aslope.U=R.Alevel.Aslope.U)

#remove(list=vars) #CHANGED

#record time
PAR1$t1 <- Sys.time()
#################################################################################################






#################################################################################################
#3 RUN GENEVOLVE & CREATE GRAPHICS

for (iteration in 1:PAR1$number.runs){

 #reset parameters each run if rand.parameters above is not "no". If it is "no", the below just resets VAR and PAR each iteration to their original values
  BOTH   <- alter.parameters(P=PAR1,V=VAR1,max=1)  
  VAR1 <- BOTH$VAR
  PAR1 <- BOTH$PAR
  remove(BOTH)
  
 #Run GeneEvolve
  ALL <- run.gene.evolve(PAR=PAR1,VAR=VAR1,iter=iteration)

 gc()
   
 #Generate Graphics - runs function "Generate.Graphics" which is defined in script "GE.Graph.1.R"
 if (PAR1$make.graphics=="yes") {
  ge.graphics(date.names=ALL$TEMP$name.date, PARAMS=ALL$PAR, VARIANCES=VAR1, 
              PE.Var3=ALL$PE.Var3, PE.Var2=ALL$PE.Var2, rel.correlations=ALL$rel.correlations, 
              track.changes=ALL$track.changes, data.info=ALL$TEMP$info.widedata.updated)}
 }
    
#Save *.RData file if save.rdata=="yes" - names it by date
if (PAR1$save.rdata=="yes"){ save.image(file=paste(ALL$TEMP$name.date,".RData",sep="")) }
#################################################################################################



                                             



