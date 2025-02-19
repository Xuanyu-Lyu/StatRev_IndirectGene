#Updated Nov 16, 2023 to change round() in rel.correlations to 3 instead of 2





# Function create.pedigree.data - takes in total matrix and turns it into pedigree data w/ correct famids
create.pedigree.data <- function(TOT){

#mz twins
mz.twins <- TOT[,TOT["Relative.type",]==1]
Famid <- mz.twins["Father.ID",]
mz.twins  <- as.data.frame(cbind(Famid,t(mz.twins)),row.names=colnames(mz.twins))
twintype <- mz.twins$Relative.type*mz.twins$female+1
mz.twins <- cbind(mz.twins,twintype)

#dz twins
dz.twins <- TOT[,TOT["Relative.type",]==2]
Famid <- dz.twins["Father.ID",]
dz.twins  <- as.data.frame(cbind(Famid,t(dz.twins)),row.names=colnames(dz.twins))
dz.twins <- dz.twins[order(dz.twins$Famid),]
twintype <- rep(by(dz.twins$female,INDICES=dz.twins$Famid,FUN=sum),each=2)+3  #dz mm = 3, os = 4, ff=5
dz.twins <- cbind(dz.twins,twintype)

#sibs of twins
sibs <- TOT[,TOT["Relative.type",]==3]
Famid <- sibs["Father.ID",]
twintype <- 0
sibs  <- as.data.frame(cbind(Famid,t(sibs),twintype),row.names=colnames(sibs))

#spouses of twins
twins <- TOT[,TOT["Relative.type",]<3]
spouses <- TOT[,TOT["Relative.type",]==6 & TOT["Spouse.ID",] != 0]
index <- match(spouses["Spouse.ID",],twins["Subj.ID",])
Famid <- twins["Father.ID",index]
spouses <- as.data.frame(cbind(Famid,t(spouses),twintype),row.names=colnames(spouses))
twins <- as.data.frame(t(twins),row.names=colnames(twins))

#parents of twins
twin.fathers <- TOT[,TOT["Relative.type",]==7 & TOT["female",]==0]
twin.mothers <- TOT[,TOT["Relative.type",]==7 & TOT["female",]==1]
Famid <- twin.fathers["Subj.ID",]
twin.fathers <- as.data.frame(cbind(Famid,t(twin.fathers),twintype),row.names=colnames(twin.fathers))
twin.mothers <- as.data.frame(cbind(Famid,t(twin.mothers),twintype),row.names=colnames(twin.mothers)) #mothers are in the same order as fathers

#children of twins
maletwin.children <-  TOT[,TOT["Relative.type",]==4]
Famid <- maletwin.children["Fathers.Father.ID",]
maletwin.children <- as.data.frame(cbind(Famid,t(maletwin.children),twintype),row.names=colnames(maletwin.children))
femaletwin.children <-  TOT[,TOT["Relative.type",]==5]
Famid <- femaletwin.children["Mothers.Father.ID",]
femaletwin.children <- as.data.frame(cbind(Famid,t(femaletwin.children),twintype),row.names=colnames(femaletwin.children))

#start Pedigree Data Set
PedigreeData <- rbind(mz.twins,dz.twins,sibs,spouses,twin.fathers,twin.mothers,maletwin.children,femaletwin.children)
PedigreeData <- PedigreeData[order(PedigreeData$Famid),-which(names(mz.twins)=="marriageable")]   #order dataset by famid and remove marriageable variable
twins.sibs <- rbind(mz.twins,dz.twins,sibs)

#figure out when each person's death age is
PedigreeData$death.age <- (rbeta(nrow(PedigreeData),5.8,1.6)*60)+29

#change Relative.type
PedigreeData[PedigreeData$Relative.type==1,"reltype"] <- "mz"
PedigreeData[PedigreeData$Relative.type==2,"reltype"] <- "dz"
PedigreeData[PedigreeData$Relative.type==3,"reltype"] <- "sib"
PedigreeData[(PedigreeData$Relative.type==4 | PedigreeData$Relative.type==5),"reltype"] <- "child"
PedigreeData[PedigreeData$Relative.type==6,"reltype"] <- "spouse"
PedigreeData[PedigreeData$Relative.type==7,"reltype"] <- "parent"

#figure out whether each family is of type MZ male, MZ female, DZ male, DZ female, or DZ os
famids <- unique(PedigreeData$Famid)
num.families <- length(famids)
num.per.family <- table(PedigreeData$Famid)
index <- match(famids,PedigreeData$Famid) #returns index for first rel type in each family; usually the twin, but not always
PedigreeData$fam.twintype <- rep(PedigreeData[index,"twintype"],times=num.per.family)
PedigreeData[PedigreeData$fam.twintype==1,"famtype"] <- "mzm"
PedigreeData[PedigreeData$fam.twintype==2,"famtype"] <- "mzf"
PedigreeData[PedigreeData$fam.twintype==3,"famtype"] <- "dzm"
PedigreeData[PedigreeData$fam.twintype==4,"famtype"] <- "dzos"
PedigreeData[PedigreeData$fam.twintype==5,"famtype"] <- "dzf"

#remove families with no twins
PedigreeData <- PedigreeData[PedigreeData$fam.twintype >0 & PedigreeData$fam.twintype <6,]

return(PedigreeData)
}

############################################################################






#Function make.aged.data - takes ped data (each line is individual) and ages each individual appropriately within families
make.aged.data <- function(x,twin.age.range,data.age.range,varnames,rep,corrections,BETA.matrix){

# create ages of individuals in datasets as function of range.twin.age 
#twin ages - these are the ages of the twins, and for now, are the same within families
num.per.family <- table(x$Famid)
twn.age <- runif(length(x[x$twintype != 0,1])/2,twin.age.range[rep*2-1],twin.age.range[rep*2])
twn.age.fam <- rep(twn.age,times=as.vector(num.per.family))
x$cur.age <- twn.age.fam

#create ages of siblings - ages relative to twins
age <-  rnorm(nrow(x),0,3)  #make normally dist var of age - how much older/younger sibs are
age[age<0] <- age[age<0]-1  #make at least a year spacing between twins & sibs
age[age>=0] <- age[age>=0]+1  #make at least a year spacing between twins & sibs
x[x$Relative.type==3,"cur.age"] <- x[x$Relative.type==3,"cur.age"] + age[x$Relative.type==3]

#create ages of parents, making parents older by X amount from oldest sib
max.sib.ages <- aggregate(x$cur.age,list(id=x$Famid),max)  #oldest twin/sib per family
parent.ages <- rgamma(length(num.per.family),shape=4,scale=3)+16
parent.ages[parent.ages<18] <- 18  #make sure parents are at least 18 years older 
parent.ages <- max.sib.ages[,2]+parent.ages
parent.ages <- rep(parent.ages,times=num.per.family)
x[(x$Relative.type==7 & x$female==0),"cur.age"] <- parent.ages[(x$Relative.type==7 & x$female==0)]
x[(x$Relative.type==7 & x$female==1),"cur.age"] <- (parent.ages[(x$Relative.type==7 & x$female==1)] +
  rnorm(length(x[(x$Relative.type==7 & x$female==1),"cur.age"]),-3,3))

#create ages of spouses - ages relative to twins; make wives 3 years younger on avg, husbands 3 years older
age <-  rnorm(nrow(x),0,3)  #make normally dist var of age - how much older/younger spouses are
x[(x$Relative.type==6 & x$female==0),"cur.age"] <- x[(x$Relative.type==6 & x$female==0),"cur.age"] +
  age[(x$Relative.type==6 & x$female==0)] +  rnorm(length(x[(x$Relative.type==6 & x$female==0),"cur.age"]),3,3)
x[(x$Relative.type==6 & x$female==1),"cur.age"] <- x[(x$Relative.type==6 & x$female==1),"cur.age"] +
  age[(x$Relative.type==6 & x$female==1)] +  rnorm(length(x[(x$Relative.type==6 & x$female==1),"cur.age"]),-3,3)

#create ages of twin children - ages relative to twins
age <-  x$cur.age - rgamma(nrow(x),shape=4,scale=3)-14  #make gamma dist var of age - how much older/younger offspring are
x[(x$Relative.type==4 | x$Relative.type==5),"cur.age"] <- age[(x$Relative.type==4 | x$Relative.type==5)]

#remove spouses who haven't yet met and family members who are dead, under data.age.range[1], or over data.age.range[2] (set by parameter at beginning of script)
x <- x[(x$cur.age<x$age.at.mar & x$Relative.type==6)==FALSE,]
x <- x[x$cur.age<x$death.age,]  #remove people who are dead
x <- x[x$cur.age>data.age.range[1],]  #remove people who are too young
x <- x[x$cur.age<data.age.range[2],]  #remove people who would be too old to participate, even if alive

#remove families with no twins
famids <- unique(x$Famid)
num.families <- length(famids)
num.per.family <- table(x$Famid)
index <- match(famids,x$Famid) #returns index for first rel type in each family; usually the twin, but not always
dummy <- rep(x[index,"twintype"],times=num.per.family)
x <- x[dummy >0 & dummy <6,]

#sort the data
x <- x[order(x$Famid,x$Relative.type,x$cur.age),]

#create new phenotypes depending on ages
x$AGE <- (x$cur.age-corrections$age.mean)*corrections$age.sd
x$A.by.AGE <- x$AGE*x$A.slopes.age
x$cur.phenotype <- as.vector(BETA.matrix %*% as.matrix(t(x[,varnames])))

return(x)
}

############################################################################







#Function make.wide.data - takes ped data (each line is individual) and turns into wide data (each line is family)
make.wide.data <- function(pdata){

#1 ARRANGE DATA WITHIN FAMILIES
#create newrel & reduce pdata
newrel <-(pdata$Relative.type==1)*1 +(pdata$Relative.type==2)*1 +
         (pdata$Relative.type==3)*3 +(pdata$Relative.type==4)*4 +
         (pdata$Relative.type==5)*4 + (pdata$Relative.type>5)*pdata$Relative.type

pdata <- cbind(pdata[,c("Famid","Subj.ID","Father.ID","Mother.ID","Spouse.ID","female",
                            "cur.age","cur.phenotype","twintype","fam.twintype")],newrel) #variables of interest for now
names(pdata)[c(8,11)] <- c("phenotype","newrel")

#define "twin1" & "twin 2" @ random within families
dd <- which(pdata$newrel==1)
twinord <- rep(c(1,2),times=(ceiling(length(dd)/2)))[1:length(dd)]
pdata[(pdata$newrel==1),"newrel"] <- twinord

#now males first in DZOS twins (newrel=1 for male and 2 for female DZOS's)
x.os <- pdata[pdata$twintype==4,]
os.sex <- rep(c(0,1),length(unique(x.os$Famid)))
os.famid <- rep(unique(x.os$Famid),each=2)
os.twinord <- rep(c(1,2),length(unique(x.os$Famid)))
os.dat <- data.frame(Famid=os.famid,female=os.sex,new=os.twinord)
x.os <- merge(x.os,os.dat,by.x=c("Famid","female"),by.y=c("Famid","female"),all.x=TRUE,all.y=TRUE,sort=TRUE)
x.os <- na.omit(x.os[,c("Subj.ID","new")])
pdata[match(x.os$Subj.ID,pdata$Subj.ID),"newrel"] <- x.os$new   
pdata <- pdata[order(pdata$Famid,pdata$newrel),]

#define brothers within families
bros.per.family <- table(pdata[pdata$newrel==3 & pdata$female==0,]$Famid)
bronums <- 30+stack(sapply(bros.per.family,sample,replace=FALSE))[,1]
pdata[pdata$newrel==3  & pdata$female==0,"newrel"] <- bronums

#define sisters within families
sis.per.family <- table(pdata[pdata$newrel==3 & pdata$female==1,]$Famid)
sisnums <- 40+stack(sapply(sis.per.family,sample,replace=FALSE))[,1]
pdata[pdata$newrel==3  & pdata$female==1,"newrel"] <- sisnums

#define spouses within families
sp.id <- pdata[pdata$newrel==6,"Spouse.ID"]
index <- match(sp.id,pdata$Subj.ID)
spousenums <- 50+pdata[index,"newrel"]
pdata[pdata$newrel==6,"newrel"] <- spousenums
pdata <- pdata[(is.na(pdata$newrel)==FALSE),] #remove spouses with no living twins

#define parents within families: fathers=60,mothers=61
if (length(pdata[pdata$newrel==7,1])>0){
pdata[pdata$newrel==7,"newrel"] <- pdata[pdata$newrel==7,"female"]+60}

#define children within families
kid.fa.id <- pdata[pdata$newrel==4,"Father.ID"]
twinid <- (pdata$newrel<3)*pdata$Subj.ID
kid.fa.index <- match(kid.fa.id,twinid)
kid.fa.index[is.na(kid.fa.index)] <- 0

kid.mo.id <- pdata[pdata$newrel==4,"Mother.ID"]
kid.mo.index <- match(kid.mo.id,twinid)
kid.mo.index[is.na(kid.mo.index)] <- 0

index <- kid.fa.index+kid.mo.index
index[index==0] <- NA

tw1malekid <- 70*(pdata[index,"newrel"]==1 & pdata[pdata$newrel==4,"female"]==0)
tw1femalekid <- 80*(pdata[index,"newrel"]==1 & pdata[pdata$newrel==4,"female"]==1)
tw2malekid <- 90*(pdata[index,"newrel"]==2 & pdata[pdata$newrel==4,"female"]==0)
tw2femalekid <- 100*(pdata[index,"newrel"]==2 & pdata[pdata$newrel==4,"female"]==1)

kidnums <- tw1malekid+tw1femalekid+tw2malekid+tw2femalekid
pdata[pdata$newrel==4,"newrel"] <- kidnums
pdata <- pdata[(is.na(pdata$newrel)==FALSE),] #remove children with no living twin parent

#define male children of twin1 in families
if (length(pdata[pdata$newrel==70,1])>0){
tw1male.per.family <- table(pdata[pdata$newrel==70,]$Famid)
if (length(unique((pdata[pdata$newrel==70,"Famid"])))>1) tw1malenums <- 70+stack(sapply(tw1male.per.family,sample,replace=FALSE))[,1]
if (length(unique((pdata[pdata$newrel==70,"Famid"])))==1) tw1malenums <- 70+sapply(tw1male.per.family,sample,replace=FALSE)
pdata[pdata$newrel==70,"newrel"] <- tw1malenums}

#define female children of twin1 in families
if (length(pdata[pdata$newrel==80,1])>0){
tw1female.per.family <- table(pdata[pdata$newrel==80,]$Famid)
if (length(unique((pdata[pdata$newrel==80,"Famid"])))>1) tw1femalenums <- 80+stack(sapply(tw1female.per.family,sample,replace=FALSE))[,1]
if (length(unique((pdata[pdata$newrel==80,"Famid"])))==1) tw1femalenums <- 80+sapply(tw1female.per.family,sample,replace=FALSE)
pdata[pdata$newrel==80,"newrel"] <- tw1femalenums}

#define male children of twin2 in families
if (length(pdata[pdata$newrel==90,1])>0){
tw2male.per.family <- table(pdata[pdata$newrel==90,]$Famid)
if (length(unique((pdata[pdata$newrel==90,"Famid"])))>1) tw2malenums <- 90+stack(sapply(tw2male.per.family,sample,replace=FALSE))[,1]
if (length(unique((pdata[pdata$newrel==90,"Famid"])))==1) tw2malenums <- 90+sapply(tw2male.per.family,sample,replace=FALSE)
pdata[pdata$newrel==90,"newrel"] <- tw2malenums}

#define female children of twin2 in families
if (length(pdata[pdata$newrel==100,1])>0){
tw2female.per.family <- table(pdata[pdata$newrel==100,]$Famid)
if (length(unique((pdata[pdata$newrel==100,"Famid"])))>1) tw2femalenums <- 100+stack(sapply(tw2female.per.family,sample,replace=FALSE))[,1]
if (length(unique((pdata[pdata$newrel==100,"Famid"])))==1) tw2femalenums <- 100+sapply(tw2female.per.family,sample,replace=FALSE)
pdata[pdata$newrel==100,"newrel"] <- tw2femalenums}

#sort this dataset
pdata <- pdata[order(pdata$Famid,pdata$newrel),]

#KEY FOR PDATA RELATIVE TYPE (newrel):
#     1 = twin 1 (mz or dz)
#     2 = twin 2 (mz or dz)
#     31-39 = Male Siblings (in rand order)
#     41-49 = Female Siblings (in rand order)
#     51 = Spouse of twin 1
#     52 = Spouse of twin 2
#     60 = Twins' father
#     61 = Twins' mother
#     71-79 = Male children of twin 1
#     81-89 = Female children of twin 1
#     91-99 = Male children of twin 2
#     101-109 = Female children of twin 2


#2 CREATE WIDE FORMATED DATA
#first make "long format" cascade data 
famids <- unique(pdata$Famid)
num.families <- length(famids)
long.data <- as.data.frame(matrix(0,nrow=num.families*18,ncol=4))
long.data[,1] <- 1:nrow(long.data)
long.data[,2] <- rep(famids,each=18)
long.data[,3] <- rep(c(1,2,60,61,31,32,41,42,51,52,71,72,81,82,91,92,101,102),num.families)
dimnames(long.data)[[2]] <- c("order","Famid","newrel","phenotype")
long.data <- merge(long.data,pdata,by.x=c("Famid","newrel"),by.y=c("Famid","newrel"),all.x=TRUE,all.y=TRUE,sort=FALSE)
names(long.data) <- c("Famid","newrel","order","phenotype.x","Subj.ID","Father.ID","Mother.ID","Spouse.ID","female","cur.age","phenotype.y","twintype","fam.twintype")
long.data <- long.data[is.na(long.data$order)==FALSE,]  #remove extra relatives (e.g., 3+ bros or sisters, 3+ male children of twin 1, etc)
long.data <- long.data[order(long.data$order),]
long.data <- round(long.data[,c("Famid","Subj.ID","newrel","phenotype.y","female","cur.age","fam.twintype","twintype")],3)

#now make "wide format" data
famids <- unique(long.data$Famid)
index <- match(famids,long.data$Famid) #returns index for first rel type in each family
fam.twintype <- ifelse(is.na(long.data[index,"fam.twintype"]),long.data[index+1,"fam.twintype"],long.data[index,"fam.twintype"])
fam.twintype <- as.character(fam.twintype)
fam.twintype[fam.twintype=="1"] <- "mzm"
fam.twintype[fam.twintype=="2"] <- "mzf"
fam.twintype[fam.twintype=="3"] <- "dzm"
fam.twintype[fam.twintype=="4"] <- "dzos"
fam.twintype[fam.twintype=="5"] <- "dzf"

wide.data <- cbind(fam.twintype, reshape(long.data[,c("Famid","newrel","phenotype.y")],idvar="Famid",direction="wide",timevar="newrel"),
              reshape(long.data[,c("Famid","newrel","cur.age")],idvar="Famid",direction="wide",timevar="newrel")[,-1])

names(wide.data) <- c("twintype","Famid","Tw1","Tw2","Fa","Mo","bro1","bro2","sis1","sis2","sp.tw1","sp.tw2",
                "son1.tw1","son2.tw1","dau1.tw1","dau2.tw1", "son1.tw2","son2.tw2","dau1.tw2","dau2.tw2",
                "Tw1.age","Tw2.age","Fa.age","Mo.age","bro1.age","bro2.age","sis1.age","sis2.age","sp.tw1.age",
                "sp.tw2.age", "son1.tw1.age","son2.tw1.age","dau1.tw1.age","dau2.tw1.age", "son1.tw2.age",
                "son2.tw2.age","dau1.tw2.age","dau2.tw2.age")

return(wide.data)
}









#Function data.info - get basic info about sample sizes in datasets
widedata.info <- function(x){
  x <- as.data.frame(x)

n.mzm <- length(na.omit(x[x$twintype=="mzm","Tw1"])) + length(na.omit(x[x$twintype=="mzm","Tw2"]))
n.mzf <- length(na.omit(x[x$twintype=="mzf","Tw1"])) + length(na.omit(x[x$twintype=="mzf","Tw2"]))
n.dzm <- length(na.omit(x[x$twintype=="dzm","Tw1"])) + length(na.omit(x[x$twintype=="dzm","Tw2"]))
n.dzf <- length(na.omit(x[x$twintype=="dzf","Tw1"])) + length(na.omit(x[x$twintype=="dzf","Tw2"]))
n.dzos <- length(na.omit(x[x$twintype=="dzos","Tw1"])) + length(na.omit(x[x$twintype=="dzos","Tw2"]))
n.twins <- n.mzm+n.mzf+n.dzm+n.dzos+n.dzf 

n.tw.fa <-   length(na.omit(x[,"Fa"])) 
n.tw.mo <-   length(na.omit(x[,"Mo"]))

n.sps <-   length(na.omit(x[,"sp.tw1"]))+length(na.omit(x[,"sp.tw2"]))
            
n.tw.ch <- length(na.omit(x[,"son1.tw1"])) + length(na.omit(x[,"son2.tw1"])) + length(na.omit(x[,"dau1.tw1"])) +length(na.omit(x[,"dau2.tw1"])) +
  length(na.omit(x[,"son1.tw2"])) + length(na.omit(x[,"son2.tw2"])) + length(na.omit(x[,"dau1.tw2"])) +length(na.omit(x[,"dau2.tw2"]))

n.sibs <-    length(na.omit(x[,"bro1"])) + length(na.omit(x[,"bro2"])) + length(na.omit(x[,"sis1"])) + length(na.omit(x[,"sis2"]))

n.tot <- nrow(x)*c(2,2,4,2,8)

ours <- c(n.twins,n.tw.fa+n.tw.mo,n.sibs,n.sps,n.tw.ch)
pct.nonmiss <- ours/n.tot
names(pct.nonmiss) <- c("twins","parents.twins","sibs.twins","spouses.twins","children.twins")

#make a vector that has comparable numbers to va30k
va30k <- c(14761,3360,4800,3184,4391)
twin.scaler <- va30k[1]/ours[1]
scaled.ours <- ours*twin.scaler
pct.nonmiss.va30k <- va30k/(n.tot*twin.scaler)

#make final list to return
DATA.INFO <- list(n.mzm=n.mzm,n.mzf=n.mzf,n.dzm=n.dzm,n.dzf=n.dzf,n.dzos=n.dzos,n.twins=n.twins,n.tw.fa=n.tw.fa,n.tw.mo=n.tw.mo,
                  n.sps=n.sps,n.tw.ch=n.tw.ch,n.sibs=n.sibs,n.tot=n.tot,n.data=ours,pct.nonmissing=pct.nonmiss,va30k=va30k,pct.nonmiss.va30k=pct.nonmiss.va30k)

return(DATA.INFO)
}




#Function make.missing - induce a pattern of missingness in dataset for each column
make.missing.mask <- function(x,missing,col.grouping=1:ncol(x)){

#This function tries to make exactly the % missing in each column as defined in 'missing' 
#if more NAs already exist than are specified, no additional missing values are created
#col.grouping is a vector telling which columns should be grouped together for purposes of missingness. Defaults to each
#column being independent
  
  if(length(missing) != length(unique(col.grouping))) stop("The 'missing' argument must have as many elements as there are unique elements in col.grouping")
  if(max(missing)>1 | min(missing) < 0) stop("The 'missing' argument must contain percentages (between 0 & 1) of missing data  for each column of x")
  if(length(col.grouping) != ncol(x)) stop("The 'col.grouping' argument must have as many elements as there are columns in dataset x")

missing.mask <- matrix(1,nrow=nrow(x),ncol=ncol(x))
missing.mask[is.na(x)] <- NA  

  k <- 0
  for (i in unique(col.grouping)){
    k <- k+1
    index <- which(col.grouping==i)
    sub.x <- as.matrix(x[,index])
    tot.elements <- nrow(sub.x)*ncol(sub.x)
    already.missing <- sum(is.na(sub.x))
    non.missing <- tot.elements-already.missing
    num.needed.miss <- floor(tot.elements*(missing[k]))
    make.miss <- num.needed.miss-already.missing
    make.miss <- (make.miss>0)*make.miss

    y <- sample(c(rep(1,non.missing-make.miss),rep(NA,make.miss)),size=non.missing,replace=FALSE)
   
   ind2 <-  which(is.na(sub.x)==FALSE,arr.ind=TRUE)
   missing.mask[,index][ind2] <- y
   }
  dimnames(missing.mask)[[2]] <- names(x)
  return(missing.mask)
}








#Function find.rel.correlations
#NOTE: this does not remove the mean from the columns; mean differences between generations will produce 'artificially' high correlations
#Note also that the correct way to do this is using intraclass correlations - but it makes bugger all difference

find.rel.correlations <- function(widedata){

twintype <- widedata$twintype  
widedata <- as.matrix(widedata[,2:ncol(widedata)])  
  
mzs <- widedata[twintype=="mzm" | twintype=="mzf",]
dzs <- widedata[twintype=="dzm" | twintype=="dzos" | twintype=="dzf",]
sibs <- na.omit(rbind(widedata[,c("Tw1","bro1")],widedata[,c("Tw1","bro2")],widedata[,c("Tw1","sis1")],widedata[,c("Tw1","sis2")],
              widedata[,c("Tw2","bro1")],widedata[,c("Tw2","bro2")],widedata[,c("Tw2","sis1")],widedata[,c("Tw2","sis2")],
              widedata[,c("bro1","bro2")],widedata[,c("bro1","sis1")],widedata[,c("bro1","sis2")],widedata[,c("bro2","sis1")],widedata[,c("bro2","sis2")],
              widedata[,c("son1.tw1","son2.tw1")],widedata[,c("son1.tw1","dau1.tw1")],widedata[,c("son1.tw1","dau2.tw1")],
              widedata[,c("son2.tw1","dau1.tw1")],widedata[,c("son2.tw1","dau2.tw1")],widedata[,c("son2.tw1","dau1.tw1")],
              widedata[,c("son1.tw2","son2.tw2")],widedata[,c("son1.tw2","dau1.tw2")],widedata[,c("son1.tw2","dau2.tw2")],
              widedata[,c("son2.tw2","dau1.tw2")],widedata[,c("son2.tw2","dau2.tw2")],widedata[,c("son2.tw2","dau1.tw2")]))
par.child <- na.omit(rbind(widedata[,c("Tw1","son1.tw1")],widedata[,c("Tw1","son2.tw1")],widedata[,c("Tw1","dau1.tw1")],widedata[,c("Tw1","dau2.tw1")],
                   widedata[,c("Tw2","son1.tw2")],widedata[,c("Tw2","son2.tw2")],widedata[,c("Tw2","dau1.tw2")],widedata[,c("Tw2","dau2.tw2")],
                   widedata[,c("Tw1","Fa")],widedata[,c("Tw1","Mo")],widedata[,c("Tw2","Fa")],widedata[,c("Tw2","Mo")],
                   widedata[,c("bro1","Fa")],widedata[,c("bro1","Mo")],widedata[,c("bro2","Fa")],widedata[,c("bro2","Mo")],
                   widedata[,c("sis1","Fa")],widedata[,c("sis1","Mo")],widedata[,c("sis2","Fa")],widedata[,c("sis2","Mo")])) #something screwy here
granp <- na.omit(rbind(widedata[,c("Fa","son1.tw1")],widedata[,c("Fa","son2.tw1")],widedata[,c("Fa","dau1.tw1")],widedata[,c("Fa","dau2.tw1")],
               widedata[,c("Fa","son1.tw2")],widedata[,c("Fa","son2.tw2")],widedata[,c("Fa","dau1.tw2")],widedata[,c("Fa","dau2.tw2")],
               widedata[,c("Mo","son1.tw1")],widedata[,c("Mo","son2.tw1")],widedata[,c("Mo","dau1.tw1")],widedata[,c("Mo","dau2.tw1")],
               widedata[,c("Mo","son1.tw2")],widedata[,c("Mo","son2.tw2")],widedata[,c("Mo","dau1.tw2")],widedata[,c("Mo","dau2.tw2")]))   
mz.avnc <- na.omit(rbind(mzs[,c("Tw1","son1.tw2")],mzs[,c("Tw1","son2.tw2")],mzs[,c("Tw1","dau1.tw2")],mzs[,c("Tw1","dau2.tw2")],
                 mzs[,c("Tw2","son1.tw1")],mzs[,c("Tw2","son2.tw1")],mzs[,c("Tw2","dau1.tw1")],mzs[,c("Tw2","dau2.tw1")]))
avnc <- na.omit(rbind(dzs[,c("Tw1","son1.tw2")],dzs[,c("Tw1","son2.tw2")],dzs[,c("Tw1","dau1.tw2")],dzs[,c("Tw1","dau2.tw2")],
              dzs[,c("Tw2","son1.tw1")],dzs[,c("Tw2","son2.tw1")],dzs[,c("Tw2","dau1.tw1")],dzs[,c("Tw2","dau2.tw1")],
              widedata[,c("bro1","son1.tw2")],widedata[,c("bro1","son2.tw2")],widedata[,c("bro1","dau1.tw2")],widedata[,c("bro1","dau2.tw2")],
              widedata[,c("bro1","son1.tw1")],widedata[,c("bro1","son2.tw1")],widedata[,c("bro1","dau1.tw1")],widedata[,c("bro1","dau2.tw1")],
              widedata[,c("bro2","son1.tw2")],widedata[,c("bro2","son2.tw2")],widedata[,c("bro2","dau1.tw2")],widedata[,c("bro2","dau2.tw2")],
              widedata[,c("bro2","son1.tw1")],widedata[,c("bro2","son2.tw1")],widedata[,c("bro2","dau1.tw1")],widedata[,c("bro2","dau2.tw1")],
              widedata[,c("sis1","son1.tw2")],widedata[,c("sis1","son2.tw2")],widedata[,c("sis1","dau1.tw2")],widedata[,c("sis1","dau2.tw2")],
              widedata[,c("sis1","son1.tw1")],widedata[,c("sis1","son2.tw1")],widedata[,c("sis1","dau1.tw1")],widedata[,c("sis1","dau2.tw1")],
              widedata[,c("sis2","son1.tw2")],widedata[,c("sis2","son2.tw2")],widedata[,c("sis2","dau1.tw2")],widedata[,c("sis2","dau2.tw2")],
              widedata[,c("sis2","son1.tw1")],widedata[,c("sis2","son2.tw1")],widedata[,c("sis2","dau1.tw1")],widedata[,c("sis2","dau2.tw1")]))
mz.cous <- na.omit(rbind(mzs[,c("son1.tw1","son1.tw2")],mzs[,c("son1.tw1","son2.tw2")],mzs[,c("son1.tw1","dau1.tw2")],mzs[,c("son1.tw1","dau2.tw2")],
                 mzs[,c("son2.tw1","son1.tw2")],mzs[,c("son2.tw1","son2.tw2")],mzs[,c("son2.tw1","dau1.tw2")],mzs[,c("son2.tw1","dau2.tw2")],
                 mzs[,c("dau1.tw1","dau1.tw2")],mzs[,c("dau1.tw1","son2.tw2")],mzs[,c("dau1.tw1","dau1.tw2")],mzs[,c("dau1.tw1","dau2.tw2")],
                 mzs[,c("dau2.tw1","son1.tw2")],mzs[,c("dau2.tw1","dau2.tw2")],mzs[,c("dau2.tw1","dau1.tw2")],mzs[,c("dau2.tw1","dau2.tw2")]))
cous <- na.omit(rbind(dzs[,c("son1.tw1","son1.tw2")],dzs[,c("son1.tw1","son2.tw2")],dzs[,c("son1.tw1","dau1.tw2")],dzs[,c("son1.tw1","dau2.tw2")],
                 dzs[,c("son2.tw1","son1.tw2")],dzs[,c("son2.tw1","son2.tw2")],dzs[,c("son2.tw1","dau1.tw2")],dzs[,c("son2.tw1","dau2.tw2")],
                 dzs[,c("dau1.tw1","dau1.tw2")],dzs[,c("dau1.tw1","son2.tw2")],dzs[,c("dau1.tw1","dau1.tw2")],dzs[,c("dau1.tw1","dau2.tw2")],
                 dzs[,c("dau2.tw1","son1.tw2")],dzs[,c("dau2.tw1","dau2.tw2")],dzs[,c("dau2.tw1","dau1.tw2")],dzs[,c("dau2.tw1","dau2.tw2")]))
spouses <- na.omit(rbind(widedata[,c("Tw1","sp.tw1")],widedata[,c("Tw2","sp.tw2")],widedata[,c("Fa","Mo")]))

rel.correlations <- rep(NA,10)
rel.correlations[1] <- cor(mzs[,c("Tw1","Tw2")],use="pairwise.complete")[1,2]
rel.correlations[2] <- cor(widedata[twintype=="dzm" | twintype=="dzos" | twintype=="dzf",c("Tw1","Tw2")],use="pairwise.complete")[1,2]
try(rel.correlations[3] <-  cor(sibs)[1,2])
try(rel.correlations[4] <- cor(par.child)[1,2])
try(rel.correlations[5] <- cor(granp)[1,2])
try(rel.correlations[6] <- cor(mz.avnc)[1,2])
try(rel.correlations[7] <- cor(avnc)[1,2])
try(rel.correlations[8] <- cor(mz.cous)[1,2])
try(rel.correlations[9]<- cor(cous)[1,2])
try(rel.correlations[10]<-  cor(spouses)[1,2])
rel.correlations <- round(rel.correlations,3)
names(rel.correlations) <- c("MZ","DZ","sib","par.child","granp","mz.avunc","avunc","mz.cous","cous","spouse")

return(rel.correlations)
}
                      
                      
############################################################################










#Function find.rel.correlations
#NOTE: this does not remove the mean from the columns; mean differences between generations will produce 'artificially' high correlations
#Note also that the correct way to do this is using intraclass correlations - but it makes bugger all difference

find.rel.covariances <- function(widedata){

twintype <- widedata$twintype  
widedata <- as.matrix(widedata[,2:ncol(widedata)])  
  
mzs <- widedata[twintype=="mzm" | twintype=="mzf",]
dzs <- widedata[twintype=="dzm" | twintype=="dzos" | twintype=="dzf",]
sibs <- na.omit(rbind(widedata[,c("Tw1","bro1")],widedata[,c("Tw1","bro2")],widedata[,c("Tw1","sis1")],widedata[,c("Tw1","sis2")],
              widedata[,c("Tw2","bro1")],widedata[,c("Tw2","bro2")],widedata[,c("Tw2","sis1")],widedata[,c("Tw2","sis2")],
              widedata[,c("bro1","bro2")],widedata[,c("bro1","sis1")],widedata[,c("bro1","sis2")],widedata[,c("bro2","sis1")],widedata[,c("bro2","sis2")],
              widedata[,c("son1.tw1","son2.tw1")],widedata[,c("son1.tw1","dau1.tw1")],widedata[,c("son1.tw1","dau2.tw1")],
              widedata[,c("son2.tw1","dau1.tw1")],widedata[,c("son2.tw1","dau2.tw1")],widedata[,c("son2.tw1","dau1.tw1")],
              widedata[,c("son1.tw2","son2.tw2")],widedata[,c("son1.tw2","dau1.tw2")],widedata[,c("son1.tw2","dau2.tw2")],
              widedata[,c("son2.tw2","dau1.tw2")],widedata[,c("son2.tw2","dau2.tw2")],widedata[,c("son2.tw2","dau1.tw2")]))
par.child <- na.omit(rbind(widedata[,c("Tw1","son1.tw1")],widedata[,c("Tw1","son2.tw1")],widedata[,c("Tw1","dau1.tw1")],widedata[,c("Tw1","dau2.tw1")],
                   widedata[,c("Tw2","son1.tw2")],widedata[,c("Tw2","son2.tw2")],widedata[,c("Tw2","dau1.tw2")],widedata[,c("Tw2","dau2.tw2")],
                   widedata[,c("Tw1","Fa")],widedata[,c("Tw1","Mo")],widedata[,c("Tw2","Fa")],widedata[,c("Tw2","Mo")],
                   widedata[,c("bro1","Fa")],widedata[,c("bro1","Mo")],widedata[,c("bro2","Fa")],widedata[,c("bro2","Mo")],
                   widedata[,c("sis1","Fa")],widedata[,c("sis1","Mo")],widedata[,c("sis2","Fa")],widedata[,c("sis2","Mo")])) #something screwy here
granp <- na.omit(rbind(widedata[,c("Fa","son1.tw1")],widedata[,c("Fa","son2.tw1")],widedata[,c("Fa","dau1.tw1")],widedata[,c("Fa","dau2.tw1")],
               widedata[,c("Fa","son1.tw2")],widedata[,c("Fa","son2.tw2")],widedata[,c("Fa","dau1.tw2")],widedata[,c("Fa","dau2.tw2")],
               widedata[,c("Mo","son1.tw1")],widedata[,c("Mo","son2.tw1")],widedata[,c("Mo","dau1.tw1")],widedata[,c("Mo","dau2.tw1")],
               widedata[,c("Mo","son1.tw2")],widedata[,c("Mo","son2.tw2")],widedata[,c("Mo","dau1.tw2")],widedata[,c("Mo","dau2.tw2")]))   
mz.avnc <- na.omit(rbind(mzs[,c("Tw1","son1.tw2")],mzs[,c("Tw1","son2.tw2")],mzs[,c("Tw1","dau1.tw2")],mzs[,c("Tw1","dau2.tw2")],
                 mzs[,c("Tw2","son1.tw1")],mzs[,c("Tw2","son2.tw1")],mzs[,c("Tw2","dau1.tw1")],mzs[,c("Tw2","dau2.tw1")]))
avnc <- na.omit(rbind(dzs[,c("Tw1","son1.tw2")],dzs[,c("Tw1","son2.tw2")],dzs[,c("Tw1","dau1.tw2")],dzs[,c("Tw1","dau2.tw2")],
              dzs[,c("Tw2","son1.tw1")],dzs[,c("Tw2","son2.tw1")],dzs[,c("Tw2","dau1.tw1")],dzs[,c("Tw2","dau2.tw1")],
              widedata[,c("bro1","son1.tw2")],widedata[,c("bro1","son2.tw2")],widedata[,c("bro1","dau1.tw2")],widedata[,c("bro1","dau2.tw2")],
              widedata[,c("bro1","son1.tw1")],widedata[,c("bro1","son2.tw1")],widedata[,c("bro1","dau1.tw1")],widedata[,c("bro1","dau2.tw1")],
              widedata[,c("bro2","son1.tw2")],widedata[,c("bro2","son2.tw2")],widedata[,c("bro2","dau1.tw2")],widedata[,c("bro2","dau2.tw2")],
              widedata[,c("bro2","son1.tw1")],widedata[,c("bro2","son2.tw1")],widedata[,c("bro2","dau1.tw1")],widedata[,c("bro2","dau2.tw1")],
              widedata[,c("sis1","son1.tw2")],widedata[,c("sis1","son2.tw2")],widedata[,c("sis1","dau1.tw2")],widedata[,c("sis1","dau2.tw2")],
              widedata[,c("sis1","son1.tw1")],widedata[,c("sis1","son2.tw1")],widedata[,c("sis1","dau1.tw1")],widedata[,c("sis1","dau2.tw1")],
              widedata[,c("sis2","son1.tw2")],widedata[,c("sis2","son2.tw2")],widedata[,c("sis2","dau1.tw2")],widedata[,c("sis2","dau2.tw2")],
              widedata[,c("sis2","son1.tw1")],widedata[,c("sis2","son2.tw1")],widedata[,c("sis2","dau1.tw1")],widedata[,c("sis2","dau2.tw1")]))
mz.cous <- na.omit(rbind(mzs[,c("son1.tw1","son1.tw2")],mzs[,c("son1.tw1","son2.tw2")],mzs[,c("son1.tw1","dau1.tw2")],mzs[,c("son1.tw1","dau2.tw2")],
                 mzs[,c("son2.tw1","son1.tw2")],mzs[,c("son2.tw1","son2.tw2")],mzs[,c("son2.tw1","dau1.tw2")],mzs[,c("son2.tw1","dau2.tw2")],
                 mzs[,c("dau1.tw1","dau1.tw2")],mzs[,c("dau1.tw1","son2.tw2")],mzs[,c("dau1.tw1","dau1.tw2")],mzs[,c("dau1.tw1","dau2.tw2")],
                 mzs[,c("dau2.tw1","son1.tw2")],mzs[,c("dau2.tw1","dau2.tw2")],mzs[,c("dau2.tw1","dau1.tw2")],mzs[,c("dau2.tw1","dau2.tw2")]))
cous <- na.omit(rbind(dzs[,c("son1.tw1","son1.tw2")],dzs[,c("son1.tw1","son2.tw2")],dzs[,c("son1.tw1","dau1.tw2")],dzs[,c("son1.tw1","dau2.tw2")],
                 dzs[,c("son2.tw1","son1.tw2")],dzs[,c("son2.tw1","son2.tw2")],dzs[,c("son2.tw1","dau1.tw2")],dzs[,c("son2.tw1","dau2.tw2")],
                 dzs[,c("dau1.tw1","dau1.tw2")],dzs[,c("dau1.tw1","son2.tw2")],dzs[,c("dau1.tw1","dau1.tw2")],dzs[,c("dau1.tw1","dau2.tw2")],
                 dzs[,c("dau2.tw1","son1.tw2")],dzs[,c("dau2.tw1","dau2.tw2")],dzs[,c("dau2.tw1","dau1.tw2")],dzs[,c("dau2.tw1","dau2.tw2")]))
spouses <- na.omit(rbind(widedata[,c("Tw1","sp.tw1")],widedata[,c("Tw2","sp.tw2")],widedata[,c("Fa","Mo")]))


mz.inlaw <- na.omit(rbind(mzs[,c("Tw1","sp.tw2")],mzs[,c("Tw2","sp.tw1")]))
sib.inlaw <- na.omit(rbind(widedata[,c("bro1","sp.tw1")],widedata[,c("bro1","sp.tw2")],widedata[,c("bro2","sp.tw1")],widedata[,c("bro2","sp.tw2")],
                     widedata[,c("sis1","sp.tw1")],widedata[,c("sis1","sp.tw2")],widedata[,c("sis2","sp.tw1")],widedata[,c("sis2","sp.tw2")],
                     dzs[,c("Tw1","sp.tw2")],dzs[,c("Tw2","sp.tw1")]))
mz.sps.inlaw <- na.omit(mzs[,c("sp.tw1","sp.tw2")])
sib.sps.inlaw <- na.omit(dzs[,c("sp.tw1","sp.tw2")])
par.inlaw <-  na.omit(rbind(widedata[,c("Fa","sp.tw1")],widedata[,c("Fa","sp.tw2")],widedata[,c("Mo","sp.tw1")],widedata[,c("Mo","sp.tw2")]))
mz.avnc.inlaw <- na.omit(rbind(mzs[,c("sp.tw1","son1.tw2")],mzs[,c("sp.tw1","son2.tw2")],mzs[,c("sp.tw1","dau1.tw2")],mzs[,c("sp.tw1","dau2.tw2")],
                 mzs[,c("sp.tw2","son1.tw1")],mzs[,c("sp.tw2","son2.tw1")],mzs[,c("sp.tw2","dau1.tw1")],mzs[,c("sp.tw2","dau2.tw1")]))
avnc.inlaw <- na.omit(rbind(dzs[,c("sp.tw1","son1.tw2")],dzs[,c("sp.tw1","son2.tw2")],dzs[,c("sp.tw1","dau1.tw2")],dzs[,c("sp.tw1","dau2.tw2")],
              dzs[,c("sp.tw2","son1.tw1")],dzs[,c("sp.tw2","son2.tw1")],dzs[,c("sp.tw2","dau1.tw1")],dzs[,c("sp.tw2","dau2.tw1")]))




rel.covariances <- rep(NA,17)
rel.covariances[1] <- cov(mzs[,c("Tw1","Tw2")],use="pairwise.complete")[1,2]
rel.covariances[2] <- cov(dzs[,c("Tw1","Tw2")],use="pairwise.complete")[1,2]
try(rel.covariances[3] <-  cov(sibs)[1,2])
try(rel.covariances[4] <- cov(par.child)[1,2])
try(rel.covariances[5] <- cov(granp)[1,2])
try(rel.covariances[6] <- cov(mz.avnc)[1,2])
try(rel.covariances[7] <- cov(avnc)[1,2])
try(rel.covariances[8] <- cov(mz.cous)[1,2])
try(rel.covariances[9]<- cov(cous)[1,2])
try(rel.covariances[10]<-  cov(spouses)[1,2])
try(rel.covariances[11]<-  cov(mz.inlaw)[1,2])
try(rel.covariances[12]<-  cov(sib.inlaw)[1,2])
try(rel.covariances[13]<-  cov(mz.sps.inlaw)[1,2])
try(rel.covariances[14]<-  cov(sib.sps.inlaw)[1,2])
try(rel.covariances[15]<-  cov(par.inlaw)[1,2])
try(rel.covariances[16]<-  cov(mz.avnc.inlaw)[1,2])
try(rel.covariances[17]<-  cov(avnc.inlaw)[1,2])

rel.covariances <- round(rel.covariances,4)
names(rel.covariances) <- c("MZ","DZ","sib","par.child","granp","mz.avunc","avunc","mz.cous","cous","spouse",
                            "mz.inlaw","sib.inlaw","mz.sps.inlaw","sib.sps.inlaw","par.inlaw","mz.avnc.inlaw","avnc.inlaw")
return(rel.covariances)
}
                      
                      
############################################################################

