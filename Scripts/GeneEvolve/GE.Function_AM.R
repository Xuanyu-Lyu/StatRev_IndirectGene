

############################################################################
#ASSORTATIVE MATING FUNCTION

assort.mate <- function(effects.mate,latent.AM,currun,stop.inbreeding=TRUE,mu.constant=FALSE){
  
#For debugging:
#  latent.AM=VAR$latent.AM;currun=1;stop.inbreeding=TRUE;mu.constant=PAR$mu.constant
  

#probability that each individual marries for generation currun; capped at .95
prob.marry <- min(.95, rnorm(1,mean=.85,sd=.04))
size.mate <- ncol(effects.mate)

### MARRIAGEABLE PEOPLE - this section creates an equal number of males and females who will be paired off below
#We need to make the number of marriable people the same in females & males
num.females.mate <- sum(effects.mate["female",])
num.males.mate  <- -sum(effects.mate["female",]-1)
number.married <- sum((runif(size.mate)< rep(prob.marry,size.mate))*1 )
if (odd(number.married)==TRUE){ number.married <- number.married +  sample(c(-1,1),1)} #make number.married always even

#make males.married & females.married the same number; the smallest of number.married/2, males.mate or females.mate
couples <- min(num.males.mate,num.females.mate,number.married/2) 

#samplers to be used for males & females
female.sampler <- c(rep(1,couples),rep(0,num.females.mate-couples))
male.sampler <- c(rep(1,couples),rep(0,num.males.mate-couples))

#include whether each person is marriageable
effects.mate["marriageable",(effects.mate["female",]==1)] <- sample(female.sampler,size=num.females.mate,replace=FALSE)
effects.mate["marriageable",(effects.mate["female",]==0)] <- sample(male.sampler,size=num.males.mate,replace=FALSE)


### MATING VIA PRIMARY ASSORTMENT - ANCESTORS
#split up the effects.mate matrix
males.effects.mate <- effects.mate[,((effects.mate["female",]==0)& (effects.mate["marriageable",]==1))]
females.effects.mate <- effects.mate[,((effects.mate["female",]==1)& (effects.mate["marriageable",]==1))]

#order the males & females by their mating phenotypic value
males.effects.mate <- males.effects.mate[,order(males.effects.mate["mating.phenotype",])]
females.effects.mate <- females.effects.mate[,order(females.effects.mate["mating.phenotype",])]

#making the template correlation matrix and ordering the two variables
#as it is (empirical=TRUE), the correlation is *exactly* the AM coefficient each generation; may change to FALSE to make it more realistic
#nevertheless, there is a stochastic element to it due to the sorting
AM.mate <- latent.AM[currun]

#This is where we change the correlation, AM.mate, to be what it needs to be to ensure that mu = latent.AM each generation. Basically, if the V(P) increases, the spousal correlation increases as well. Note that latent.AM is the correlation when mu.constant=FALSE, and it is mu when latent.AM=TRUE. However, AM.mate is alway the correlation.
if(mu.constant){AM.mate <- latent.AM[currun]*sd(males.effects.mate["mating.phenotype",])*sd(females.effects.mate["mating.phenotype",])}
template.AM.dist <- mvrnormal(n=couples,mu=c(0,0),Sigma=matrix(c(1,AM.mate,AM.mate,1),nrow=2),empirical=TRUE)
rank.template.males <- rank(template.AM.dist[,1])
rank.template.females <- rank(template.AM.dist[,2])

#marry the males and females according to the assortative mating coefficient earlier inputed
males.effects.mate <-  males.effects.mate[,rank.template.males]
females.effects.mate <-  females.effects.mate[,rank.template.females]


###AVOID INBREEDING
if (stop.inbreeding==TRUE){
  
#for now, we just drop the inbreeding pairs; an alternative is to remate them
no.sib.inbreeding <- males.effects.mate["Father.ID",]!=females.effects.mate["Father.ID",]
no.cousin.inbreeding <- {males.effects.mate["Fathers.Father.ID",]!=females.effects.mate["Fathers.Father.ID",]&
                      males.effects.mate["Fathers.Father.ID",]!=females.effects.mate["Mothers.Father.ID",]&
                      males.effects.mate["Mothers.Father.ID",]!=females.effects.mate["Fathers.Father.ID",]&
                      males.effects.mate["Mothers.Father.ID",]!=females.effects.mate["Mothers.Father.ID",]&
                      males.effects.mate["Fathers.Mother.ID",]!=females.effects.mate["Fathers.Mother.ID",]&
                      males.effects.mate["Fathers.Mother.ID",]!=females.effects.mate["Mothers.Mother.ID",]&
                      males.effects.mate["Mothers.Mother.ID",]!=females.effects.mate["Fathers.Mother.ID",]&
                      males.effects.mate["Mothers.Mother.ID",]!=females.effects.mate["Mothers.Mother.ID",]}

no.inbreeding <- (no.sib.inbreeding & no.cousin.inbreeding)


### CREATE FINAL MARRIED PAIRS WHO AREN'T INBRED - ANCESTORS
#for now, we just drop inbred pairs & create the male & female files (both in order of married pairs, so 1st column married 1st column)
#note that there can be no inbreeding when currun==1, and only sib-sib inbreeding when currun==2
{if (currun==1){
     males.effects.mate <- males.effects.mate[,no.sib.inbreeding]
     females.effects.mate <- females.effects.mate[,no.sib.inbreeding]}
else {
     males.effects.mate <- males.effects.mate[,no.inbreeding]
     females.effects.mate <- females.effects.mate[,no.inbreeding]}}

#End avoid inbreeding
}



### CREATE COR.SPOUSES & "SPOUSE.ID" VARIABLE
females.married <- males.married <- dim(males.effects.mate)[2]
cor.spouses <- cor(males.effects.mate["cur.phenotype",],females.effects.mate["cur.phenotype",])

males.effects.mate["Spouse.ID",] <- females.effects.mate["Subj.ID",]
females.effects.mate["Spouse.ID",] <- males.effects.mate["Subj.ID",]

male.index <- match(males.effects.mate["Subj.ID",],effects.mate["Subj.ID",])
effects.mate["Spouse.ID",male.index] <- males.effects.mate["Spouse.ID",]
female.index <- match(females.effects.mate["Subj.ID",],effects.mate["Subj.ID",])
effects.mate["Spouse.ID",female.index] <- females.effects.mate["Spouse.ID",]

#make a list of what is returned from the function
list(females.married=females.married,males.married=males.married,males.effects.mate=males.effects.mate,
females.effects.mate=females.effects.mate,size.mate=size.mate,effects.mate=effects.mate, cor.spouses=cor.spouses)


#END FUNCTION
}

############################################################################
