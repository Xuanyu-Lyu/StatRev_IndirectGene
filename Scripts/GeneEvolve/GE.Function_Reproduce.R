
############################################################################
#REPRODUCE FUNCTION


reproduce <- function(children.vector,number.couples,males.effects.mate,females.effects.mate,number.genes=PAR$number.genes,
                      type.mz=FALSE){

#next 3 lines just for debugging
#children.vector <- mz.vector
#number.couples <- females.married
#type.mz <- TRUE

#number of children per marriage & population growth 
#number of children for each married pair - this is stochastic (poisson process),
#so in any given generation, pop size may be larger or smaller than expected
#However, since we've defined expected population sizes up front, population sizes
#will NOT drift far away from what we've specified (this was a problem beforehand)

#number of children needed to have population growth as expected this generation
{if (type.mz==TRUE) number.children <- sum(children.vector)/2
else number.children <- sum(children.vector)}
  
#which rows have the locus information?
loc1.row <- which(dimnames(males.effects.mate)[[1]]=="Pat.Loc1") 
lastloc.row <- loc1.row + (number.genes*2) -1


### PATERNAL MEIOSIS
#make matrix of paternal effect, with fathers having as many columns as they have children
males.effects.mate <- males.effects.mate[,rep(1:number.couples,children.vector)]

#vectorize the paternal genes in 1xnumber.genes*2*number.children size vector
{if (type.mz==TRUE)
    { #make matrix of paternal effect, with fathers having ONE column for pair of MZ twins (only one sperm made them both!)
        paternal.gamete.effect <- males.effects.mate[,seq(1,ncol(males.effects.mate),2)]
      #vectorize the paternal genes in 1xnumber.genes*2*number.mz size vector
        pat.vector <- as.vector(paternal.gamete.effect[loc1.row:lastloc.row,])}
else   pat.vector <- as.vector(males.effects.mate[loc1.row:lastloc.row,])}

#figure out the indices of the above vector that randomly chooses one of the two genes per locus to pass on per child
chosen <- sample(c(1,(number.genes+1)),(number.genes*number.children),replace=TRUE)
add1 <- rep((0:(number.genes-1)),number.children)
add2 <- rep(seq(from=0,by=(number.genes*2),length.out=number.children),each=number.genes)
pat.chosen.index <- chosen+add1+add2

#choose these genes and place in a matrix
{if (type.mz==TRUE)
   {patgenes.init <-  matrix(pat.vector[pat.chosen.index],nrow=number.genes,ncol=number.children)
    patgenes <- patgenes.init[,rep(1:number.children,each=2)]}
else
patgenes <-  matrix(pat.vector[pat.chosen.index],nrow=number.genes,ncol=number.children)}



### MATERNAL MEIOSIS
#make matrix of maternal effect, with fathers having as many columns as they have children
females.effects.mate <- females.effects.mate[,rep(1:number.couples,children.vector)]

#vectorize the paternal genes in 1xnumber.genes*2*number.children size vector
{if (type.mz==TRUE)
    { #make matrix of paternal effect, with fathers having ONE column for pair of MZ twins (only one sperm made them both!)
        maternal.gamete.effect <- females.effects.mate[,seq(1,ncol(males.effects.mate),2)]
      #vectorize the paternal genes in 1xnumber.genes*2*number.mz size vector
        mat.vector <- as.vector(maternal.gamete.effect[loc1.row:lastloc.row,])}
else   mat.vector <- as.vector(females.effects.mate[loc1.row:lastloc.row,])}

#figure out the indices of the above vector that randomly chooses one of the two genes per locus to pass on per child
chosen <- sample(c(1,(number.genes+1)),(number.genes*number.children),replace=TRUE)
add1 <- rep((0:(number.genes-1)),number.children)
add2 <- rep(seq(from=0,by=(number.genes*2),length.out=number.children),each=number.genes)
mat.chosen.index <- chosen+add1+add2

#choose these genes and place in a matrix
{if (type.mz==TRUE)
   {matgenes.init <-  matrix(mat.vector[mat.chosen.index],nrow=number.genes,ncol=number.children)
    matgenes <- matgenes.init[,rep(1:number.children,each=2)]}
else
matgenes <-  matrix(mat.vector[mat.chosen.index],nrow=number.genes,ncol=number.children)}

#return
list(males.effects.mate=males.effects.mate,females.effects.mate=females.effects.mate,patgenes=patgenes,matgenes=matgenes,
     children.vector=children.vector)

#END FUNCTION
}

############################################################################







