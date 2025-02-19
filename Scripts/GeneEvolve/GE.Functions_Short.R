#LIST OF SHORT FUNCTIONS


#FIRST FUNCTION - FROM THE "MASS" LIBRARY (this way you don't have to load MASS each time for a single function)
mvrnormal <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) {
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, sym = TRUE, EISPACK = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
        t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
        nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
        drop(X)
    else t(X)}


#SECOND FUNCTION:
#this function creates a  "gene.matrix" with dimensions row=PAR$number.genes col=# individuals * 2(for each chromosome)
#for now, we can use the "common" model, where there are two common alleles at each locus (.1 < MAF <.5),
#or a "rare" model, where there are 10,000 equally uncommon alleles at each locus

make.gene.matrix <- function(gene.number=PAR$number.genes,gm=PAR$gene.model,psize=PAR$start.pop.size,allele.num=PAR$number.alleles){
  gene.matrix.pat <- matrix(0,nrow=gene.number,ncol=psize)
  gene.matrix.mat <- matrix(0,nrow=gene.number,ncol=psize)
  x <- runif(gene.number,min=.1,max=.9)
  if (gene.number > 51) {stop("You must input 50 genes or fewer")} 
  
  if (gm == "common") {
     probvect <- matrix(rbind(x,1-x),nrow=2)
     if (allele.num > 100) {stop("You must input fewer than 100 alleles per locus in the common model")}
          gene.template <- matrix(ceiling(runif(gene.number*2,min=0,max=2000)),nrow=2, ncol=gene.number)   
    
     for (i in 1:gene.number){gene.matrix.pat[i,] <-   sample(gene.template[,i],size=psize,replace=TRUE,prob=probvect[,i])}
     for (i in 1:gene.number){gene.matrix.mat[i,] <-   sample(gene.template[,i],size=psize,replace=TRUE,prob=probvect[,i])}

     freq.descrip <- probvect}

  else if (gm == "rare") {
     #this makes a choice of allele labels for each locus (column), and makes certain that the same label is never reused
     #(if psize > 10^4, this will reuse labels, but do so a minimum number of times)
      gt.reps <- max(1,ceiling((gene.number*psize*2*2)/10^4))
      gene.template <- rep(1:10^4,gt.reps)
      gene.matrix <- matrix(sample(gene.template,gene.number*psize*2,replace=FALSE),nrow=gene.number)
      gene.matrix.pat <- gene.matrix[,1:psize]
      gene.matrix.mat <- gene.matrix[,(psize+1):(psize*2)]

      freq.descrip <- NULL
      gene.template <- NULL}

  else {stop("You must enter either 'common' or 'rare' in PAR$gene.model entry")}
  list(data.pat=gene.matrix.pat,data.mat=gene.matrix.mat,pq=freq.descrip,labels=gene.template)
}


#THIRD FUNCTION - FROM THE GTOOLS LIBRARY
odd <- function(x){x != as.integer(x/2) * 2}


#FOURTH FUNCTION - CORRELATION COEFFICIENTS OF SIBLINGS 
#note that, because sibships differ in size, this isn't the proper way to figure intraclass correlations
#but I've checked it against ICC1.lme and there is little real difference
#NEW - actually, just use the pearson correlation rather than intraclass - huge improvement in time; little difference in value
#y <- effects.mate; effectnames <- PAR$parnames[5:length(PAR$parnames)]

multi.intra <- function(y,effectnames){
y <- t(y[c("Father.ID",effectnames),])

#ensures that there is a family with 9 sibs
newmat <- matrix(rep(c(1,rep(NA,length(effectnames))),9),byrow=TRUE,nrow=9)
dimnames(newmat)[2] <- dimnames(y)[2]
y <- as.data.frame(rbind(y,newmat))

#create matrix
y <- y[order(y$Father.ID),]
indx <- as.vector(table(y$Father.ID))
y$time <- unlist(sapply(indx,sample,replace=FALSE))
y1 <- as.matrix(reshape(y,idvar="Father.ID",direction="wide",v.names=effectnames,timevar="time"))
result <- numeric(length=length(effectnames))

for (j in 2:(length(effectnames)+1)){
y2 <- y1[,c(1,seq(j,ncol(y1),by=length(effectnames)))]  
y3 <- na.omit(rbind(y2[,c(2,3)],y2[,c(2,4)],y2[,c(2,5)],y2[,c(2,6)],y2[,c(2,7)],y2[,c(2,8)],y2[,c(2,9)],y2[,c(2,10)],y2[,c(3,4)],y2[,c(3,5)],y2[,c(3,6)],
            y2[,c(3,7)],y2[,c(3,8)],y2[,c(3,9)],y2[,c(3,10)],y2[,c(4,5)],y2[,c(4,6)],y2[,c(4,7)],y2[,c(4,8)],y2[,c(4,9)],y2[,c(4,10)],y2[,c(5,6)],y2[,c(5,7)],
            y2[,c(5,8)],y2[,c(5,9)],y2[,c(5,10)],y2[,c(6,7)],y2[,c(6,8)],y2[,c(6,9)],y2[,c(6,10)],y2[,c(7,8)],y2[,c(7,9)],y2[,c(7,10)],y2[,c(8,9)],
            y2[,c(8,10)],y2[,c(9,10)]))
result[j-1] <- round(cor(y3)[1,2],3)
}
result}



#FIFTH FUNCTION - CREATE MATRIX TO TRACK CHANGES OVER TIME
change.tracker <- function(run,x,y,mean.names,effects.names,beta1,var.names,cov.names,popsize,cor.spouses){
x[mean.names,run+1] <- round(apply(y[effects.names,],1,mean)*beta1,3)
cov.varcomp <- round(cov(t(y[effects.names,]*beta1)),3)
x[var.names,run+1] <- diag(cov.varcomp)
x[cov.names,run+1]  <- cov.varcomp[lower.tri(cov.varcomp)]
x["popsize",run+1] <- popsize
x["cor.spouses",run+1] <- round(cor.spouses,3)
x}




#SIXTH FUNCTION - FIGURE OUT EXPECTED VA CHANGE (DUE TO AM) OVER TIME
expectedVA <- function(z=track.changes,n=PAR$number.genes){
fe <- 1-(1/(2*n))
va0 <- z["var.A",1]
r <- z["cor.spouses",1]
number.gen <- which(dimnames(z)[[2]]=="reltp7")
va <- numeric(length=number.gen)
va[1] <- va0
ve <- sum(z[c("var.F","var.S","var.U","var.MZ","var.TW","var.SEX","var.AGE"),1])

for (i in 2:number.gen) {va[i] <- .5*va0 + (.5*va[i-1]) + (.5*r*fe)*(va[i-1]^2/(va[i-1]+ve))}
va.equilibrium <- va0/(1-(r*fe*(va[i-1]/(va[i-1]+ve))))

list(VA=va, VA.EQ=va.equilibrium)}




#SEVENTH FUNCTION - RANDOM PARAMETER VALUES;
#Use this to create random parameter values as specified by user. If rand.parameters =="no", this function doesn't do anything

alter.parameters <- function(P,V,max=1){
    if (P$rand.parameters != 'no') {
    xxx <- round(runif(1,0,max),3)

    if (P$rand.parameters =="A")  V$A <- xxx
    if (P$rand.parameters =="D")  V$D <- xxx
    if (P$rand.parameters =="F")  V$F <- xxx
    if (P$rand.parameters =="S")  V$S <- xxx
    
    if (P$rand.parameters =="all")  {V$A <- runif(1,0,.5);V$D <- runif(1,0,.5);V$F <- runif(1,0,.5);V$S <- runif(1,0,.5)
                   tot <- V$A+V$D+V$F+V$S; V$A <- round(V$A/tot,3);V$D <- round(V$D/tot,3);V$F <- round(V$F/tot,3); V$S <- round(V$S/tot,3)}

    if (P$rand.parameters =="AM.mod"){
    if (runif(1)<.5) P$am.model <- rep(1,14) else P$am.model <- c(0,0,0,1,1,1,1,1,1,1,0,0,0,0)}}

   return(BOTH=list(VAR=V,PAR=P))}

    
