######################################################################################
# Compatible with GeneEvolve version 0.69 +                                          #
# SCRIPT FOR GENERATING GRAPHICS FROM PEDEVOLVE SIMULATIONS                          #
# by Matthew C Keller, updated 11/18/2007                                            #
# updated Nov 16, 2023 to get the rounding to 3 instead of 2 & fix a labels issue in plots
######################################################################################






############################################################################
#GENERATE.GRAPHICS
ge.graphics <- function(date.names,PARAMS,VARIANCES,PE.Var2,PE.Var3,rel.correlations,track.changes,data.info){

#Open PDF document
name <- paste("GE_Graph-",date.names,".pdf",sep="")
pdf(name,width=11, height=8.5,pointsize=12, paper='special')
############################################################################




############################################################################
#FOR DEBUGGING 
#If debugging, don't run "GENERATE GRAPHICS" and uncommment out the following
#date.names=ALL$TEMP$name.date
#PARAMS=ALL$PAR
#VARIANCES=VAR1
#PE.Var3=ALL$PE.Var3
#PE.Var2=ALL$PE.Var2
#rel.correlations=ALL$rel.correlations #need to find where this is created and round to 3 instead of 2
#track.changes=ALL$track.changes
#data.info=ALL$TEMP$info.widedata.updated

############################################################################






############################################################################
#0.5 GET INFORMATION ABOUT THE PHENOTYPE, VT, AND AM MODEL

bmat <- (c(VARIANCES$A,VARIANCES$AA,VARIANCES$D,VARIANCES$F,VARIANCES$S,VARIANCES$U,VARIANCES$MZ,VARIANCES$TW,VARIANCES$SEX,VARIANCES$AGE,VARIANCES$A.by.SEX,VARIANCES$A.by.AGE,VARIANCES$A.by.S,VARIANCES$A.by.U) != 0)*1

vt.paths <- bmat*PARAMS$vt.model
am.paths <- bmat*PARAMS$am.multiplier

ph.tf <- rep(FALSE,length(PARAMS$am.multiplier))
am.tf <- rep(FALSE,length(PARAMS$am.multiplier))
vt.tf <- rep(FALSE,length(PARAMS$vt.model))
ph.tf[bmat!=0] <- TRUE
am.tf[am.paths!=0] <- TRUE
vt.tf[vt.paths!=0] <- TRUE

ph.variables <- PARAMS$varnames[ph.tf]
am.variables <- PARAMS$varnames[am.tf]
vt.variables <- PARAMS$varnames[vt.tf]
############################################################################






############################################################################
#1 OUTPUT BASIC INFORMATION ABOUT RUN

line.nums <- 25
op <- par(mar=c(2,2,line.nums,1),xpd=TRUE,oma=c(1,1,1,1))

input <- round(unlist(VARIANCES)[1:PARAMS$number.parameters],3)
names(input)[names(input)=="mvt"] <- "m"
names(input)[names(input)=="pvt"] <- "p"
data.input <- c(data.info$n.mzm+data.info$n.mzf,data.info$n.dzm+data.info$n.dzf+data.info$n.dzos,data.info$n.data[2:5])

plot(1:14,1:14,type="n",axes=FALSE,bty='n',)
mtext(paste("Results Summary from GeneEvolve Simulation ",PARAMS$t4),line=line.nums-1,cex=1.5)

mtext("Basic Parameters:",line=line.nums-3,cex=1,at=1,adj=0)
mtext(paste("Output Dir:",as.character(PARAMS$output.directory)),line=line.nums-4,cex=.8,adj=0,at=1.25)
mtext(paste("Number of generations:",as.character(PARAMS$number.generations)),line=line.nums-5,cex=.8,adj=0,at=1.25)
mtext(paste("Population size at start:",as.character(PARAMS$start.pop.size)),line=line.nums-6,cex=.8,adj=0,at=1.25)
mtext(paste("Number of genes:",as.character(PARAMS$number.genes)),line=line.nums-7,cex=.8,adj=0,at=1.25)
mtext(paste("Components going into Phenotype:",sprintf(toString(ph.variables))),line=line.nums-8,cex=.8,adj=0,at=1.25)
mtext(paste("Components going into Parenting Phenotype:",sprintf(toString(vt.variables))),line=line.nums-9,cex=.8,adj=0,at=1.25)
mtext(paste("Components going into Mating Phenotype:",sprintf(toString(am.variables)),"; spousal corr on this phenotype = ",as.character(VARIANCES$latent.AM[1])),line=line.nums-10,cex=.8,adj=0,at=1.25)

mtext("Sample Sizes in ETFD Dataset:",line=line.nums-12,cex=1,at=1,adj=0)
mtext(c("MZ","DZ","Parents of Twins","Sibs of Twins","Spouses of Twins","Children of Twins"),line=line.nums-13, cex=.7,at=c(2,4,6,8,10,12),adj=0.5)
mtext(data.input,line=line.nums-14, cex=.7,at=c(2,4,6,8,10,12),adj=0.5)

mtext(c("Variance Components & Paths, User Input:","(Note: U+MZ+TW=E for non-twins; m & p = maternal & paternal VT paths)"),line=line.nums-16,cex=c(1,.7),at=c(1,6),adj=0)
mtext(names(input),line=line.nums-17, cex=.7,at=(1:14),adj=0.5)
mtext(input,line=line.nums-18, cex=.7,at=(1:14),adj=0.5)

mtext("Time",line=line.nums-20,cex=1,adj=0,at=1)
mtext(paste("Simulation started: ",as.character(PARAMS$t1)),line=line.nums-21,at=1.25,cex=.8,adj=0)
mtext(paste("Simulation ended: ",as.character(Sys.time())),line=line.nums-22,at=1.25,cex=.8,adj=0)
mtext(paste("Minutes taken - looping through generations: ",round(difftime(PARAMS$t3,PARAMS$t1,units='mins'),2)),line=line.nums-23,at=1.25,cex=.8,adj=0)
mtext(paste("Minutes taken - creating pedigree datasets: ",round(difftime(PARAMS$t4,PARAMS$t2,units='mins'),2)),line=line.nums-24,at=1.25,cex=.8,adj=0)
mtext(paste("Minutes taken - finding relative correlations: ",round(difftime(PARAMS$t5,PARAMS$t4,units='mins'),2)),line=line.nums-25,at=1.25,cex=.8,adj=0)
mtext(paste("Minutes taken - TOTAL: ",round(difftime(Sys.time(),PARAMS$t1,units='mins'),2)),line=line.nums-26,at=1.25,cex=.8,adj=0)

mtext("Warnings in script:",line=line.nums-28,cex=1,adj=0,at=1)
mtext(PARAMS$startpop.warn,line=line.nums-29,at=1.25,cex=.8,adj=0)


if (max(PE.Var3)>10) {
  mtext("The variance of a component is greater than 10x its starting value. Graphing will fail. Your starting values have created an unstable model.",
        line=line.nums-33,at=1.25,cex=.8,adj=0)}

############################################################################






###############################################################################
#2 GRAPHING PEDIGREE CORRELATIONS
par(mar=c(3,2,6,1),xpd=TRUE,oma=c(2,2,2,2))

locs <-rel.correlations + .03
locs[locs<0] <- .05

mycolors <- c('darkblue','mediumpurple','turquoise','goldenrod4','lightgoldenrod',
              'brown','orangered','darkgreen','darkseagreen','snow2')
       
bpl <- barplot(rel.correlations, col=mycolors,ylim=c(-.1,1), xlab="Relative Type",ylab="Correlation",
        names.arg=c("MZ","DZ","Sib","Par.Cd","Gp.Gcd","MZAv","Avnc","MZcs","Cs","Sps"),cex.names=1)
text(bpl,locs,round(rel.correlations,3))

varcovar <- c(PE.Var2[c(17:30,33)],PE.Var2[31]*track.changes["cor.spouses",ncol(track.changes)])/PE.Var2[31]
lim <- par("usr")
labels.at <- seq(lim[1],lim[2],length.out=length(varcovar))

mtext("Correlations between 10 Relative Types",line=5,cex=1.5)
mtext("True Standardized Variance Components in Dataset",line=3,cex=.75)
mtext(c(names(input),"Covs","AM"),line=2, cex=.7,at=labels.at,adj=0.5)
mtext(round(varcovar,3),line=1, cex=.7,at=labels.at,adj=0.5)
################################################################################




try({



############################################################################
#3 GRAPHING TRACK CHANGES

#Axis break stuff
ymax <- max(ceiling(PE.Var3[,1:(ncol(PE.Var3)-1)]*11))/10

#make new track.changes that includes "E"
par(mar=c(4,4,8,1),oma=c(1,1,1,1),xpd=NA)
names2 <- rownames(PE.Var3)
start <- PE.Var3[names2,1]
end <- PE.Var3[names2,(ncol(PE.Var3)-1)]
dat <- PE.Var3[names2,ncol(PE.Var3)]
names3 <- names2[1:nrow(PE.Var3)-1]        #all names except r(sps) - meaning that we won't graph r(sps)
names3 <- names3[PE.Var3[1:(nrow(PE.Var3)-1),ncol(PE.Var3)] != 0]  
use.colors <- c(c('darkblue','red','green','orange','yellow','purple','brown','cyan4','olivedrab4','hotpink1','lightblue')[1:(length(names3)-1)],"black")

#Plot it
plot(PE.Var3["V(P)",1:(ncol(PE.Var3)-1)],ylim=c(0,ymax),type='n',lwd=2,col='black',
     xlab='Generation Number',ylab='Variance or Covariance',bty='n',axes=TRUE)
lim <- par("usr")
labels.at <- seq(lim[1],lim[2],length.out=length(names2)+1)

mtext("Change in Variance Across Generations - Includes V(P)",line=6,cex=1.5)
mtext("Unstandardized Variance Components",line=5,cex=.7)
mtext(names2,side=3,outer=FALSE,line=4, cex=.7,at=labels.at[2:length(labels.at)])
mtext(c("Start:",round(start,3)),side=3,outer=FALSE,line=3, cex=.7,at=labels.at[1:length(labels.at)])
mtext(c("End:",round(end,3)),side=3,outer=FALSE,line=2, cex=.7,at=labels.at[1:length(labels.at)])
mtext(c("Data:",round(dat,3)),side=3,outer=FALSE,line=1, cex=.7,at=labels.at[1:length(labels.at)])


#legend(x=c(lim[1],(lim[2]-.1*(lim[2]-lim[1]))),
#       y=c(lim[4],(lim[4]-.15*(lim[4]-lim[3]))),
#       legend=names3,lty=1,col=use.colors,lwd=2,ncol=min(5,length(names3)),bg='snow2',yjust=1,text.width=1)

legend(x=lim[1], y=lim[4],legend=names3,lty=1,col=use.colors,lwd=2,ncol=min(5,length(names3)),bg='snow2',yjust=1,text.width=1)


for (i in 1:length(names3)){
lines(PE.Var3[names3[i],1:(ncol(PE.Var3)-1)],col=use.colors[i],lwd=2)}
###############################################################################






############################################################################
#4 GRAPHING TRACK CHANGES

PE.Var4 <- PE.Var3[- (nrow(PE.Var3)-1),]      #remove V(P), which is next to last row

#Axis break stuff
ymax <- max(ceiling(PE.Var4[,1:(ncol(PE.Var4)-1)]*11))/10

#make new track.changes that includes "E"
par(mar=c(4,4,8,1),oma=c(1,1,1,1),xpd=NA)
names5 <- rownames(PE.Var4)
start2 <- PE.Var4[names5,1]
end2 <- PE.Var4[names5,(ncol(PE.Var4)-1)]
dat2 <- PE.Var4[names5,ncol(PE.Var4)]
names6 <- names5[1:(nrow(PE.Var4)-1)]    #remove r(sps)
names6 <- names6[PE.Var4[1:(nrow(PE.Var4)-1),ncol(PE.Var4)] != 0]
use.colors <- c('darkblue','red','green','orange','yellow','purple','brown','cyan4','olivedrab4','hotpink1','lightblue')[1:(length(names6))]

#Plot it
plot(PE.Var4["V(A)",1:(ncol(PE.Var4)-1)],ylim=c(0,ymax),type='n',lwd=2,col='black',
     xlab='Generation Number',ylab='Variance or Covariance',bty='n',axes=TRUE)
lim <- par("usr")
labels.at2 <- seq(lim[1],lim[2],length.out=length(input)+1)

mtext("Change in Variance Across Generations - Does Not Include V(P)",line=6,cex=1.5)
mtext("Unstandardized Variance Components",line=5,cex=.7)
mtext(names5,side=3,outer=FALSE,line=4, cex=.7,at=labels.at[2:length(labels.at2)])
mtext(c("Start:",round(start2,3)),side=3,outer=FALSE,line=3, cex=.7,at=labels.at[1:length(labels.at2)])
mtext(c("End:",round(end2,3)),side=3,outer=FALSE,line=2, cex=.7,at=labels.at[1:length(labels.at2)])
mtext(c("Data:",round(dat2,3)),side=3,outer=FALSE,line=1, cex=.7,at=labels.at[1:length(labels.at2)])

#legend(x=c(lim[1],(lim[2]-.1*(lim[2]-lim[1]))),
#       y=c(lim[4],(lim[4]-.15*(lim[4]-lim[3]))),
#       legend=names6,lty=1,col=use.colors,lwd=2,ncol=min(5,length(names6)),bg='snow2',yjust=1,text.width=1)

legend(x=lim[1],y=lim[4],legend=names6,lty=1,col=use.colors,lwd=2,ncol=min(5,length(names6)),bg='snow2',yjust=1,text.width=1)

for (i in 1:length(names6)){
lines(PE.Var4[names6[i],1:(ncol(PE.Var4)-1)],col=use.colors[i],lwd=2)}
###############################################################################







############################################################################
#5 GRAPH CHANGES IN VA DUE TO AM

Plot.expectedVA <- function(y=PE.Var3, expected,input=TEMP$input){
  
  #create the three variance components
  va <- y["V(A)",1:length(expected$VA)]
  eva <- expected$VA
  ve <- y["V(P)",1:length(expected$VA)] - va
  
  #plotting parameters
  ymax <- max(ceiling(y[1:(nrow(y)-2),1:(ncol(y)-1)]*11))/10
  par(mar=c(4,4,7,1),oma=c(1,1,1,1),xpd=NA)
  names2 <- rownames(y)
  start <- y[names2,1]
  end <- y[names2,(ncol(y)-1)]
  names3 <- names2[y[,ncol(y)] != 0]
  names3 <- c("V(A)","E[V(A)]")
  use.colors <- c('darkblue','red')
  
  #Plot it
  plot(va,ylim=c(0,ymax),type='n',lwd=2,col='darkblue',xlab='Generation Number',ylab='Variance',bty='n',axes=TRUE)
  lim <- par("usr")
  labels.at <- seq(lim[1],lim[2],length.out=length(input)+2)
  
  mtext("Change in V(A) Across Generations due to Assortative Mating",line=5,cex=1.5)
  mtext("Unstandardized Variance Components",line=4,cex=.75)
  mtext(names2,side=3,outer=FALSE,line=3, cex=.75,at=labels.at[2:length(labels.at)])
  mtext(c("Start:",round(start,3)),side=3,outer=FALSE,line=2, cex=.75,at=labels.at[1:length(labels.at)])
  mtext(c("End:",round(end,3)),side=3,outer=FALSE,line=1, cex=.75,at=labels.at[1:length(labels.at)])
  
#  legend(x=c(lim[1],(lim[2]-.1*(lim[2]-lim[1]))),
#         y=c(lim[4],(lim[4]-.1*(lim[4]-lim[3]))),
#         legend=names3,lty=1,col=use.colors,lwd=2,ncol=min(5,length(names3)),bg='snow2',yjust=1,text.width=1)
  
  legend(x=lim[1],y=lim[4],legend=names3,lty=1,col=use.colors,lwd=2,ncol=min(5,length(names3)),bg='snow2',yjust=1,text.width=1)
  
  lines(va,col=use.colors[1],lwd=2)
  lines(eva,col=use.colors[2],lwd=2,lty=2)
  #lines(ve,col=use.colors[3],lwd=2)
}


ex <- expectedVA(z=track.changes,n=PARAMS$number.genes)
Plot.expectedVA(y=PE.Var3, expected=ex,input=input)



############################################################################




})




########################
#END GRAPHICS SCRIPT

par(op)

dev.off()

}
########################








