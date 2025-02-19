# Step 1: simulate a group of data
source("Scripts/GeneEvolve/SIMULATE_DAT_GEN.R")

# Initialize the genetic information
num.cvs = 100

RUN.MARKERS <- FALSE #whether to only consider GRMs built from CVs (FALSE) or both CVs and SNPs (TRUE)

MIN.MAF <- .1; MAX.MAF <- .50#AM Simulation wildcards

MAF.VECTOR <- runif(num.cvs,MIN.MAF,MAX.MAF) #Can change the distribution of MAFs here

GENTP.VAR <- MAF.VECTOR*(1-MAF.VECTOR)*2

ALPHA.VECTOR <- sample(c(-1,1),num.cvs,replace=TRUE)*sqrt(1/(num.cvs*GENTP.VAR)) #Can change the distribution of effect sizes here - fixed f'n of MAF

CV.INFO <- data.frame(MAF=MAF.VECTOR,alpha=ALPHA.VECTOR)

results <- list()
for(k in 1:10){
    data_list <- AM.SIMULATE(
        CV.INFO = CV.INFO, 
        H2.T0 = .5, 
        NUM.GENERATIONS = 20, 
        POP.SIZE = 5000, 
        MATE.COR = .4, 
        AVOID.INB = TRUE, 
        SAVE.EACH.GEN = TRUE, 
        SAVE.COVS = TRUE, 
        SEED = k*10, 
        VF.T0 = 0, 
        PROP.H2.LATENT = .7, 
        Unequal_AM = FALSE)

    data_df <- data_list$HISTORY 
    results[[k]] <- data_df
    cat("Simulation ", k, " completed\n")
    
}
#saveRDS(results, file = "Scripts/PythonRDR/Her.5-AM.4-Lat.5-VF0-results.rds")
as.data.frame(results[[1]]["XO"]$XO[20]) |> head()
as.data.frame(results[[1]]["XL"]$XL[20]) |> head()
as.data.frame(results[[1]]$PHEN[[20]]) |> names()

# get the phenotype and genetic information for the offspring
df_phe_gen_offspring <- cbind(as.data.frame(results[[1]]$PHEN[[20]])[,c("ID","Y", "Father.ID","Mother.ID")], as.data.frame(results[[1]]["XO"]$XO[20]))
colnames(df_phe_gen)[1:4] <- c("ID","Y", "Father.ID","Mother.ID")

# get the genetic information for the parents
df_gen_parent <- cbind(as.data.frame(results[[1]]$PHEN[[19]])[,c("ID")], as.data.frame(results[[1]]["XO"]$XO[19]))
colnames(df_gen_parent)[1] <- "ID"

# get the genetic information for the parents, based on offspring row
for(i in 1: nrow(df_phe_gen_offspring)){
    mother_gene <- df_gen_parent[df_gen_parent$ID == df_phe_gen_offspring$Mother.ID[i],][,-1]
    father_gene <- df_gen_parent[df_gen_parent$ID == df_phe_gen_offspring$Father.ID[i],][,-1]
    if(i == 1){
        df_gene_mother <- mother_gene
        df_gene_father <- father_gene
    } else {
        df_gene_mother <- rbind(df_gene_mother, mother_gene)
        df_gene_father <- rbind(df_gene_father, father_gene)
    }
}

# save offspring, mother and father genetic information
write.table(df_phe_gen_offspring, file = "Scripts/PythonRDR/Her.5-AM.4-Lat.5-VF0-offspring.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(df_gene_mother, file = "Scripts/PythonRDR/Her.5-AM.4-Lat.5-VF0-mother.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(df_gene_father, file = "Scripts/PythonRDR/Her.5-AM.4-Lat.5-VF0-father.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

