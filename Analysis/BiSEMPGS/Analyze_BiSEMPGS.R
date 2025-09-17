# a script to analyze BiSEMPGS data
library(ggplot2)
library(patchwork)

sample_sizes <- c(8000)

getDf <- function(summary_list, fixed = FALSE) {
    status_codes <- sapply(summary_list, function(x) x$statusCode)
    df <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
    colnames(df) <- summary_list[[2]]$parameters$name
    #colnames(df) 
    # Loop over the elements in the summary_list
    for(i in 1:length(summary_list)) {
        for(j in 1:nrow(summary_list[[i]]$parameters)){
            if (!is.null(summary_list[[i]]$parameters$Estimate[j])) {
                df[i,j] <- summary_list[[i]]$parameters$Estimate[j]
            } else {
                print(paste("NULL value at i =", i, "and j =", j))
            }
        }
    }

    df$status_codes <- status_codes
    df <- df[df$status_codes %in% c("OK", "OK/green"),]
    # exclude a that hit the lower bound
    if(!fixed){
        #df <- df[df$a11!=0.3 & df$a22!=0.3,]
        #df <- df[df$a11>0.4 & df$a22>0.4,]
        #df <- df[df$VY11>1 & df$VY22>1,]
       # df <- df[df$VE11>0 & df$VE12>0 & df$VE22>0,]
    }
        # Ensure required parameter columns exist (avoid NULL columns -> subscript out of bounds)
    required_cols <- c('f11','f12','f21','f22',
                       'delta11','delta22',
                       'a11','a22',
                       'k11','k12','k22')
    for(col in required_cols){
        if(!(col %in% names(df))){
            df[[col]] <- NA_real_
        } else {
            # coerce to numeric to avoid factor/list issues
            df[[col]] <- as.numeric(df[[col]])
        }
    }
    # Initialize new columns for genetic nurture effects
    df$phi11 <- NA
    df$phi12 <- NA
    df$phi21 <- NA
    df$phi22 <- NA
    df$rho11 <- NA
    df$rho12 <- NA
    df$rho21 <- NA
    df$rho22 <- NA
    
    # Compute genetic nurture effects for each row
    for (i in 1:nrow(df)) {
        # Create matrices from the estimates
        f_mat <- matrix(c(df$f11[i], df$f12[i], df$f21[i], df$f22[i]), nrow=2, byrow=TRUE)
        delta_mat <- matrix(c(df$delta11[i], 0, 0, df$delta22[i]), nrow=2, byrow=TRUE)
        a_mat <- matrix(c(df$a11[i], 0, 0, df$a22[i]), nrow=2, byrow=TRUE)
        k_mat <- matrix(c(.5, df$k12[i], df$k12[i], .5), nrow=2, byrow=TRUE)
        j_mat = k_mat
        # Compute pure observed genetic nurture effects: phi = 2*f*delta*k
        phi_mat <- f_mat %*% delta_mat %*% k_mat %*% t(delta_mat) %*% t(f_mat) *4
        df$phi11[i] <- phi_mat[1,1]
        df$phi12[i] <- phi_mat[1,2]
        df$phi21[i] <- phi_mat[2,1]
        df$phi22[i] <- phi_mat[2,2]
        
        # Compute pure latent genetic nurture effects: rho = 2*f*a*j
        rho_mat <- f_mat %*% a_mat %*% j_mat %*% t(a_mat) %*% t(f_mat) *4
        df$rho11[i] <- rho_mat[1,1]
        df$rho12[i] <- rho_mat[1,2]
        df$rho21[i] <- rho_mat[2,1]
        df$rho22[i] <- rho_mat[2,2]
    }
    
    # exclude the outliers that is three standard deviations away from the mean, only applied to numerical variables
    #df <- df[apply(df[,1:(ncol(df)-15)], 2, function(x) all(abs(x - mean(x, na.rm = TRUE)) < 8*sd(x, na.rm = TRUE))),]

    # general lower bound variables
    varname_vector <- c()
    return(df)
}

# sum_list1 <- readRDS("Analysis/BiSEMPGS/phenoVT_geneticAM_8000_summary_list.rds")

# sum_df1 <- getDf(sum_list1, fixed = FALSE)

# # write the df
# write.csv(sum_df1, "Analysis/BiSEMPGS/phenoVT_geneticAM_8000_summary_list.csv", row.names = FALSE)

# sum_list2 <- readRDS("Analysis/BiSEMPGS/phenoVT_phenoAM_8000_summary_list.rds")
# sum_df2 <- getDf(sum_list2, fixed = FALSE)
# # write the df
# write.csv(sum_df2, "Analysis/BiSEMPGS/phenoVT_phenoAM_8000_summary_list.csv", row.names = FALSE)

# sum_list3 <- readRDS("Analysis/BiSEMPGS/socialVT_phenoAM_8000_summary_list.rds")
# sum_df3 <- getDf(sum_list3, fixed = FALSE)
# # write the df
# write.csv(sum_df3, "Analysis/BiSEMPGS/socialVT_phenoAM_8000_summary_list.csv", row.names = FALSE)

# sum_list4 <- readRDS("Analysis/BiSEMPGS/phenoVT_socialAM_8000_summary_list.rds")
# sum_df4 <- getDf(sum_list4, fixed = FALSE)
# # write the df
# write.csv(sum_df4, "Analysis/BiSEMPGS/phenoVT_socialAM_8000_summary_list.csv", row.names = FALSE)

sum_list5 <- readRDS("Analysis/BiSEMPGS/t1pheVT_t2socVT_uniphenoAM_8000_summary_list.rds")
sum_df5 <- getDf(sum_list5, fixed = FALSE)
# write the df
#write.csv(sum_df5, "Analysis/BiSEMPGS/t1pheVT_t2socVT_uniphenoAM_8000_summary_list.csv", row.names = FALSE)


sum_list6 <- readRDS("Analysis/BiSEMPGS/01_t1pheVTnoAM_t2socVTnoAM_8000_summary_list.rds")
sum_df6 <- getDf(sum_list6, fixed = FALSE)

sum_list7 <- readRDS("Analysis/BiSEMPGS/02_t1noVTpheAM_t2noVTnoAM_8000_summary_list.rds")
sum_df7 <- getDf(sum_list7, fixed = FALSE)

# # use psych::describe to get summary on each condition and save the summary table as tsv files
# library(psych)
# describe_df1 <- psych::describe(sum_df1)
# #write.table(describe_df1, "Analysis/BiSEMPGS/phenoVT_geneticAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
# describe_df2 <- psych::describe(sum_df2)
# #write.table(describe_df2, "Analysis/BiSEMPGS/phenoVT_phenoAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
# describe_df3 <- psych::describe(sum_df3)
# #write.table(describe_df3, "Analysis/BiSEMPGS/socialVT_phenoAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
# describe_df4 <- psych::describe(sum_df4)
#write.table(describe_df4, "Analysis/BiSEMPGS/phenoVT_socialAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
describe_df5 <- psych::describe(sum_df5)

describe_df6 <- psych::describe(sum_df6)

describe_df7 <- psych::describe(sum_df7)

rownames(describe_df5)

# describe_df1 <- as.data.frame(describe_df1)
# describe_df1$variable <- rownames(describe_df1)
# describe_df2 <- as.data.frame(describe_df2)
# describe_df2$variable <- rownames(describe_df2)
# describe_df3 <- as.data.frame(describe_df3)
# describe_df3$variable <- rownames(describe_df3)
# describe_df4 <- as.data.frame(describe_df4)
# describe_df4$variable <- rownames(describe_df4)
describe_df5 <- as.data.frame(describe_df5)
describe_df5$variable <- rownames(describe_df5)
describe_df6 <- as.data.frame(describe_df6)
describe_df6$variable <- rownames(describe_df6)
describe_df7 <- as.data.frame(describe_df7)
describe_df7$variable <- rownames(describe_df7)

# # add a column for the condition
# describe_df1$condition <- "phenoVT_geneticAM"
# describe_df2$condition <- "phenoVT_phenoAM"
# describe_df3$condition <- "socialVT_phenoAM"
# describe_df4$condition <- "phenoVT_socialAM"
describe_df5$condition <- "t1pheVT_t2socVT_uniphenoAM"
describe_df6$condition <- "01_t1pheVTnoAM_t2socVTnoAM"
describe_df7$condition <- "02_t1noVTpheAM_t2noVTnooAM"

# combine the dataframes
describe_df <- rbind(describe_df5, describe_df6, describe_df7)

# put the condition column first and variable column second
describe_df <- describe_df[, c("condition", "variable", setdiff(names(describe_df), c("condition", "variable")))]

# write the describe_df to a tsv file
write.table(describe_df, "Analysis/BiSEMPGS/summary_table.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

