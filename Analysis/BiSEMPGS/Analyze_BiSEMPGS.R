# a script to analyze BiSEMPGS data
library(ggplot2)
library(patchwork)

conditionNames <- c("phenotypic_transmission", "social_transmission")
sample_sizes <- c(4000, 8000, 16000)

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
  
    # exclude the outliers that is three standard deviations away from the mean, only applied to numerical variables
    #df <- df[apply(df[,1:(ncol(df)-15)], 2, function(x) all(abs(x - mean(x, na.rm = TRUE)) < 8*sd(x, na.rm = TRUE))),]

    # general lower bound variables
    varname_vector <- c()
    return(df)
}

sum_list1 <- readRDS("Analysis/BiSEMPGS/phenoVT_geneticAM_8000_summary_list.rds")

sum_df1 <- getDf(sum_list1, fixed = FALSE)

# write the df
write.csv(sum_df1, "Analysis/BiSEMPGS/phenoVT_geneticAM_8000_summary_list.csv", row.names = FALSE)

sum_list2 <- readRDS("Analysis/BiSEMPGS/phenoVT_phenoAM_8000_summary_list.rds")
sum_df2 <- getDf(sum_list2, fixed = FALSE)
# write the df
write.csv(sum_df2, "Analysis/BiSEMPGS/phenoVT_phenoAM_8000_summary_list.csv", row.names = FALSE)

sum_list3 <- readRDS("Analysis/BiSEMPGS/socialVT_phenoAM_8000_summary_list.rds")
sum_df3 <- getDf(sum_list3, fixed = FALSE)
# write the df
write.csv(sum_df3, "Analysis/BiSEMPGS/socialVT_phenoAM_8000_summary_list.csv", row.names = FALSE)

sum_list4 <- readRDS("Analysis/BiSEMPGS/phenoVT_socialAM_8000_summary_list.rds")
sum_df4 <- getDf(sum_list4, fixed = FALSE)
# write the df
write.csv(sum_df4, "Analysis/BiSEMPGS/phenoVT_socialAM_8000_summary_list.csv", row.names = FALSE)


# use psych::describe to get summary on each condition and save the summary table as tsv files
library(psych)
describe_df1 <- psych::describe(sum_df1)
#write.table(describe_df1, "Analysis/BiSEMPGS/phenoVT_geneticAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
describe_df2 <- psych::describe(sum_df2)
#write.table(describe_df2, "Analysis/BiSEMPGS/phenoVT_phenoAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
describe_df3 <- psych::describe(sum_df3)
#write.table(describe_df3, "Analysis/BiSEMPGS/socialVT_phenoAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
describe_df4 <- psych::describe(sum_df4)
#write.table(describe_df4, "Analysis/BiSEMPGS/phenoVT_socialAM_8000_summary_table.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

rownames(describe_df1)

describe_df1 <- as.data.frame(describe_df1)
describe_df1$variable <- rownames(describe_df1)
describe_df2 <- as.data.frame(describe_df2)
describe_df2$variable <- rownames(describe_df2)
describe_df3 <- as.data.frame(describe_df3)
describe_df3$variable <- rownames(describe_df3)
describe_df4 <- as.data.frame(describe_df4)
describe_df4$variable <- rownames(describe_df4)

# add a column for the condition
describe_df1$condition <- "phenoVT_geneticAM"
describe_df2$condition <- "phenoVT_phenoAM"
describe_df3$condition <- "socialVT_phenoAM"
describe_df4$condition <- "phenoVT_socialAM"

# combine the dataframes
describe_df <- rbind(describe_df1, describe_df2, describe_df3, describe_df4)

# put the condition column first and variable column second
describe_df <- describe_df[, c("condition", "variable", setdiff(names(describe_df), c("condition", "variable")))]

# write the describe_df to a tsv file
write.table(describe_df, "Analysis/BiSEMPGS/summary_table.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)