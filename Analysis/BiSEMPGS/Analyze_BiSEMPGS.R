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

sum_list <- readRDS("Analysis/BiSEMPGS/phenotypic_transmission_16000_summary_list.rds")

df_phenoVT_16000 <- getDf(sum_list, fixed = FALSE)

psych::describe(df_phenoVT_16000)


0.54^2 + 0.71^2

0.44^2 + 0.72^2
