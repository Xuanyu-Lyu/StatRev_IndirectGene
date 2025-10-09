# BiSEMPGS Data Extraction Script
# This script reads RDS files containing BiSEMPGS results and extracts parameter estimates
# into CSV files for further processing

library(psych)

# Function to extract parameter estimates from BiSEMPGS summary list
getDf <- function(summary_list, condition_name, fixed = FALSE) {
    cat("Processing condition:", condition_name, "\n")
    
    # Get status codes for each run
    status_codes <- sapply(summary_list, function(x) x$statusCode)
    
    # Initialize dataframe with parameter columns
    df <- data.frame(matrix(ncol = nrow(summary_list[[2]]$parameters), nrow = length(summary_list)))
    colnames(df) <- summary_list[[2]]$parameters$name
    
    # Extract parameter estimates
    for(i in 1:length(summary_list)) {
        for(j in 1:nrow(summary_list[[i]]$parameters)){
            if (!is.null(summary_list[[i]]$parameters$Estimate[j])) {
                df[i,j] <- summary_list[[i]]$parameters$Estimate[j]
            } else {
                cat("NULL value at run", i, "parameter", j, "\n")
            }
        }
    }
    
    # Add status codes and filter for successful runs
    df$status_codes <- status_codes
    df <- df[df$status_codes %in% c("OK", "OK/green"),]
    
    # Add condition identifier
    df$condition <- condition_name
    
    # Ensure required parameter columns exist
    required_cols <- c('f11','f12','f21','f22',
                       'delta11','delta22',
                       'a11','a22',
                       'k11','k12','k22')
    
    for(col in required_cols){
        if(!(col %in% names(df))){
            df[[col]] <- NA_real_
        } else {
            # Coerce to numeric to avoid factor/list issues
            df[[col]] <- as.numeric(df[[col]])
        }
    }
    
    # Add replication number
    df$replication <- 1:nrow(df)
    
    cat("Successfully processed", nrow(df), "runs for", condition_name, "\n")
    return(df)
}

# Main execution
main <- function() {
    cat("Starting BiSEMPGS parameter extraction...\n")
    
    # Define conditions and their corresponding RDS files (conditions 5-8)
    conditions <- list(
        #"05_t1pheVTnoAM_t2socVTnoAM_PGSall" = "Analysis/BiSEMPGS/05_t1pheVTnoAM_t2socVTnoAM_PGSall_8000_summary_list.rds",
        "06_t1noVTpheAM_t2pheVTpheAM_PGSall" = "Analysis/BiSEMPGS/06_t1noVTpheAM_t2pheVTpheAM_PGSall_8000_summary_list.rds",
        "07_t1noVTsocAM_t2pheVTsocAM_PGSall" = "Analysis/BiSEMPGS/07_t1noVTsocAM_t2pheVTsocAM_PGSall_8000_summary_list.rds",
        "08_t1noVTgenAM_t2pheVTgenAM_PGSall" = "Analysis/BiSEMPGS/08_t1noVTgenAM_t2pheVTgenAM_PGSall_8000_summary_list.rds"
    )
    
    # Process each condition
    all_dataframes <- list()
    
    for(condition_name in names(conditions)) {
        rds_file <- conditions[[condition_name]]
        
        cat("\n", strrep("=", 60), "\n")
        cat("Processing:", condition_name, "\n")
        cat("Reading:", rds_file, "\n")
        
        # Check if file exists
        if(!file.exists(rds_file)) {
            cat("WARNING: File not found:", rds_file, "\n")
            next
        }
        
        # Read RDS file
        tryCatch({
            summary_list <- readRDS(rds_file)
            
            # Extract parameters
            condition_df <- getDf(summary_list, condition_name, fixed = FALSE)
            
            # Save individual condition CSV
            output_file <- paste0("Analysis/BiSEMPGS/", condition_name, "_parameters.csv")
            write.csv(condition_df, output_file, row.names = FALSE)
            cat("Saved parameters to:", output_file, "\n")
            
            # Store for combined dataframe
            all_dataframes[[condition_name]] <- condition_df
            
        }, error = function(e) {
            cat("ERROR processing", condition_name, ":", e$message, "\n")
        })
    }
    
    # Combine all conditions into one dataframe
    if(length(all_dataframes) > 0) {
        cat("\nCombining all conditions...\n")
        combined_df <- do.call(rbind, all_dataframes)
        
        # Save combined CSV
        combined_output <- "Analysis/BiSEMPGS/all_conditions_parameters.csv"
        write.csv(combined_df, combined_output, row.names = FALSE)
        cat("Saved combined parameters to:", combined_output, "\n")
        
        # Create and save summary statistics
        cat("Creating summary statistics...\n")
        
        # Remove non-numeric columns for summary
        numeric_cols <- sapply(combined_df, is.numeric)
        numeric_df <- combined_df[, numeric_cols]
        
        # Create summary by condition
        summary_stats <- list()
        for(condition in unique(combined_df$condition)) {
            condition_data <- combined_df[combined_df$condition == condition, ]
            condition_numeric <- condition_data[, numeric_cols]
            
            if(nrow(condition_numeric) > 0) {
                describe_result <- psych::describe(condition_numeric)
                describe_df <- as.data.frame(describe_result)
                describe_df$variable <- rownames(describe_df)
                describe_df$condition <- condition
                
                # Reorder columns
                describe_df <- describe_df[, c("condition", "variable", 
                                             setdiff(names(describe_df), c("condition", "variable")))]
                
                summary_stats[[condition]] <- describe_df
            }
        }
        
        # Combine summary statistics
        if(length(summary_stats) > 0) {
            combined_summary <- do.call(rbind, summary_stats)
            
            # Save summary statistics
            summary_output <- "Analysis/BiSEMPGS/parameters_summary_statistics.tsv"
            write.table(combined_summary, summary_output, sep = "\t", 
                       row.names = FALSE, col.names = TRUE, quote = FALSE)
            cat("Saved summary statistics to:", summary_output, "\n")
        }
        
        #cat("\n", "="*60, "\n")
        cat("Parameter extraction complete!\n")
        cat("Total conditions processed:", length(all_dataframes), "\n")
        cat("Total runs across all conditions:", nrow(combined_df), "\n")
        #cat("="*60, "\n")
        
    } else {
        cat("No conditions were successfully processed.\n")
    }
}

# Run the main function
main()