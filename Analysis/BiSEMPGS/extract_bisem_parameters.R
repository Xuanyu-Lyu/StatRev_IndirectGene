# BiSEMPGS Data Extraction Script
# This script reads RDS files containing BiSEMPGS results and extracts parameter estimates
# into CSV files for further processing

library(psych)

# Function to extract parameter estimates from BiSEMPGS summary list
getDf <- function(summary_list, condition_name, fixed = FALSE) {
    cat("Processing condition:", condition_name, "\n")
    
    # Debug: Check the structure of summary_list
    cat("Summary list class:", class(summary_list), "\n")
    cat("Summary list length:", length(summary_list), "\n")
    
    # If it's not a list or has unexpected structure, print first few elements
    if(!is.list(summary_list) || length(summary_list) == 0) {
        cat("Unexpected data structure. First few elements:\n")
        print(head(summary_list))
        stop("summary_list is not in expected format (list of result objects)")
    }
    
    # Check first element structure
    if(length(summary_list) > 0) {
        cat("First element class:", class(summary_list[[1]]), "\n")
        if(is.list(summary_list[[1]])) {
            cat("First element names:", names(summary_list[[1]]), "\n")
        }
    }
    
    # Validate that we have the expected structure
    if(!is.list(summary_list[[1]]) || !("statusCode" %in% names(summary_list[[1]]))) {
        stop("First element doesn't have expected structure (missing statusCode)")
    }
    
    if(!("parameters" %in% names(summary_list[[1]]))) {
        stop("First element doesn't have expected structure (missing parameters)")
    }
    
    # Get status codes for each run
    status_codes <- sapply(summary_list, function(x) {
        if(is.list(x) && "statusCode" %in% names(x)) {
            return(x$statusCode)
        } else {
            return("UNKNOWN")
        }
    })
    
    # Find a valid parameters structure for column names
    valid_params_index <- which(sapply(summary_list, function(x) {
        is.list(x) && "parameters" %in% names(x) && 
        is.data.frame(x$parameters) && nrow(x$parameters) > 0
    }))[1]
    
    if(is.na(valid_params_index)) {
        stop("No valid parameters found in any result")
    }
    
    cat("Using parameters structure from index:", valid_params_index, "\n")
    
    # Initialize dataframe with parameter columns
    df <- data.frame(matrix(ncol = nrow(summary_list[[valid_params_index]]$parameters), nrow = length(summary_list)))
    colnames(df) <- summary_list[[valid_params_index]]$parameters$name
    
    # Extract parameter estimates
    for(i in 1:length(summary_list)) {
        if(is.list(summary_list[[i]]) && "parameters" %in% names(summary_list[[i]])) {
            params <- summary_list[[i]]$parameters
            if(is.data.frame(params) && nrow(params) > 0) {
                for(j in 1:nrow(params)){
                    if (!is.null(params$Estimate[j])) {
                        df[i,j] <- params$Estimate[j]
                    } else {
                        cat("NULL value at run", i, "parameter", j, "\n")
                    }
                }
            }
        } else {
            cat("Invalid structure at run", i, "\n")
        }
    }
    
    # Add status codes and filter for successful runs
    df$status_codes <- status_codes
    
    # Check what status codes are present
    cat("Status codes found:", unique(status_codes), "\n")
    cat("Status code counts:\n")
    print(table(status_codes))
    
    # Filter for successful runs (status code 1 typically indicates success in OpenMx)
    df <- df[df$status_codes %in% c(0, 1, 2) | df$status_codes %in% c("OK", "OK/green"),]
    
    # Check if we have any successful runs
    if(nrow(df) == 0) {
        cat("WARNING: No successful runs found for", condition_name, "\n")
        cat("Available status codes:", unique(status_codes), "\n")
        # Return empty dataframe with proper structure
        df <- data.frame(matrix(ncol = ncol(df) + 1, nrow = 0))
        colnames(df) <- c(colnames(df)[1:(ncol(df)-1)], "condition")
        return(df)
    }
    
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
        "05_t1pheVTnoAM_t2socVTnoAM_PGSall" = "Analysis/BiSEMPGS/05_t1pheVTnoAM_t2socVTnoAM_PGSall_8000_summary_list.rds",
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
            
            # Debug: Print basic info about loaded data
            cat("Loaded RDS file successfully\n")
            cat("Data class:", class(summary_list), "\n")
            cat("Data length/size:", length(summary_list), "\n")
            
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
            cat("Error details:\n")
            print(e)
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