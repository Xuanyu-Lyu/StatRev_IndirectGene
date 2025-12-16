source("Analysis/BiSEMPGS/04-OpenMxFunctions.R")
library(crayon)
# conditionNames <- c("phenoVT_phenoAM",
#                     "socialVT_phenoAM",
#                     "phenoVT_socialAM",
#                     "phenoVT_geneticAM",
#                     "socialphenoVT_phenoAM",
#                     "t1pheVT_t2socVT_uniphenoAM",
#                     "01_t1pheVTnoAM_t2socVTnoAM",
#                     "02_t1noVTpheAM_t2noVTnoAM"
#                     #"03_t1noVTsocAM_t2noVTnoAM"
#                     #"04_t1noVTgenAM_t2noVTnoAM"
#                     )
#conditionNames <- c("01_t1pheVTnoAM_t2socVTnoAM", "02_t1noVTpheAM_t2noVTnoAM", "03_t1noVTsocAM_t2noVTnoAM", "04_t1noVTgenAM_t2noVTnoAM")
conditionNames <- c("05_t1pheVTnoAM_t2socVTnoAM_PGSall", "06_t1noVTpheAM_t2pheVTpheAM_PGSall", "07_t1noVTsocAM_t2pheVTsocAM_PGSall", "08_t1noVTgenAM_t2pheVTgenAM_PGSall")

v_sample <- c(3.2e4)
folder_pattern <- paste0("nfam", v_sample)
data_pattern <- paste0("_nfam", v_sample, ".txt")
save_pattern <- paste0("_", v_sample)

# folder_pattern <- c("nfam16000", "nfam32000", "nfam48000", "nfam64000", "nfam80000")
# data_pattern <- c("_16000.txt", "_32000.txt", "_48000.txt",  "_64000.txt", "_80000.txt")
# save_pattern <- c("_16000", "_32000", "_48000", "_64000", "_80000")

for (i in 1){
    for (j in 1:length(folder_pattern)){
        data_path <- paste0("Data/", conditionNames[i], "/", folder_pattern[j])
        l_files <- list.files(data_path, pattern = data_pattern[j])
        save_path <- paste0("Analysis/BiSEMPGS/", conditionNames[i], save_pattern[j], "_summary_list.rds")
        
        # Load existing results if they exist
        if (file.exists(save_path)){
            summary_list <- readRDS(save_path)
            cat(yellow("Loaded existing results for", conditionNames[i], "\t n=", folder_pattern[j], "\n"))
        } else {
            summary_list <- list()
        }
        
        for (k in c(184:length(l_files))){ # skip 261 and 708 which failed to converge
            fit <- fitBiSEMPGS_m2_tol(paste0(data_path, "/", l_files[k]), 
                                      feaTol = 1e-6, 
                                      optTol = 1e-9,
                                      jitterMean = 0.5,
                                      jitterVar = .1,
                                      exhaustive = TRUE,
                                      extraTries = 2)
            summary_list[[l_files[k]]] <- fit
            cat(magenta("\n", conditionNames[i], "\tn=", folder_pattern[j], "\tModel", l_files[k], "has been fitted\n"))
            
            # Save every 50 iterations
            if (k %% 50 == 0){
                saveRDS(summary_list, save_path)
                cat(cyan("\nIntermediate save at iteration", k, "for", conditionNames[i], "\t n=", folder_pattern[j], "\n"))
            }
        }
        
        # Final save
        saveRDS(summary_list, save_path)
        cat(green("\nFinal summary list for", conditionNames[i], "\t n=", folder_pattern[j], "has been saved\n"))

    }
}


# true_a <- list(c(0.438178, 0.464758),
#             c(0.5656854,0.464758),
#             c(0.669328,0.464758),
#             c(0.7589466,0.464758))

# true_a <- list(c(0.6928203, 0.5366563),
#             c(0.7483315,0.5366563),
#             c(0.7745967,0.5366563),
#             c(0.7874008,0.5366563),
#             c(0.7937254,0.5366563))

# if(TRUE){
#     for (i in 1:5){
#     for (j in 5){
#         data_path <- paste0("Data/Paper/", conditionNames[i], "/", folder_pattern[j])
#         l_files <- list.files(data_path, pattern = data_pattern[j])
#         # check if the data has been fitted
#         # if (file.exists(paste0("Analysis/Paper/", conditionNames[i], "/m2_paper_version2", save_pattern[j], "_summary_list_fixedA.rds"))){
#         #     cat("Summary list for", conditionNames[i], "\t n=", folder_pattern[j], "has been fitted\n")
#         #     next
#         # }
#         summary_list <- list()
#         for (k in 1:length(l_files)){
#             fit <- fitBiSEMPGS_m2_tol_fixH2(paste0(data_path, "/", l_files[k]), 
#                                       avalue = true_a[[i]],
#                                       feaTol = 1e-6, 
#                                       optTol = 1e-9,
#                                       jitterMean = 0.5,
#                                       jitterVar = .1,
#                                       exhaustive = FALSE,
#                                       extraTries = 5)
#             summary_list[[l_files[k]]] <- fit
#             cat(magenta("\n", conditionNames[i], "\tn=", folder_pattern[j], "\tModel", l_files[k], "has been fitted\n"))
#         }
#         save_path <- paste0("Analysis/Paper/", conditionNames[i], "/m2_paper_version2", save_pattern[j], "_summary_list_fixedA.rds")
#         saveRDS(summary_list, save_path)
#         cat(green("\nSummary list for", conditionNames[i], "\t n=", folder_pattern[j], "has been saved\n"))

#     }
# }
# }

