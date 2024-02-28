library(rlang)
library(withr)
library(rex)
library(dplyr)
library(visNetwork)
library(Matrix)
library(combinat)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename
matches <- dplyr::matches
mutate <- dplyr::mutate
desc <- dplyr::desc
library(tidyr)
library(magrittr)
subtract <- magrittr::subtract
library(stringr)
library(glue)
library(igraph)
library(openxlsx)
library(reshape2)
as_data_frame <- igraph::as_data_frame
library(assertthat)
library(limma)
library(this.path)
library(entropy)
column_to_rownames <- tibble::column_to_rownames
rownames_to_column <- tibble::rownames_to_column
library(gtools)
library(tools)
library(peakRAM)

script_dir <- this.dir()
source(file.path(script_dir, "sim_mislabeled_dataset.R"))
source(file.path(script_dir, "MislabelSolver.R"))

cmd_args <- commandArgs(trailingOnly = TRUE)
params_grid_file <- cmd_args[1]

# Column of relabels for each incremental improvement
# 1. unambiguous majority
# 2. majority with cycle detection
# 3. majority with comprehensivesim_name
# 4. iterative ensemble with local search

local({
    n_subjects = 200
    n_samples_per_subject = 10
    n_swap_cats = 5
    fraction_mislabel = 0.5
    fraction_anchor = 0.1
    fraction_ghost = 0.1
    seed = 1986
    output_dir = "/Users/charlesdeng/Workspace/mislabeling/simulations"
    sim_name = "testtest"
    args_list = list(n_subjects = n_subjects, 
                     n_samples_per_subject = n_samples_per_subject, 
                     n_swap_cats = n_swap_cats, 
                     fraction_mislabel = fraction_mislabel, 
                     fraction_anchor = fraction_anchor, 
                     fraction_ghost = fraction_ghost, 
                     seed = seed, 
                     output_dir = output_dir, 
                     sim_name = sim_name)
})

## all is 5352.4
## to 118 is 1982 (swapping samples)
## to 149 is 1094.1 (baseline + majority)
## to 163 is 747 (comprehensive)
## to 180 is 3673.5 (ensemble)

## Just comprehensive: 1861.5
## Just majority:
## Just one iter local: 

run_sim <- function(
        n_subjects, 
        n_samples_per_subject, 
        n_swap_cats,
        fraction_mislabel, 
        fraction_anchor, 
        fraction_ghost,
        seed,
        output_dir,
        sim_name) {
    
    ## Create tracker to show we tried this sim
    file.create(file.path(output_dir, "01_tried", sim_name))
    
    print(glue("Running simulation {sim_name}"))
    
    n_subjects_per_group <- as.integer(n_subjects/2)
    n_samples_per_group <- n_subjects_per_group * n_samples_per_subject
    n_samples <- n_samples_per_group * 2
    n_mislabels <- as.integer(fraction_mislabel * n_samples)
    n_anchor_samples <- as.integer(fraction_anchor * (n_samples - n_mislabels))
    n_ghost_samples <- as.integer(fraction_ghost * n_samples)
    
    sample_meta_data_params <- list(
        n_subjects_per_group = n_subjects_per_group,
        n_samples_per_group = n_samples_per_group,
        n_swap_cats = n_swap_cats,
        n_mislabels = n_mislabels,
        seed = seed
    )
    results_list <- do.call(sim_mislabeled_sample_meta_data, sample_meta_data_params)
    sample_meta_data <- results_list$sample_meta_data
    swap_cats <- results_list$swap_cats
    
    sample_genotype_data <- sample_meta_data %>%
        select(Sample_ID, Subject_ID, Genotype_Group_ID)
    ## Delete genotype information for ghost samples
    ghost_samples <- sample(sample_meta_data %>% pull(Sample_ID), n_ghost_samples, replace=FALSE)
    sample_genotype_data[sample_genotype_data$Sample_ID %in% ghost_samples, "Genotype_Group_ID"] <- NA_character_
    anchor_samples <- sample(sample_meta_data %>% filter(!Mislabeled) %>% pull(Sample_ID), n_anchor_samples, replace=FALSE)
    
    mislabel_solver <- new("MislabelSolver", sample_genotype_data, swap_cats, anchor_samples)
    results_df <- sample_meta_data %>% 
        select(Sample_ID, Subject_ID, Genotype_Group_ID, E_Sample_ID, E_Subject_ID, Mislabeled) %>%
        left_join(swap_cats, by="Sample_ID") %>% 
        rename(Init_Sample_ID = Sample_ID,
               Init_Subject_ID = Subject_ID,
               True_Sample_ID = E_Sample_ID,
               True_Subject_ID = E_Subject_ID) %>% 
        mutate(
            Genotype_Group_ID = ifelse(Init_Sample_ID %in% ghost_samples, NA_character_, Genotype_Group_ID),
            Is_Anchor = Init_Sample_ID %in% anchor_samples
        )
    
    n_total_subject_mislabels <- sum((results_df$Init_Subject_ID != results_df$True_Subject_ID) & !is.na(results_df$Genotype_Group_ID))
    n_total_sample_mislabels <- sum((results_df$Init_Sample_ID != results_df$True_Sample_ID) & !is.na(results_df$Genotype_Group_ID))
    print(glue("{n_total_subject_mislabels} total subject mislabels"))
    print(glue("{n_total_sample_mislabels} total sample mislabels"))
    
    start_time = Sys.time()
    
    ## Create solve results
    # 1. Baseline majority search
    mislabel_solver <- solve_majority_search(mislabel_solver, unambiguous_only=TRUE)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Init_Component_ID, Sample_ID, Subject_ID, Solved, Component_ID) %>% 
        rename(Baseline_Component_ID = Component_ID,
               Sample_ID_baseline = Sample_ID,
               Subject_ID_baseline = Subject_ID,
               Solved_baseline = Solved)
    results_df <- results_df %>% full_join(curr_results_df, by="Init_Sample_ID")
    print(glue("Baseline run for {sim_name}"))
    
    n_baseline_subject_mislabels <- sum((results_df$Subject_ID_baseline != results_df$True_Subject_ID) & !is.na(results_df$Genotype_Group_ID))
    n_baseline_sample_mislabels <- sum((results_df$Sample_ID_baseline != results_df$True_Sample_ID) & !is.na(results_df$Genotype_Group_ID))
    print(glue("{n_baseline_subject_mislabels} subject mislabels after baseline"))
    print(glue("{n_baseline_sample_mislabels} sample mislabels after baseline"))
    
    # 2. Majority search with cycles
    mislabel_solver <- solve_majority_search(mislabel_solver)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Sample_ID, Subject_ID, Solved, Component_ID) %>% 
        rename(Majority_Component_ID = Component_ID,
               Sample_ID_majority = Sample_ID,
               Subject_ID_majority = Subject_ID,
               Solved_majority = Solved)
    results_df <- results_df %>% full_join(curr_results_df, by="Init_Sample_ID")
    print(glue("Majority run for {sim_name}"))
    
    n_majority_subject_mislabels <- sum((results_df$Subject_ID_majority != results_df$True_Subject_ID) & !is.na(results_df$Genotype_Group_ID))
    n_majority_sample_mislabels <- sum((results_df$Sample_ID_majority != results_df$True_Sample_ID) & !is.na(results_df$Genotype_Group_ID))
    print(glue("{n_majority_subject_mislabels} subject mislabels after majority"))
    print(glue("{n_majority_sample_mislabels} sample mislabels after majority"))
    
    # 3. Majority search with comprehensive
    mislabel_solver <- solve_comprehensive_search(mislabel_solver)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Sample_ID, Subject_ID, Solved, Component_ID) %>% 
        rename(Majority_Comprehensive_Component_ID = Component_ID,
            Sample_ID_majority_comprehensive = Sample_ID,
               Subject_ID_majority_comprehensive = Subject_ID,
               Solved_majority_comprehensive = Solved)
    results_df <- results_df %>% full_join(curr_results_df, by="Init_Sample_ID")
    print(glue("Majority w comprehensive run for {sim_name}"))
    
    n_comprehensive_subject_mislabels <- sum((results_df$Subject_ID_majority_comprehensive != results_df$True_Subject_ID) & !is.na(results_df$Genotype_Group_ID))
    n_comprehensive_sample_mislabels <- sum((results_df$Sample_ID_majority_comprehensive != results_df$True_Sample_ID) & !is.na(results_df$Genotype_Group_ID))
    print(glue("{n_comprehensive_subject_mislabels} subject mislabels after comprehensive"))
    print(glue("{n_comprehensive_sample_mislabels} sample mislabels after comprehensive"))
    
    if (max(table(results_df$Majority_Comprehensive_Component_ID)) > 5000) {
        file.create(file.path(output_dir, "02_comp_over_5000", sim_name))
    }
    
    # 4. Majority search iterative ensemble
    mislabel_solver <- solve(mislabel_solver)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Sample_ID, Subject_ID, Solved, Component_ID) %>% 
        rename(Ensemble_Component_ID = Component_ID,
               Sample_ID_ensemble = Sample_ID,
               Subject_ID_ensemble = Subject_ID,
               Solved_ensemble = Solved)
    results_df <- results_df %>% full_join(curr_results_df, by="Init_Sample_ID")
    print(glue("Ensemble run for {sim_name}"))
    
    n_ensemble_subject_mislabels <- sum((results_df$Subject_ID_ensemble != results_df$True_Subject_ID) & !is.na(results_df$Genotype_Group_ID))
    n_ensemble_sample_mislabels <- sum((results_df$Sample_ID_ensemble != results_df$True_Sample_ID) & !is.na(results_df$Genotype_Group_ID))
    print(glue("{n_ensemble_subject_mislabels} subject mislabels after ensemble"))
    print(glue("{n_ensemble_sample_mislabels} sample mislabels after ensemble"))
    
    end_time = Sys.time()
    
    run_time <- end_time - start_time
    writeLines(as.character(run_time), file.path(output_dir, "04_runtimes", sim_name))
    
    output_path <- file.path(output_dir, "05_solve_results", paste0(sim_name, ".csv"))
    write.csv(results_df, output_path)

    print(glue("Job complete for {sim_name}"))
    return(results_df)
}

## For all sims
## If over 5000 after comprehensive search, write sim name to directory "over_5000" and quit the sim
## Try local search
## If failed because of timeout, write sim name to directory "timeout_2hours" and quit the sim
## If succeeded, 
## Failures due to memory are setdiff over all sims with every outcome above

params_grid <- readRDS(params_grid_file)
for (i in 1:nrow(params_grid)) {
    args_list <- as.list(params_grid[i, ])
    sim_name <- args_list$sim_name
    if (file.exists(file.path(args_list$output_dir, "01_tried", sim_name))) {next}
    tryCatch(
        expr = {
            withTimeout(do.call(run_sim, args_list), timeout=2 * 60 * 60)
        },
        error = function(e){ 
            errorMessage <- paste("Error for sim_name:", sim_name, "\n", "Error message:", conditionMessage(e))
            if (grepl("reached elapsed time limit", errorMessage)) {
                ## Create tracker to show this sim timed out
                file.create(file.path(args_list$output_dir, "03_timeout_2hrs", sim_name))
            }
            print(errorMessage)
        }
    )
}


