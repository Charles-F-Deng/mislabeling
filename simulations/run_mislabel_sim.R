library(rlang)
library(withr)
library(rex)
library(dplyr)
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
as_data_frame <- igraph::as_data_frame
library(assertthat)
library(limma)
column_to_rownames <- tibble::column_to_rownames
rownames_to_column <- tibble::rownames_to_column

setwd("/Users/charlesdeng/Workspace/mislabeling/simulations")
source("sim_mislabeled_dataset.R")
source("MislabelSolver.R")

## Example values for testing
# Number of subjects (10 - 10000, by=100 for first 1000, by 1000 after) (20)
# Number of samples_per_subject (2 to 10, by=1) (10)
# Number of swap_cat (1 to 5 groups, 5 to 20) (10)
# Pct mislabeled (1 mislabel (first 10% by 1%, next by 5%, then by 10) to 100%) (20)
# Pct anchor samples (0 to 10%, by=2) (5)
# Pct ghost samples (0 to 10%, by=2) (5)
# Pct deletions (can be left out of the grid)
# Contaminated samples (like a ghost sample)
# 5 seeds (5)

# Column of relabels for each incremental improvement
# 1. unambiguous majority
# 2. majority with cycle detection
# 3. majority with comprehensive
# 4. iterative ensemble with local search

run_sim <- function(
        n_subjects, 
        n_samples_per_subject, 
        n_swap_cats,
        fraction_mislabel, 
        fraction_anchor, 
        fraction_ghost,
        seed,
        output_dir) {
    
    n_subjects_per_group <- as.integer(n_subjects/2)
    n_samples_per_group <- n_subjects_per_group * n_samples_per_subject
    n_samples <- n_samples_per_group * 2
    n_mislabels <- as.integer(fraction_mislabel * n_samples)
    n_anchor_samples <- as.integer(fraction_anchor * (n_samples - n_mislabels))
    n_ghost_samples <- as.integer(fraction_ghost * n_samples)
    
    elist_params <- list(
        n_subjects_per_group = n_subjects_per_group,
        n_samples_per_group = n_samples_per_group,
        n_swap_cats = n_swap_cats,
        n_mislabels = n_mislabels,
        seed = seed,
        ## The rest of these inputs don't matter for the purpose of evaluating the MislabelSolver
        n_features = 10,    
        n_sv = 1,
        fraction_de_case = 0.5,
        case_sd = 0.1,
        subject_sd = 0.1,
        sex_sd = 0.1,
        age_sd = 0.1,
        sv_sd = 0.1,
        resid_sd = 0.05
    )
    elist <- do.call(sim_mislabeled_data, elist_params)
    
    sample_genotype_data <- elist$targets %>%
        select(Sample_ID, Subject_ID, Genotype_Group_ID)
    ## Delete genotype information for ghost samples
    ghost_samples <- sample(rownames(sample_genotype_data), n_ghost_samples, replace=FALSE)
    sample_genotype_data[ghost_samples, "Genotype_Group_ID"] <- NA_character_
    rownames(sample_genotype_data) <- NULL
    swap_cats <- elist$swap_cats
    anchor_samples <- sample(elist$targets %>% filter(!Mislabeled) %>% pull(Sample_ID), n_anchor_samples)
    
    mislabel_solver_initial <- new("MislabelSolver", sample_genotype_data, swap_cats, anchor_samples)
    results_df <- elist$targets %>% 
        select(Sample_ID, Subject_ID, Genotype_Group_ID, E_Sample_ID, E_Subject_ID, Mislabeled) %>% 
        rename(Init_Sample_ID = Sample_ID,
               Init_Subject_ID = Subject_ID,
               True_Sample_ID = E_Sample_ID,
               True_Subject_ID = E_Subject_ID)
    rownames(results_df) <- NULL
    join_cols <- c("Init_Sample_ID", "Init_Subject_ID", "Genotype_Group_ID", "Init_Component_ID")
    
    ## Create solve results
    # 1. Baseline majority search
    mislabel_solver <- mislabel_solver_initial
    mislabel_solver <- solve_majority_search(mislabel_solver, unambiguous_only = TRUE)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Init_Component_ID, Sample_ID, Subject_ID, Solved) %>% 
        rename(Sample_ID_baseline = Sample_ID,
               Subject_ID_baseline = Subject_ID,
               Solved_baseline = Solved)
    results_df <- results_df %>% full_join(curr_results_df)
    
    # 2. Majority search with cycles
    mislabel_solver <- mislabel_solver_initial
    mislabel_solver <- solve_majority_search(mislabel_solver)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Sample_ID, Subject_ID, Solved) %>% 
        rename(Sample_ID_majority = Sample_ID,
               Subject_ID_majority = Subject_ID,
               Solved_majority = Solved)
    results_df <- results_df %>% full_join(curr_results_df)
    
    # 3. Majority search with comprehensive
    mislabel_solver <- mislabel_solver_initial
    mislabel_solver <- solve_comprehensive_search(mislabel_solver)
    mislabel_solver <- solve_majority_search(mislabel_solver)
    mislabel_solver <- solve_comprehensive_search(mislabel_solver)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Sample_ID, Subject_ID, Solved) %>% 
        rename(Sample_ID_majority_comprehensive = Sample_ID,
               Subject_ID_majority_comprehensive = Subject_ID,
               Solved_majority_comprehensive = Solved)
    results_df <- results_df %>% full_join(curr_results_df)
    
    # 4. Majority search iterative ensemble
    mislabel_solver <- mislabel_solver_initial
    mislabel_solver <- solve(mislabel_solver)
    curr_results_df <- mislabel_solver@.solve_state$relabel_data %>% 
        select(Init_Sample_ID, Sample_ID, Subject_ID, Solved) %>% 
        rename(Sample_ID_ensemble = Sample_ID,
               Subject_ID_ensemble = Subject_ID,
               Solved_ensemble = Solved)
    results_df <- results_df %>% full_join(curr_results_df)
    
    run_sim_params <- c(n_subjects, n_samples_per_subject, n_swap_cats, fraction_mislabel,
                fraction_anchor, fraction_ghost, seed)
    output_file <- paste0(paste(run_sim_params, collapse = "-"), ".csv")
    output_path <- file.path(output_dir, output_file)
    write.csv(results_df, output_path)
}

run_sim(10, 4, 3, 0.10, 0.05, 0.05, 1986, output_dir="/Users/charlesdeng/Workspace/mislabeling/simulations")

