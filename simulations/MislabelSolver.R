# library(visNetwork)
# library(gtools)
# library(tidyr)
# library(igraph)
# library(withr)
# library(dplyr)
# library(assertthat)
# library(reshape2)
# library(glue)

## TODO
# swap_cat_shape - only need it for plotting. appropriate to put it in relabel_data or does it belong somewhere else
# Add a panel of tests to run automatically and compare against baselines
# relabels: fix validation to ensure you can only relabel unsolved samples
# local_search: change find_neighbors to take in parameters unsolved_relabel_data. Then switch local search to be component-by-component
# plotting: for the unsolved option, include samples where label is not found
# plotting: redo plot_corrections. Also, maybe add swaps created on top of original graph
# For label not found option, I think right now its possible to make duplicates
# create a summary function
# Majority search doesn't currently allow for ghosts or deletions. Fix this in find_relabel_samples_from_putative_subjects
# Revisit cycle detection algorithm
# solve_comprehensive_search: what happens when number of genotype groups exceeds number of subjects?
# misc: Don't run set.seed(1), switch everything to with_seed
# misc: Standardize names for unsolved_relabel_data, relabel_data, unsolved_ghost_data, etc.

## DONE
# structure: have user pass in SwapCat_ID as part sample_genotype_data, refactor to be called sample_metadata. Then remove swap_cat passes
# initialize: assign component_ID even if the component is solved to begin with
# plotting: use red vertices for samples where label is not found
# comprehensive_search: fix bug of n > r
# find_neighbors: Memory optimization by doing sparse matrices for distance matrix
# find_neighbors: Memory and runtime optimization by computing distance with matrix multiplication
# solve_local_search: time optimization by doing a bulk search for local search

EMPTY_RELABELS <- data.frame(relabel_from=character(0), relabel_to=character(0))
VISNETWORK_SWAPCAT_SHAPES <- c("dot", "square", "triangle", "diamond", "star")
LABEL_NOT_FOUND <- "LABELNOTFOUND"
MAX_GENOTYPES_COMP_SEARCH <- 7
MAX_NEIGHBORS_LOCAL_SEARCH <- 1e6

setClass("MislabelSolver",
         representation(
             sample_genotype_data = "data.frame",
             swap_cats = "data.frame",
             anchor_samples = "character",
             .solve_state = "list"
         ),
         prototype(
             swap_cats = NULL,
             anchor_samples = character(0)
         )
)

setMethod("initialize", "MislabelSolver",
          function(.Object, sample_genotype_data, swap_cats=NULL, anchor_samples=character(0)) {
              ## Convert and validate inputs
              sample_genotype_data <- as.data.frame(lapply(sample_genotype_data, as.character))
              .validate_sample_genotype_data(sample_genotype_data)
              
              if (is.null(swap_cats)) {
                  swap_cats <- sample_genotype_data[, "Sample_ID", drop=FALSE]
                  swap_cats$SwapCat_ID <- "SwapCat1"
              }
              swap_cats <- as.data.frame(lapply(swap_cats, as.character))
              .validate_swap_cats(sample_genotype_data, swap_cats)
              ## Provided there are enough shapes, assign a unique shape to each SwapCat_ID
              all_swap_cat_ids <- names(sort(table(swap_cats$SwapCat_ID), decreasing=TRUE))
              swap_cat_shapes <- data.frame(
                  SwapCat_ID = all_swap_cat_ids,
                  SwapCat_Shape = "dot"
              )
              if (length(all_swap_cat_ids) <= length(VISNETWORK_SWAPCAT_SHAPES)) {
                  swap_cat_shapes$SwapCat_Shape = VISNETWORK_SWAPCAT_SHAPES[1:length(all_swap_cat_ids)]
              }
              swap_cats <- swap_cats %>% 
                  left_join(swap_cat_shapes, by="SwapCat_ID")
              
              ## Wrangle 'anchor_samples'
              anchor_samples <- unique(as.character(anchor_samples))
              .validate_anchor_samples(sample_genotype_data, anchor_samples)
              
              ## Initialize object 'solve_state'
              relabel_data <- sample_genotype_data %>% 
                  mutate(
                      Init_Sample_ID = Sample_ID,
                      Init_Subject_ID = Subject_ID,
                      Solved = FALSE
                  ) %>% 
                  left_join(swap_cats, by="Sample_ID")
              unsolved_relabel_data <- relabel_data %>% filter(!is.na(Genotype_Group_ID))
              unsolved_ghost_data <- relabel_data %>% filter(is.na(Genotype_Group_ID))
              putative_subjects <- data.frame(Genotype_Group_ID = character(0),
                                              Subject_ID = character(0))
              ambiguous_subjects <- list()
              solve_state <- list(
                  relabel_data = relabel_data,
                  unsolved_relabel_data = unsolved_relabel_data,
                  unsolved_ghost_data = unsolved_ghost_data,
                  putative_subjects = putative_subjects,
                  ambiguous_subjects = ambiguous_subjects
              )
              
              .Object@sample_genotype_data <- sample_genotype_data
              .Object@swap_cats <- swap_cats
              .Object@anchor_samples <- anchor_samples
              .Object@.solve_state <- solve_state
              
              .Object <- .update_solve_state(.Object, initialization=TRUE)
              return(.Object)
          }
)

setGeneric("solve_majority_search", function(object, unambiguous_only=FALSE) {
    standardGeneric("solve_majority_search")
})

setGeneric("solve_comprehensive_search", function(object) {
    standardGeneric("solve_comprehensive_search")
})

setGeneric("solve_local_search", function(object, ...) {
    standardGeneric("solve_local_search")
})

setGeneric("solve_local_search_bulk", function(object, ...) {
    standardGeneric("solve_local_search_bulk")
})

setGeneric("solve", function(object, ...) {
    standardGeneric("solve")
})

setGeneric("write_corrections", function(object) {
    standardGeneric("write_corrections")
})

setMethod("solve_majority_search", "MislabelSolver",
          function(object, unambiguous_only=FALSE) {
              set.seed(1)
              print("Starting majority search")
              if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                  print("0 samples relabeled")
                  return(object)
              }
              
              ## 1. Update putative subjects
              votes <- table(object@.solve_state$unsolved_relabel_data$Genotype_Group_ID, 
                             object@.solve_state$unsolved_relabel_data$Subject_ID)
              votes_by_genotype <- data.frame(
                  Genotype_Group_ID = rownames(votes),
                  Max_Subject_ID = colnames(votes)[apply(votes, 1, which.max)],
                  n = rowSums(votes),
                  n_Max_Subject_ID = apply(votes, 1, max)
              ) %>%
                  filter(
                      n_Max_Subject_ID >= 2,
                      n_Max_Subject_ID > n/2
                  ) %>% 
                  rename(Subject_ID = Max_Subject_ID) %>% 
                  select(Genotype_Group_ID, Subject_ID)
              votes_by_subject <- data.frame(
                  Subject_ID = colnames(votes),
                  Max_Genotype_Group_ID = rownames(votes)[apply(votes, 2, which.max)],
                  n = colSums(votes),
                  n_Max_Genotype_Group_ID = apply(votes, 2, max)
              ) %>% 
                  filter(
                      n_Max_Genotype_Group_ID >= 2,
                      n_Max_Genotype_Group_ID > n/2
                  ) %>% 
                  rename(Genotype_Group_ID = Max_Genotype_Group_ID) %>% 
                  select(Subject_ID, Genotype_Group_ID)
              new_putative_subjects <- inner_join(votes_by_subject, votes_by_genotype, by=c("Genotype_Group_ID", "Subject_ID")) %>% 
                  anti_join(object@.solve_state$putative_subjects, by=c("Genotype_Group_ID", "Subject_ID"))
              object <- .update_putative_subjects(object, new_putative_subjects)
              
              ## 2. Find relabel cycles
              unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
              putative_subjects <- object@.solve_state$putative_subjects
              relabels <- .find_relabel_cycles_from_putative_subjects(unsolved_relabel_data, 
                                                                      putative_subjects,
                                                                      unambiguous_only=unambiguous_only,
                                                                      allow_unknowns=FALSE)
              
              ## 3. Relabel samples and update solve state
              object <- .relabel_samples(object, relabels)
              print(paste0(nrow(relabels), " samples relabeled"))
              return(object)
          }
)

setMethod("solve_comprehensive_search", "MislabelSolver",
          function(object) {
              set.seed(1)
              print("Starting comprehensive search")
              if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                  print("0 samples relabeled")
                  return(object)
              }
              
              ## 1. Update putative subjects
              component_ids <- sort(unique(object@.solve_state$unsolved_relabel_data$Component_ID))
              for (component_id in component_ids) {
                  cc_unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data %>% filter(Component_ID == component_id)
                  cc_unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data %>% filter(Component_ID == component_id)
                  cc_sample_ids <- c(cc_unsolved_relabel_data$Sample_ID, cc_unsolved_ghost_data$Sample_ID)
                  cc_swap_cat_ids <- unique(c(cc_unsolved_relabel_data$SwapCat_ID, cc_unsolved_ghost_data$SwapCat_ID))
                  cc_genotypes <- unique(cc_unsolved_relabel_data$Genotype_Group_ID)
                  cc_subjects <- unique(cc_unsolved_relabel_data$Subject_ID)
                  
                  ## For now, pass out of components where number of Genotype_Group(s) is greater than number of Subject_ID(s)
                  if (length(cc_genotypes) > length(cc_subjects)) {
                      # print(component_id)
                      next
                  }
                  
                  ## Lock genotypes that already have a putative subject assigned, and find all possible permutations for free genotypes
                  locked_genotypes <- intersect(putative_subjects$Genotype_Group_ID, cc_genotypes)
                  locked_subjects <- intersect(putative_subjects$Subject_ID, cc_subjects)
                  free_genotypes <- setdiff(cc_genotypes, locked_genotypes)
                  free_subjects <- setdiff(cc_subjects, locked_subjects)
                  
                  if (length(free_genotypes) > MAX_GENOTYPES_COMP_SEARCH) {
                      next
                  }
                  
                  if (length(free_genotypes) > 0 & length(free_subjects) > 0) {
                      n <- length(free_subjects)
                      r <- length(free_genotypes)
                      perm_genotypes <- permutations(n, r, free_subjects)
                      colnames(perm_genotypes) <- sort(free_genotypes)
                      n_perm <- nrow(perm_genotypes)
                      for (locked_genotype_id in locked_genotypes) {
                          locked_subject_id <- putative_subjects[putative_subjects$Genotype_Group_ID == locked_genotype_id, "Subject_ID"][[1]]
                          new_perm_col <- matrix(data=locked_subject_id, ncol=1, nrow=n_perm, dimnames=list(NULL, locked_genotype_id))
                          perm_genotypes <- cbind(perm_genotypes, new_perm_col)
                      }
                  } else {
                      locked_putative_subjects <- putative_subjects[putative_subjects$Genotype_Group_ID %in% locked_genotypes, ]
                      perm_genotypes <- t(locked_putative_subjects$Subject_ID)
                      colnames(perm_genotypes) <- locked_putative_subjects$Genotype_Group_ID
                  }
                  perm_genotypes <- as.matrix(perm_genotypes, dimnames=c("Permutation_ID", "Genotype_Group_ID"))
                  n_perms <- nrow(perm_genotypes)
                  permutation_ids <- paste0("Permutation", formatC(1:n_perms, width=nchar(n_perms), format="d", flag="0"))
                  rownames(perm_genotypes) <- permutation_ids
                  
                  ## For each Genotype_Group_ID/Subject_ID permutation, determine
                  ## 1. The number of existing samples to relabel
                  ## 2. The number of ghost samples needed to add
                  ## 3. The number of indels required after ghost samples are included
                  label_counts <- cc_unsolved_relabel_data %>% 
                      select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
                      group_by(Subject_ID, SwapCat_ID) %>% 
                      summarize(n_labels = n(), .groups="drop") 
                  ghost_label_counts <- cc_unsolved_ghost_data %>% 
                      select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
                      group_by(Subject_ID, SwapCat_ID) %>% 
                      summarize(n_ghost_labels = n(), .groups="drop") 
                  genotype_counts <- cc_unsolved_relabel_data %>% 
                      select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
                      group_by(Genotype_Group_ID, SwapCat_ID) %>% 
                      summarize(n_in_genotype = n(), .groups="drop")
                  genotype_subject_concordant_counts <- cc_unsolved_relabel_data %>% 
                      select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
                      group_by(Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
                      summarize(n_samples_correct = n(), .groups="drop")
                  
                  ## Create a "long" version of perm_genotypes
                  long_perm_genotypes <- melt(perm_genotypes)
                  colnames(long_perm_genotypes) <- c("Permutation_ID", "Genotype_Group_ID", "Subject_ID")
                  long_perm_genotypes <- as.data.frame(lapply(long_perm_genotypes, as.character))
                  
                  ## Create an empty matrix of stats for each permutation
                  count_cols <- c("n_labels", "n_ghost_labels", "n_in_genotype", "n_samples_correct")
                  permutation_stats <- matrix(0, nrow=n_perms, ncol=length(count_cols))
                  colnames(permutation_stats) <- count_cols
                  rownames(permutation_stats) <- permutation_ids
                  for (swap_cat_id in cc_swap_cat_ids) {
                      long_perm_genotypes$SwapCat_ID <- swap_cat_id
                      
                      ## Join in ascending order of size
                      merged_long_perm_genotypes <- long_perm_genotypes %>% 
                          left_join(label_counts, by=c("Subject_ID", "SwapCat_ID")) %>% 
                          left_join(ghost_label_counts, by=c("Subject_ID", "SwapCat_ID")) %>% 
                          left_join(genotype_counts, by=c("Genotype_Group_ID", "SwapCat_ID")) %>% 
                          left_join(genotype_subject_concordant_counts, by=c("Subject_ID", "Genotype_Group_ID", "SwapCat_ID")) %>% 
                          mutate_at(
                              vars(all_of(count_cols)),
                              ~coalesce(., 0)
                          )
                      merged_long_perm_genotypes <- merged_long_perm_genotypes %>% 
                          mutate(
                              ## In each genotype group, the number of samples that can be relabeled to a genotyped sample
                              n_samples_to_relabel = pmin(n_in_genotype, n_labels) - n_samples_correct,
                              ## In each genotype group, the number of mislabeled samples outstanding is
                              ## n_genotypes - n_correct - n_samples_to_relabel. We try to plug this gap with ghost samples
                              n_samples_to_relabel_ghost = pmin(n_in_genotype - n_samples_correct - n_samples_to_relabel, n_ghost_labels),
                              n_label_deletions = pmax(0, n_in_genotype - n_labels - n_ghost_labels),
                              n_genotype_deletions = pmax(0, n_labels - n_in_genotype)
                          ) %>% 
                          arrange(Permutation_ID)
                      
                      ## Evaluate each permutation
                      swap_cat_perm_stats <- merged_long_perm_genotypes %>% 
                          group_by(Permutation_ID) %>% 
                          summarize(
                              n_samples_correct = sum(n_samples_correct),
                              n_samples_to_relabel = sum(n_samples_to_relabel),
                              n_samples_to_relabel_ghost = sum(n_samples_to_relabel_ghost),
                              n_genotype_deletions = sum(n_genotype_deletions),
                              n_label_deletions = sum(n_label_deletions),
                              .groups = "drop"               
                          ) %>%
                          mutate(
                              n_samples_to_relabel = n_samples_to_relabel + pmin(n_genotype_deletions, n_samples_to_relabel_ghost),
                              n_genotype_deletions = pmax(0, n_genotype_deletions - n_samples_to_relabel_ghost),
                              ## The weighting scheme is arbitrary right now
                              perm_score = n_samples_to_relabel + 1.5 * n_samples_to_relabel_ghost + 2 * (n_genotype_deletions + n_label_deletions)
                          ) %>% 
                          tibble::column_to_rownames("Permutation_ID")
                      
                      permutation_stats <- permutation_stats + swap_cat_perm_stats
                  }
                  
                  permutation_stats <- permutation_stats %>%
                      as.data.frame() %>% 
                      tibble::rownames_to_column("Permutation_ID") %>% 
                      arrange(perm_score)
                  
                  ## To find a single solution, take top row
                  ## TODO: record any ties
                  best_permutation <- perm_genotypes[permutation_stats$Permutation_ID[1], , drop=FALSE]
                  
                  new_putative_subjects <- best_permutation %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      tibble::rownames_to_column()
                  colnames(new_putative_subjects) <- c("Genotype_Group_ID", "Subject_ID")
                  ## Specify when a Subject_ID in the component doesn't have a Genotype_Group_ID, or vice versa
                  if (length(cc_genotypes) > length(cc_subjects)) {
                      unmatched_genotypes <- setdiff(cc_genotypes, new_putative_subjects$Genotype_Group_ID)
                      new_putative_subjects <- rbind(new_putative_subjects, 
                                                     data.frame(Genotype_Group_ID = unmatched_genotypes, Subject_ID = NA_character_))
                  }
                  if (length(cc_genotypes) < length(cc_subjects)) {
                      unmatched_subjects <- setdiff(cc_subjects, new_putative_subjects$Subject_ID)
                      new_putative_subjects <- rbind(new_putative_subjects, 
                                                     data.frame(Genotype_Group_ID = NA_character_, Subject_ID = unmatched_subjects))
                  }
                  new_putative_subjects <- new_putative_subjects %>%  
                      anti_join(object@.solve_state$putative_subjects, by=c("Genotype_Group_ID", "Subject_ID"))
                  object <- .update_putative_subjects(object, new_putative_subjects)
              }
              
              ## Find relabel cycles
              unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
              unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
              putative_subjects <- object@.solve_state$putative_subjects
              relabels <- .find_relabel_cycles_from_putative_subjects(unsolved_relabel_data, putative_subjects, 
                                                                      unsolved_ghost_data, allow_unknowns=TRUE)
              
              ## Relabel samples and update solve state
              object <- .relabel_samples(object, relabels) 
              print(paste0(nrow(relabels), " samples relabeled"))
              return(object)
          }
)

# 1. Create a 2D table of Genotype_Group_ID x Subject_ID counts (2D table)
# 2. Compute baseline scaled genotypes for each Genotype_Group (named vector)
# 3. Find all neighbors
# 4. Write a function that takes in 2 samples to swap, and finds the delta entropy and returns
# 5. mapply will give the delta column
setMethod("solve_local_search", "MislabelSolver", 
          function(object, n_iter=1, include_ghost=FALSE, filter_concordant_vertices=FALSE) {
              set.seed(1)
              print("Starting local search")
              calc_scaled_entropy <- function(x) {
                  return(sum(x*log(x/sum(x)), na.rm=TRUE))
              }
              
              for (i in 1:n_iter) {
                  print(glue("Iteration {i}:: local search, 'include_ghost'={include_ghost}, 'filter_concordant_vertices'={filter_concordant_vertices}"))
                  unsolved_all_data <- rbind(object@.solve_state$unsolved_relabel_data,
                                             object@.solve_state$unsolved_ghost_data)
                  if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                      print("0 samples relabeled")
                      return(object)
                  }
                  votes <- table(object@.solve_state$unsolved_relabel_data$Genotype_Group_ID, 
                                 object@.solve_state$unsolved_relabel_data$Subject_ID)
                  base_entropies <- apply(votes, MARGIN=1, calc_scaled_entropy)
                  
                  calc_swapped_delta_entropy <- function(swap_from_subject, swap_from_genotype,
                                                         swap_to_subject, swap_to_genotype) {
                      delta <- 0
                      if (!is.na(swap_from_genotype)) {
                          genotype_base_entropy <- base_entropies[swap_from_genotype]
                          genotype_votes_vec <- votes[swap_from_genotype, ]
                          genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[swap_from_subject] - 1
                          genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[swap_to_subject] + 1
                          genotype_new_entropy <- calc_scaled_entropy(genotype_votes_vec)
                          delta <- delta + genotype_new_entropy - genotype_base_entropy
                      }
                      if (!is.na(swap_to_genotype)) {
                          genotype_base_entropy <- base_entropies[swap_to_genotype]
                          genotype_votes_vec <- votes[swap_to_genotype, ]
                          genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[swap_to_subject] - 1
                          genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[swap_from_subject] + 1
                          genotype_new_entropy <- calc_scaled_entropy(genotype_votes_vec)
                          delta <- delta + genotype_new_entropy - genotype_base_entropy
                      }
                      return(delta)
                  }
                  
                  neighbors <- .find_neighbors(object, include_ghost, filter_concordant_vertices) %>% 
                      left_join(
                          unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")], 
                          by=c("Sample_A"="Sample_ID")
                      ) %>%
                      dplyr::rename(Subject_A = Subject_ID, Genotype_Group_A = Genotype_Group_ID) %>% 
                      left_join(
                          unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID", "Component_ID")], 
                          by=c("Sample_B"="Sample_ID")
                      ) %>% 
                      dplyr::rename(Subject_B = Subject_ID, Genotype_Group_B = Genotype_Group_ID) 
                  
                  print(glue("{nrow(neighbors)} candidate swaps being evaluated..."))
                  all_component_ids <- sort(unique(object@.solve_state$unsolved_relabel_data$Component_ID))
                  relabels <- data.frame(matrix(data=NA, nrow=length(all_component_ids), ncol=2, dimnames=list(c(), c("relabel_from", "relabel_to"))))
                  curr_idx <- 1
                  for (curr_component_id in all_component_ids) {
                      cc_relabel_data <- unsolved_all_data %>% filter(Component_ID == curr_component_id)
                      cc_neighbors <- neighbors %>% filter(Component_ID == curr_component_id)
                      
                      if (nrow(cc_neighbors) == 0) {next}
                      cc_neighbor_objectives <- cc_neighbors %>% 
                          mutate(
                              delta = mapply(calc_swapped_delta_entropy, 
                                             swap_from_subject=Subject_A, 
                                             swap_from_genotype=Genotype_Group_A,
                                             swap_to_subject=Subject_B,
                                             swap_to_genotype=Genotype_Group_B)
                          )
                      cc_relabels <- cc_neighbor_objectives %>%
                          filter(delta > 0, delta == max(delta)) 
                      if (nrow(cc_relabels) == 0) {next}
                      cc_relabels <- cc_relabels %>% 
                          sample_n(1) %>%
                          transmute(
                              relabel_from=Sample_A,
                              relabel_to=Sample_B
                          )
                      relabels[curr_idx, c("relabel_from", "relabel_to")] <- cc_relabels
                      curr_idx <- curr_idx + 1
                  }
                  
                  relabels <- relabels[!is.na(relabels[, 1]) & !is.na(relabels[, 2]), ]
                  relabels <- rbind(relabels, data.frame(relabel_from=relabels$relabel_to, relabel_to=relabels$relabel_from))
                  object <- .relabel_samples(object, relabels)
                  print(paste0(nrow(relabels), " samples relabeled"))
              }
              
              return(object)
          }
)

# 1. Create a 2D table of Genotype_Group_ID x Subject_ID counts (2D table)
# 2. Compute baseline scaled genotypes for each Genotype_Group (named vector)
# 3. Find all neighbors
# 4. Iterate through components
# 5. If component has more than MAX_NEIGHBORS_LOCAL_SEARCH neighbors, sample down to MAX_NEIGHBORS_LOCAL_SEARCH
# 6. mapply to give the delta column for each neighbor
# 7. Find minimum number of mislabels in the component
# 8. Find number of swaps to look for (n swaps fixes at most 2n mislabels, so that must be less that frac * min_mislabels_in_component)
# 9. Iteratively find a swap, then drop all swaps that involve either of the 2 genotypes that you used
# 10. rbind all swaps found for the component
setMethod("solve_local_search_bulk", "MislabelSolver", 
          function(object, n_iter=1, frac_per_iter=0.10, include_ghost=FALSE, filter_concordant_vertices=FALSE) {
              set.seed(1)
              print("Starting bulk local search")
              calc_scaled_entropy <- function(x) {
                  return(sum(x*log(x/sum(x)), na.rm=TRUE))
              }
              
              for (i in 1:n_iter) {
                  if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                      return(object)
                  }
                  
                  unsolved_all_data <- rbind(object@.solve_state$unsolved_relabel_data,
                                             object@.solve_state$unsolved_ghost_data)
                  
                  print(glue("Iteration {i}:: bulk local search, 'frac_per_iter'={frac_per_iter}, 'include_ghost'={include_ghost}, 'filter_concordant_vertices'={filter_concordant_vertices}"))
                  
                  votes <- table(object@.solve_state$unsolved_relabel_data$Genotype_Group_ID, 
                                 object@.solve_state$unsolved_relabel_data$Subject_ID)
                  base_entropies <- apply(votes, MARGIN=1, calc_scaled_entropy)
                  min_mislabels <- apply(votes, MARGIN=1, \(x) sum(x) - max(x))
                  
                  calc_swapped_delta_entropy <- function(swap_from_subject, swap_from_genotype,
                                                         swap_to_subject, swap_to_genotype) {
                      delta <- 0
                      if (!is.na(swap_from_genotype)) {
                          genotype_base_entropy <- base_entropies[swap_from_genotype]
                          genotype_votes_vec <- votes[swap_from_genotype, ]
                          genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[swap_from_subject] - 1
                          genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[swap_to_subject] + 1
                          genotype_new_entropy <- calc_scaled_entropy(genotype_votes_vec)
                          delta <- delta + genotype_new_entropy - genotype_base_entropy
                      }
                      if (!is.na(swap_to_genotype)) {
                          genotype_base_entropy <- base_entropies[swap_to_genotype]
                          genotype_votes_vec <- votes[swap_to_genotype, ]
                          genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[swap_to_subject] - 1
                          genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[swap_from_subject] + 1
                          genotype_new_entropy <- calc_scaled_entropy(genotype_votes_vec)
                          delta <- delta + genotype_new_entropy - genotype_base_entropy
                      }
                      return(delta)
                  }
                  
                  neighbors <- .find_neighbors(object, include_ghost, filter_concordant_vertices) %>% 
                      left_join(
                          unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")], 
                          by=c("Sample_A"="Sample_ID")
                      ) %>%
                      dplyr::rename(Subject_A = Subject_ID, Genotype_Group_A = Genotype_Group_ID) %>% 
                      left_join(
                          unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID", "Component_ID")], 
                          by=c("Sample_B"="Sample_ID")
                      ) %>% 
                      dplyr::rename(Subject_B = Subject_ID, Genotype_Group_B = Genotype_Group_ID) 
                  
                  print(glue("{nrow(neighbors)} total neighbors"))
                  
                  all_component_ids <- sort(unique(object@.solve_state$unsolved_relabel_data$Component_ID))
                  relabels_list <- vector("list", length=length(all_component_ids))
                  curr_idx <- 1
                  for (curr_component_id in all_component_ids) {
                      cc_relabel_data <- unsolved_all_data %>% filter(Component_ID == curr_component_id)
                      cc_neighbors <- neighbors %>% filter(Component_ID == curr_component_id)
                      cc_genotypes <- unique(cc_relabel_data$Genotype_Group_ID)
                      cc_genotypes <- cc_genotypes[!is.na(cc_genotypes)]
                      
                      if (nrow(cc_neighbors) == 0) {next}
                      
                      ## If you have too many neighbors to evaluate, downsample them
                      ## TODO: adjust min_mislabels based on how much it was downsampled
                      if (nrow(cc_neighbors) > MAX_NEIGHBORS_LOCAL_SEARCH) {
                          cc_neighbors <- cc_neighbors %>% sample_n(MAX_NEIGHBORS_LOCAL_SEARCH)
                      }
                      
                      cc_neighbor_objectives <- cc_neighbors %>% 
                          mutate(
                              delta = mapply(calc_swapped_delta_entropy, 
                                             swap_from_subject=Subject_A, 
                                             swap_from_genotype=Genotype_Group_A,
                                             swap_to_subject=Subject_B,
                                             swap_to_genotype=Genotype_Group_B)
                          )
                      
                      ## Find minimum number of mislabeled samples in the component
                      n_min_mislabels <- sum(min_mislabels[cc_genotypes])
                      
                      cc_candidate_neighbors <- cc_neighbor_objectives %>% 
                          filter(delta > 0) %>% 
                          arrange(desc(delta))
                      
                      ## Determine maximum number of swaps to look for 
                      ## If n_min_mislabels < 50, just make it a normal local search
                      n_swaps_limit <- max(1, floor(n_min_mislabels * frac_per_iter * 0.5))
                      if (n_min_mislabels <= 50) {
                          n_swaps_limit <- 1
                      }
                      cc_candidate_neighbors <- cc_candidate_neighbors %>% head(n_swaps_limit)
                      n_candidates <- nrow(cc_candidate_neighbors)
                      
                      if (n_candidates == 0) {next}
                      
                      print(glue("{n_min_mislabels} minimum mislabels in {curr_component_id}, searching through {n_candidates} candidate swaps"))
                      
                      cc_relabels_list <- vector("list", n_candidates)
                      for (j in 1:n_candidates) {
                          best_neighbor <- cc_candidate_neighbors[1, ]
                          best_neighbor_genotypes <- c(best_neighbor$Genotype_Group_A, best_neighbor$Genotype_Group_B)
                          cc_relabels_list[[j]] <- data.frame(relabel_from=best_neighbor$Sample_A, 
                                                              relabel_to=best_neighbor$Sample_B)
                          cc_candidate_neighbors <- cc_candidate_neighbors %>% 
                              filter(
                                  !(Genotype_Group_A %in% best_neighbor_genotypes),
                                  !(Genotype_Group_B %in% best_neighbor_genotypes)
                              )
                          n_candidates <- nrow(cc_candidate_neighbors)
                          if (n_candidates == 0) {break}
                          
                      }
                      cc_relabels <- do.call(rbind, cc_relabels_list) %>% 
                          filter(complete.cases(.))
                      print(glue("Found {nrow(cc_relabels)} swaps for {curr_component_id}"))
                      relabels_list[[curr_idx]] <- cc_relabels
                      curr_idx <- curr_idx + 1
                  }
                  relabels <- do.call(rbind, relabels_list)
                  relabels <- rbind(relabels, data.frame(relabel_from=relabels$relabel_to, relabel_to=relabels$relabel_from))
                  object <- .relabel_samples(object, relabels)
              }
              
              return(object)
          }
)

setMethod("solve", "MislabelSolver",
          function(object) {
              set.seed(1)
              while (TRUE) {
                  if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                      break
                  }
                  prev_relabel_data <- object@.solve_state$unsolved_relabel_data
                  
                  object <- solve_comprehensive_search(object)
                  object <- solve_majority_search(object)
                  object <- solve_comprehensive_search(object)
                  
                  comp_relabel_data <- object@.solve_state$unsolved_relabel_data
                  object <- solve_local_search(object, n_iter=1, include_ghost=TRUE, filter_concordant_vertices=TRUE)
                  
                  ## If local search found no swaps, try allowing concordant vertices
                  if (nrow(comp_relabel_data) == nrow(object@.solve_state$unsolved_relabel_data)) {
                      if (identical(comp_relabel_data, object@.solve_state$unsolved_relabel_data)) {
                          object <- solve_local_search(object, n_iter=1, include_ghost=TRUE, filter_concordant_vertices=FALSE)
                      }
                  }
                  rm(comp_relabel_data)
                  
                  if (nrow(prev_relabel_data) == nrow(object@.solve_state$unsolved_relabel_data)) {
                      if (identical(prev_relabel_data, object@.solve_state$unsolved_relabel_data)) {
                          break
                      }
                  }
                  rm(prev_relabel_data)
                  
                  object <- qdeserialize(qserialize(object))
                  gc()
              }
              
              return(object)
          }
)

setMethod("write_corrections", "MislabelSolver",
          function(object) {
              sample_summary <- object@.solve_state$relabel_data %>%
                  group_by(Genotype_Group_ID, Init_Subject_ID) %>% 
                  mutate(n_agree=n()) %>% 
                  ungroup(Init_Subject_ID) %>%
                  mutate(n_Genotype_Group=n()) %>% 
                  ungroup(Genotype_Group_ID) %>% 
                  left_join(
                      object@.solve_state$putative_subjects,
                      by = "Genotype_Group_ID",
                      suffix = c("", "_putative")
                  ) %>% 
                  transmute(
                      Component_ID = Init_Component_ID,
                      Genotype_Group_ID,
                      Ghost = is.na(Genotype_Group_ID),
                      Init_Subject_ID,
                      Init_Sample_ID,
                      Final_Subject_ID = Subject_ID,
                      Final_Sample_ID = Sample_ID,
                      Putative_Subject_ID = Subject_ID_putative,
                      n_Genotype_Group = ifelse(!Ghost, n_Genotype_Group, NA_integer_),
                      Init_Agreement = ifelse(!Ghost, paste0(n_agree-1, " out of ", n_Genotype_Group - 1), NA_character_),
                      Status = case_when(
                          Ghost ~ "ghost",
                          is.na(Final_Sample_ID) ~ "removed",
                          is.na(Putative_Subject_ID) ~ "subject_unknown",
                          Init_Sample_ID != Final_Sample_ID & (Final_Subject_ID != Putative_Subject_ID | n_Genotype_Group == 1) ~ "flagged",
                          Init_Sample_ID != Final_Sample_ID ~ "corrected",
                          n_Genotype_Group == 1 ~ "ignored_single-sample",
                          n_agree < 2 ~ "ignored",
                          TRUE ~ "validated"
                      )
                  ) %>% 
                  arrange(Component_ID, Genotype_Group_ID, Final_Subject_ID, Final_Sample_ID)
              
              # corrections_graph <- .generate_corrections_graph(object@.solve_state$relabel_data)
              # cluster_info <- components(corrections_graph)
              # cluster_sizes <- cluster_info$csize[cluster_info$membership]
              # names(cluster_sizes) <- names(cluster_info$membership)
              # sample_summary$Mislabel_Type <- NA_character_
              # sample_summary[sample_summary$Init_Sample_ID %in% names(cluster_sizes), "Mislabel_Type"] = as.character(cluster_sizes)
              # sample_summary %<>% mutate(
              #     Mislabel_Type = ifelse(
              #         is.na(Mislabel_Type),
              #         Mislabel_Type,
              #         ifelse(
              #             Mislabel_Type == "2", 
              #             "Swap",
              #             paste0(Mislabel_Type, "-cycle")
              #         ) 
              #     )
              # )
              
              ## TODO: fix this, hacky code to make abstract deadline
              ambiguities <- object@.solve_state$relabel_data %>% 
                  left_join(
                      object@.solve_state$putative_subjects,
                      by = "Genotype_Group_ID",
                      suffix = c("", "_putative")
                  ) %>% 
                  left_join(
                      object@swap_cats,
                      by = "Sample_ID"
                  ) %>% 
                  filter(Init_Subject_ID != Subject_ID_putative) %>% 
                  group_by(SwapCat_ID, Init_Subject_ID) %>% 
                  mutate(ambiguity_count = n()) %>% 
                  ungroup() %>% 
                  select(Final_Sample_ID = Init_Sample_ID, ambiguity_count)
              
              sample_summary <- sample_summary %>% 
                  left_join(
                      ambiguities,
                      by="Final_Sample_ID"
                  ) %>% 
                  mutate(
                      Status = ifelse(
                          ambiguity_count > 1,
                          "ambiguous_corrected",
                          Status
                      )
                  )
              
              # ambiguities <- sample_summary %>%
              #     filter(Status == "corrected") %>%
              #     select(Final_Subject_ID, Final_Sample_ID) %>%
              #     group_by(Final_Subject_ID) %>%
              #     mutate(n = n()) %>%
              #     filter(n > 1) %>%
              #     ungroup()
              # ambiguities$Ambiguities <- NULL
              # for (i in 1:nrow(ambiguities)) {
              #     sample_id <- ambiguities[i, "Final_Sample_ID"][[1]]
              #     subject_id <- ambiguities[i, "Final_Subject_ID"][[1]]
              #     sample_ambiguities <- ambiguities %>%
              #         filter(
              #             Final_Subject_ID == subject_id,
              #             Final_Sample_ID != sample_id
              #         ) %>%
              #         pull(Final_Sample_ID)
              #     ambiguities[i, "Ambiguities"] <- paste(sample_ambiguities, collapse=", ")
              # }
              # ambiguities %<>% select(Final_Sample_ID, Ambiguities)
              # sample_summary %<>% left_join(ambiguities, by="Final_Sample_ID")
              
              genotype_group_summary <- sample_summary %>% 
                  group_by(Genotype_Group_ID) %>% 
                  summarize(
                      Inferred_Subject_ID = names(sort(table(Final_Subject_ID), decreasing = TRUE)[1]),
                      n_Samples_validated = sum(Status == "validated"),
                      n_Samples_corrected = sum(Status == "corrected"),
                      n_Samples_removed = sum(Status == "removed"),
                      n_Samples_ignored = sum(str_detect(Status, "^ignored")),
                      n_Samples_total = n(),
                      ## TODO fix this more rigorously
                      Init_Fraction_Match = paste0(n_Samples_validated + n_Samples_ignored, " out of ", n_Samples_total)
                  ) %>% 
                  ungroup() %>% 
                  select(Genotype_Group_ID, Inferred_Subject_ID, Init_Fraction_Match, everything())
              
              component_summary <- sample_summary %>% 
                  group_by(Component_ID) %>% 
                  summarize(
                      n_Genotype_Groups = length(unique(Genotype_Group_ID)),
                      n_Subjects = length(unique(Final_Subject_ID)),
                      Size_Match = n_Genotype_Groups == n_Subjects,
                      Init_Solved = n_Genotype_Groups == 1 & n_Subjects == 1,
                      n_Samples_validated = sum(Status == "validated"),
                      n_Samples_corrected = sum(Status == "corrected"),
                      n_Samples_removed = sum(Status == "removed"),
                      n_Samples_ignored = sum(str_detect(Status, "^ignored")),
                      n_Samples_total = n()
                  )
              
              dataset_summary <- sample_summary %>%
                  summarize(
                      n_Components = length(unique(Component_ID)),
                      n_Genotype_Groups = length(unique(Genotype_Group_ID)),
                      n_Subjects = length(unique(Final_Subject_ID)),
                      n_Samples_validated = sum(Status == "validated"),
                      n_Samples_corrected = sum(Status == "corrected"),
                      n_Samples_removed= sum(Status == "removed"),
                      n_Samples_total = n()
                  )
              
              dir_name <- glue("corrections_output_", str_replace_all(Sys.time(), "[:. ]", "_"))
              dir.create(dir_name)
              
              excel_filename <- file.path(dir_name, "corrections_summary.xlsx")
              summary_list <- list(
                  "Sample" = sample_summary,
                  "Genotype_Group" = genotype_group_summary,
                  "Component" = component_summary,
                  "Dataset" = dataset_summary)
              write.xlsx(summary_list, file=excel_filename)
              
              # component_data <- rbind(unsolved_relabel_data, unsolved_ghost_data) %>% 
              #     group_by(Component_ID) %>% 
              #     summarize(
              #         n_Genotype_Group_ID = n_distinct(Genotype_Group_ID) - anyNA(Genotype_Group_ID),
              #         n_Subject_ID = n_distinct(Subject_ID),
              #         n_Sample_ID = length(Sample_ID)
              #     ) %>%
              #     mutate(Solved = n_Genotype_Group_ID <= 1 & n_Subject_ID == 1)
              # solved_components <- component_data %>% filter(Solved) %>% pull(Component_ID)
              
              init_solved_components <- object@.solve_state$relabel_data %>%
                  group_by(Init_Component_ID) %>%
                  summarize(
                      n_Genotype_Group_ID = n_distinct(Genotype_Group_ID) - anyNA(Genotype_Group_ID),
                      n_Subject_ID = n_distinct(Subject_ID),
                      n_Sample_ID = length(Sample_ID)
                  ) %>%
                  mutate(Solved = n_Genotype_Group_ID <= 1 & n_Subject_ID == 1) %>% 
                  filter(Solved) %>% pull(Init_Component_ID)
              unsolved_components <- sort(setdiff(object@.solve_state$relabel_data$Init_Component_ID, init_solved_components))
              
              genotyped_relabel_data <- object@.solve_state$relabel_data %>% filter(!is.na(Genotype_Group_ID))
              ghost_relabel_data <- object@.solve_state$relabel_data %>% filter(is.na(Genotype_Group_ID))
              for (component in unsolved_components) {
                  save_file <- file.path(dir_name, glue(component, ".png"))
                  png(filename=save_file, width=1600, height=1600)
                  layout_matrix <- matrix(1:4, ncol=2)
                  layout(layout_matrix)
                  .plot_graph(.generate_graph(genotyped_relabel_data %>% 
                                                  filter(Init_Component_ID == component) %>% 
                                                  mutate(Sample_ID = Init_Sample_ID, Subject_ID = Init_Subject_ID), 
                                              graph_type="combined",
                                              ghost_relabel_data %>% 
                                                  filter(Init_Component_ID == component) %>% 
                                                  mutate(Sample_ID = Init_Sample_ID, Subject_ID = Init_Subject_ID), 
                                              object@anchor_samples,
                                              object@swap_cats), to_write=TRUE)
                  title(main=glue(component, ": Before corrections"), cex.main=2)
                  .plot_graph(.generate_corrections_graph(object@.solve_state$relabel_data %>% filter(Init_Component_ID == component)), to_write=TRUE)
                  title(main=glue(component, ": Applied corrections"), cex.main=2)
                  .plot_graph(.generate_graph(genotyped_relabel_data %>% 
                                                  filter(Init_Component_ID == component),
                                              graph_type="combined",
                                              ghost_relabel_data %>% 
                                                  filter(Init_Component_ID == component),
                                              object@anchor_samples,
                                              object@swap_cats), to_write=TRUE)
                  title(main=glue(component, ": After corrections"), cex.main=2)
                  dev.off()
              }
          }
)


setMethod("plot", "MislabelSolver",
          function(x, 
                   y=NULL, 
                   unsolved=TRUE, 
                   query_by=c("Init_Component_ID", "Component_ID", "Subject_ID", "Genotype_Group_ID", "Sample_ID"),
                   query_val=NULL) {
              if (unsolved) {
                  relabel_data <- x@.solve_state$unsolved_relabel_data
                  ghost_data <- x@.solve_state$unsolved_ghost_data
              } else {
                  relabel_data <- x@.solve_state$relabel_data %>% filter(!is.na(Genotype_Group_ID))
                  ghost_data <- x@.solve_state$relabel_data %>% filter(is.na(Genotype_Group_ID))
              }
              anchor_samples <- x@anchor_samples
              if (!is.null(query_val)) {
                  query_by <- as.character(query_by)
                  query_by <- match.arg(query_by)
                  if (query_by == "Init_Component_ID") {
                      component_id <- rbind(relabel_data, ghost_data) %>% 
                          filter(!!sym(query_by) == query_val) %>% 
                          pull(Init_Component_ID) %>% 
                          unique()
                  } else {
                      component_id <- rbind(relabel_data, ghost_data) %>% 
                          filter(!!sym(query_by) == query_val) %>% 
                          pull(Component_ID) %>% 
                          unique()
                  }
                  if (length(component_id) == 0) {
                      warning(glue("No samples found for 'query_by' \"{query_by}\" and 'query_val' \"{query_val}\""))
                      return()
                  }
                  component_id <- component_id[[1]]
                  if (query_by == "Init_Component_ID") {
                      relabel_data <- relabel_data %>% filter(Init_Component_ID == component_id)
                      ghost_data <- ghost_data %>% filter(Init_Component_ID == component_id)
                  } else {
                      relabel_data <- relabel_data %>% filter(Component_ID == component_id)
                      ghost_data <- ghost_data %>% filter(Component_ID == component_id)
                  }
              }
              graph <- .generate_graph(relabel_data, graph_type = "combined", ghost_data, anchor_samples, populate_plotting_attributes=TRUE)
              with_seed(1, {
                  l_mds <- layout_with_mds(graph)
                  l_drl <- layout_with_drl(graph, use.seed=TRUE, seed=l_mds)
                  visIgraph(graph, layout = "layout_with_graphopt", start=l_drl) 
              })
          }
)


# setGeneric("plot_corrections", 
#            function(x, 
#                     y=NULL, 
#                     query_by=c("Component_ID", "Subject_ID", "Genotype_Group_ID", "Sample_ID"), 
#                     query_val=NULL) {
#                standardGeneric("plot_corrections")
#            })

# setMethod("plot_corrections", "MislabelSolver",
#           function(x,
#                    y=NULL,
#                    query_by=c("Init_Component_ID", "Component_ID", "Subject_ID", "Genotype_Group_ID", "Sample_ID"),
#                    query_val=NULL) {
#               relabel_data <- x@.solve_state$relabel_data
#               swap_cats <- x@swap_cats
#               if (!is.null(query_val)) {
#                   query_by <- as.character(query_by)
#                   query_by <- match.arg(query_by)
#                   if (query_by == "Component_ID") {query_by <- "Init_Component_ID"}
#                   component_id <- relabel_data %>%
#                       filter(!!sym(query_by) == query_val) %>%
#                       pull(Init_Component_ID) %>%
#                       unique()
#                   if (length(component_id) == 0) {
#                       warning(glue("No samples found for 'query_by' \"{query_by}\" and 'query_val' \"{query_val}\""))
#                       return()
#                   }
#                   component_id <- component_id[[1]]
#                   relabel_data <- relabel_data %>% filter(Init_Component_ID == component_id)
#               }
#               graph <- .generate_corrections_graph(relabel_data)
#               #graph <- .generate_corrections_graph(relabel_data, swap_cats, populate_plotting_attributes=TRUE)
#               with_seed(1, visIgraph(graph, start=l_drl))
#               with_seed(1, visIgraph(graph, layout="layout_nicely"))
#           }
# )

setMethod("summary", "MislabelSolver",
          function(object) {
              print("MislabelSolver summary statistics")
              cat("\t", "Number of Sample_ID(s): ", length(object@sample_genotype_data$Sample_ID))
              cat("\t", "test")
          }
) 

################################################################################
##################               SOLVE HELPERS                ##################
################################################################################

.update_solve_state <- function(object, initialization=FALSE) {
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        return(object)
    }
    
    ## 1. Assign a Component_ID for each unsolved Sample_ID
    combined_graph <- .generate_graph(object@.solve_state$unsolved_relabel_data, graph_type="combined", object@.solve_state$unsolved_ghost_data)
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data %>% 
        mutate(
            Component_ID = as.character(components(combined_graph)$membership[Sample_ID])
        )
    unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data %>% 
        mutate(
            Component_ID = as.character(components(combined_graph)$membership[Sample_ID])
        )
    
    ## 2. Determine which Sample_ID(s) are newly solved
    ##    A Sample_ID is solved if it belongs to a Component_ID that includes only one
    ##    Genotype_Group_ID and one Subject_ID
    component_data <- rbind(unsolved_relabel_data, unsolved_ghost_data) %>% 
        group_by(Component_ID) %>% 
        summarize(
            n_Genotype_Group_ID = n_distinct(Genotype_Group_ID) - anyNA(Genotype_Group_ID),
            n_Subject_ID = n_distinct(Subject_ID),
            n_Sample_ID = length(Sample_ID)
        ) %>%
        mutate(Solved = n_Genotype_Group_ID <= 1 & n_Subject_ID == 1)
    
    solved_components <- component_data %>% filter(Solved) %>% pull(Component_ID)
    
    unsolved_relabel_data <- unsolved_relabel_data %>% 
        mutate(Solved = Component_ID %in% solved_components)
    unsolved_ghost_data <- unsolved_ghost_data %>% 
        mutate(Solved = Component_ID %in% solved_components)
    
    ## 3. Re-rank Component_ID(s) in order of size (so that Component1 is the largest unsolved component)
    # component_data <- component_data %>% filter(!Solved)
    n_components <- nrow(component_data)
    component_data <- component_data %>% 
        arrange(Solved, desc(n_Sample_ID)) %>% 
        mutate(
            new_Component_ID = if (n_components != 0) 1:n_components else character(0),
            new_Component_ID = paste0("Component", formatC(new_Component_ID, width=nchar(n_components), format="d", flag="0"))
        )
    unsolved_relabel_data <- unsolved_relabel_data %>% 
        left_join(component_data[, c("Component_ID", "new_Component_ID")], by="Component_ID") %>% 
        mutate(Component_ID = new_Component_ID) %>% 
        select(-new_Component_ID)
    unsolved_ghost_data <- unsolved_ghost_data %>% 
        left_join(component_data[, c("Component_ID", "new_Component_ID")], by="Component_ID") %>% 
        mutate(Component_ID = new_Component_ID) %>% 
        select(-new_Component_ID)
    
    ## 4. Update putative_subjects
    ##    During initialization, lock Subject_ID/Genotype_Group_ID pairs for anchor_samples
    if (initialization) {
        anchor_putative_subjects <- object@sample_genotype_data %>% 
            filter(!is.na(Genotype_Group_ID),
                   !is.na(Subject_ID),
                   Sample_ID %in% object@anchor_samples) %>% 
            select(Subject_ID, Genotype_Group_ID) %>% 
            unique()
        object <- .update_putative_subjects(object, anchor_putative_subjects)
    }
    solved_putative_subjects <- unsolved_relabel_data %>% 
        filter(Solved) %>% 
        select(Genotype_Group_ID, Subject_ID) %>% 
        distinct()
    object <- .update_putative_subjects(object, solved_putative_subjects)
    
    ## 5. Update relabel_data, and unsolved_relabel_data
    if (initialization) {
        unsolved_relabel_data <- unsolved_relabel_data %>% 
            select(Init_Sample_ID, Init_Subject_ID, Genotype_Group_ID, Component_ID, Sample_ID, Subject_ID, Solved, SwapCat_ID, SwapCat_Shape) %>% 
            mutate(Init_Component_ID = Component_ID) %>% 
            relocate(Init_Component_ID, .before=Component_ID)
        unsolved_ghost_data <- unsolved_ghost_data %>% 
            select(Init_Sample_ID, Init_Subject_ID, Genotype_Group_ID, Component_ID, Sample_ID, Subject_ID, Solved, SwapCat_ID, SwapCat_Shape) %>% 
            mutate(Init_Component_ID = Component_ID) %>% 
            relocate(Init_Component_ID, .before=Component_ID)
        relabel_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
    } else {
        ## Update 'relabel_data' with new sample labels in 'unsolved_relabel_data' and 'unsolved_ghost_data'
        unsolved_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
        relabel_data <- object@.solve_state$relabel_data %>% 
            left_join(unsolved_data[, c("Component_ID", "Sample_ID", "Subject_ID", "Init_Sample_ID", "Solved", "SwapCat_ID", "SwapCat_Shape")], 
                      by="Init_Sample_ID", suffix = c(".x", ".y")) %>% 
            mutate(
                Component_ID = coalesce(Component_ID.y, Component_ID.x),
                Sample_ID = coalesce(Sample_ID.y, Sample_ID.x),
                Subject_ID = coalesce(Subject_ID.y, Subject_ID.x),
                Solved = coalesce(Solved.y, Solved.x),
                SwapCat_ID = coalesce(SwapCat_ID.y, SwapCat_ID.x),
                SwapCat_Shape = coalesce(SwapCat_Shape.y, SwapCat_Shape.x)
            ) %>% 
            select(-ends_with(c(".x", "y"))) %>% 
            select(Init_Sample_ID, Init_Subject_ID, Genotype_Group_ID, Init_Component_ID, Component_ID, Sample_ID, Subject_ID, Solved, SwapCat_ID, SwapCat_Shape)
    }
    unsolved_relabel_data <- unsolved_relabel_data %>% filter(!Solved)
    unsolved_ghost_data <- unsolved_ghost_data %>% filter(!Solved)
    
    ## 6. Overwrite .solve_state
    object@.solve_state$relabel_data <- relabel_data
    object@.solve_state$unsolved_relabel_data <- unsolved_relabel_data
    object@.solve_state$unsolved_ghost_data <- unsolved_ghost_data
    
    return(object)
}

.update_putative_subjects <- function(object, proposed_putative_subjects) {
    if (nrow(proposed_putative_subjects) == 0) return(object)
    ## Only add Genotype_Group_ID/Subject_ID combinations if neither the 
    ## Genotype_Group_ID nor the Subject_ID are already in putative_subjects
    existing_genotypes <- na.omit(object@.solve_state$putative_subjects$Genotype_Group_ID)
    existing_subjects <- na.omit(object@.solve_state$putative_subjects$Subject_ID)
    proposed_putative_subjects <- proposed_putative_subjects %>% 
        filter(
            !(Subject_ID %in% existing_subjects),
            !(Genotype_Group_ID %in% existing_genotypes)
        )
    putative_subjects <- rbind(object@.solve_state$putative_subjects, proposed_putative_subjects)
    .validate_putative_subjects(object@sample_genotype_data, putative_subjects)
    object@.solve_state$putative_subjects <- putative_subjects
    return(object)
}

.relabel_samples <- function(object, relabels) {
    if (nrow(relabels) == 0) return(object)
    
    ## Ensure that you can only relabel unsolved samples
    # TODO fix validation
    # .validate_relabels(object@.solve_state$unsolved_relabel_data, object@.solve_state$unsolved_ghost_data, relabels)
    relabels <- relabels %>% 
        rename("Sample_ID" = "relabel_to") %>% 
        left_join(
            object@sample_genotype_data[, c("Sample_ID", "Subject_ID")],
            by="Sample_ID"
        ) %>% mutate(
            Subject_ID = ifelse(
                is.na(Subject_ID),
                sapply(Sample_ID, \(x) unlist(strsplit(x, split=":"))[2]),
                Subject_ID
            )
        )
    ## Call it relabeled sample ID instead
    unsolved_all_data <- rbind(object@.solve_state$unsolved_relabel_data, 
                               object@.solve_state$unsolved_ghost_data)
    unsolved_all_data <- unsolved_all_data %>% 
        left_join(
            relabels, by=c("Sample_ID"="relabel_from"), suffix=c(".x", ".y")
        ) %>% 
        mutate(
            Sample_ID = ifelse(!is.na(Sample_ID.y), Sample_ID.y, Sample_ID),
            Subject_ID = ifelse(!is.na(Subject_ID.y), Subject_ID.y, Subject_ID.x)
        ) %>% 
        select(-ends_with(".x"), -ends_with(".y"))
    
    object@.solve_state$unsolved_relabel_data <- unsolved_all_data %>% filter(!is.na(Genotype_Group_ID))
    object@.solve_state$unsolved_ghost_data <- unsolved_all_data %>% filter(is.na(Genotype_Group_ID))
    object <- .update_solve_state(object)
    return(object)
}

.find_neighbors <- function(object, include_ghost=FALSE, filter_concordant_vertices=FALSE) {
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    unsolved_ghost_data <- NULL
    if (include_ghost) {
        unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
    }
    combined_graph <- .generate_graph(unsolved_relabel_data, graph_type="combined", unsolved_ghost_data)
    unsolved_all_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
    putative_subjects <- object@.solve_state$putative_subjects
    
    ## Criteria 1:  filter only pairs of vertices that are within at exactly 2 edges of each other
    adj_matrix_sparse <- Matrix::Matrix(get.adjacency(combined_graph, sparse = TRUE))
    adj_matrix_idx <- which(adj_matrix_sparse > 0, arr.ind = TRUE)
    dist_within_2_sparse <- adj_matrix_sparse %*% adj_matrix_sparse
    dist_within_2_idx <- which(dist_within_2_sparse > 0, arr.ind = TRUE)
    unique_pairs_exact1 <- data.frame(
        Row = rownames(adj_matrix_sparse)[adj_matrix_idx[, 1]],
        Col = colnames(adj_matrix_sparse)[adj_matrix_idx[, 2]]) %>%
        transmute(
            Sample_A = pmin(Row, Col),
            Sample_B = pmax(Row, Col)
        ) %>% 
        filter(Sample_A != Sample_B) %>% 
        unique()
    unique_pairs_within2 <- data.frame(
        Row = rownames(dist_within_2_sparse)[dist_within_2_idx[, 1]],
        Col = colnames(dist_within_2_sparse)[dist_within_2_idx[, 2]]) %>%
        transmute(
            Sample_A = pmin(Row, Col),
            Sample_B = pmax(Row, Col)
        ) %>% 
        filter(Sample_A != Sample_B) %>% 
        unique()
    unique_pairs <- anti_join(unique_pairs_within2, unique_pairs_exact1, by=c("Sample_A", "Sample_B"))
    
    ## Criteria 2: filter out pairs of vertices that include elements in anchor_samples
    unique_pairs <- unique_pairs[!(unique_pairs$Sample_A %in% object@anchor_samples |unique_pairs$Sample_B %in% object@anchor_samples), ]
    
    ## Criteria 3: filter only pairs of vertices that are within the same swap category
    unique_pairs <- unique_pairs %>% 
        left_join(unsolved_all_data[, c("Sample_ID", "SwapCat_ID")], by=c("Sample_A"="Sample_ID")) %>% 
        dplyr::rename("SwapCat_A"="SwapCat_ID") %>% 
        left_join(unsolved_all_data[, c("Sample_ID", "SwapCat_ID")], by=c("Sample_B"="Sample_ID")) %>% 
        dplyr::rename("SwapCat_B"="SwapCat_ID") %>% 
        filter(SwapCat_A == SwapCat_B) %>% 
        select(Sample_A, Sample_B)
    
    ## Criteria 4: filter out vertices that have at least 1 concordant edge
    if (filter_concordant_vertices) {
        concordant_edges <- E(combined_graph)[E(combined_graph)$concordant]
        concordant_vertices <- unique(c(ends(combined_graph, concordant_edges)))
        unique_pairs <- unique_pairs %>%
            filter(
                !(Sample_A %in% concordant_vertices),
                !(Sample_B %in% concordant_vertices)
            )
    }
    
    ## Criteria 5: filter out pairs of vertices where both sides will violate putative_subjects
    unique_pairs <- unique_pairs %>%
        left_join(
            unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")],
            by=c("Sample_A"="Sample_ID")) %>% 
        rename("Subject_A"="Subject_ID", "Genotype_Group_A"="Genotype_Group_ID") %>% 
        left_join(
            putative_subjects %>% filter(!is.na(Genotype_Group_ID)),
            by=c("Genotype_Group_A"="Genotype_Group_ID")) %>% 
        rename("Putative_Subject_A"="Subject_ID") %>%
        left_join(
            unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")],
            by=c("Sample_B"="Sample_ID")) %>% 
        rename("Subject_B"="Subject_ID", "Genotype_Group_B"="Genotype_Group_ID") %>% 
        left_join(
            putative_subjects %>% filter(!is.na(Genotype_Group_ID)),
            by=c("Genotype_Group_B"="Genotype_Group_ID")) %>% 
        rename("Putative_Subject_B"="Subject_ID") %>% 
        mutate(
            # The swap is invalid if
            # 1. Sample_A is already in its Putative_Subject
            # 2. Sample_B is already in its Putative_Subject
            # For swaps where the putative subjects of both samples are assigned, the swap is also invalid if
            # 3. Neither Sample_A nor Sample_B match their Putative_Subject after
            Invalid_Swap = 
                (!is.na(Putative_Subject_A) & Subject_A == Putative_Subject_A) |
                (!is.na(Putative_Subject_B) & Subject_B == Putative_Subject_B) | 
                (!is.na(Putative_Subject_A) & !is.na(Putative_Subject_B) & Putative_Subject_A != Subject_B & Putative_Subject_B != Subject_A)
        ) %>%
        filter(!Invalid_Swap) %>%
        select(Sample_A, Sample_B)
    
    return(unique_pairs)
}

## NOTE: realized that "allow_unknowns" and "allow_ghosts" only works when all samples in the component have
## been genotyped. That's fine for now since we only call this with these parameters
## after comprehensive search but it's not ideal, go back and fix this logic
.find_relabel_cycles_from_putative_subjects <- function(unsolved_relabel_data, 
                                                        putative_subjects,
                                                        unsolved_ghost_data=NULL, 
                                                        unambiguous_only=FALSE, 
                                                        allow_unknowns=FALSE) {
    allow_ghosts <- !is.null(unsolved_ghost_data)
    
    if (unambiguous_only) {
        allow_ghosts <- FALSE
        allow_unknowns <- FALSE
    }
    # 
    # unsolved_relabel_data <- unsolved_relabel_data %>% 
    #     left_join(swap_cats, by="Sample_ID")
    # 
    # if (allow_ghosts) {
    #     unsolved_ghost_data <- unsolved_ghost_data %>% 
    #         left_join(swap_cats, by="Sample_ID")
    # }
    # 
    mislabel_data <- unsolved_relabel_data %>%
        left_join(putative_subjects %>% rename(Putative_Subject_ID = Subject_ID), 
                  by="Genotype_Group_ID") %>% 
        mutate(
            Inferred_Correctly_Labeled = Putative_Subject_ID == Subject_ID,
            ## Only include samples where the current label is a Subject_ID that we've found
            ## the right Genotype_Group_ID for
            Curr_Subject_ID_Genotyped = Subject_ID %in% putative_subjects$Subject_ID
        ) %>% 
        filter(!is.na(Putative_Subject_ID) & !Inferred_Correctly_Labeled & Curr_Subject_ID_Genotyped)
    
    relabels <- EMPTY_RELABELS
    
    if (nrow(mislabel_data) == 0) {
        return(relabels)
    }
    
    mislabeled_genotype_swapcats <- mislabel_data %>% 
        select(SwapCat_ID, Genotype_Group_ID) %>% 
        unique() %>% 
        left_join(putative_subjects, by="Genotype_Group_ID")
    
    directed_edge_mats <- vector("list", length=nrow(mislabeled_genotype_swapcats))
    
    all_ghost_labels <- character(0)
    all_unknown_labels <- character(0)
    
    for (i in 1:nrow(mislabeled_genotype_swapcats)) {
        swap_cat_id <- mislabeled_genotype_swapcats[i, "SwapCat_ID"]
        genotype_group_id <- mislabeled_genotype_swapcats[i, "Genotype_Group_ID"]
        subject_id <- mislabeled_genotype_swapcats[i, "Subject_ID"]
        mislabeled_samples <- mislabel_data %>% 
            filter(SwapCat_ID == swap_cat_id, Genotype_Group_ID == genotype_group_id) %>% 
            pull(Sample_ID)
        n_mislabeled_samples <- length(mislabeled_samples)
        genotyped_labels <- unsolved_relabel_data %>% 
            filter(SwapCat_ID == swap_cat_id, Subject_ID == subject_id, Genotype_Group_ID != genotype_group_id) %>% 
            pull(Sample_ID)
        n_genotyped_labels <- length(genotyped_labels)
        ghost_labels <- character(0)
        n_ghost_labels <- 0
        if (allow_ghosts) {
            ghost_labels <- unsolved_ghost_data %>% 
                filter(SwapCat_ID == swap_cat_id, Subject_ID == subject_id) %>% 
                pull(Sample_ID)
            n_ghost_labels <- length(ghost_labels)
            ## If we already have enough labels, no need to use ghosts
            if (n_mislabeled_samples - n_genotyped_labels <= 0) {
                n_ghost_labels <- 0
                ghost_labels <- character(0)
            }
            ## If we have more ghost labels than needed, select a subset
            else if (n_mislabeled_samples - n_genotyped_labels < length(ghost_labels)) {
                n_ghost_labels <- n_mislabeled_samples - n_genotyped_labels
                ghost_labels <- ghost_labels[1:n_ghost_labels]
            }
        }
        unknown_labels <- character(0)
        n_unknowns <- 0
        if (allow_unknowns) {
            unknown_labels <- character(0)
            n_unknown_labels <- n_mislabeled_samples - n_genotyped_labels - n_ghost_labels
            if (n_unknown_labels > 0) {
                unknown_labels <- 1:n_unknown_labels
                unknown_labels <- paste0(LABEL_NOT_FOUND, ":", subject_id, ":", swap_cat_id, ":", unknown_labels)
            }
        }
        directed_edge_mats[[i]] <- expand.grid(relabel_from=mislabeled_samples, relabel_to=c(genotyped_labels, ghost_labels, unknown_labels))
        all_ghost_labels <- c(all_ghost_labels, ghost_labels)
        all_unknown_labels <- c(all_unknown_labels, unknown_labels)
    }
    
    directed_edge_df <- as.data.frame(do.call(rbind, directed_edge_mats))
    relabels_graph <- graph_from_data_frame(directed_edge_df, directed=TRUE)
    
    ## If unambiguous_only, only include samples with exactly one incoming and one outgoing edge
    if (unambiguous_only) {
        samples_with_one_incoming <- V(relabels_graph)[degree(relabels_graph, mode = "in") == 1]$name
        samples_with_one_outgoing <- V(relabels_graph)[degree(relabels_graph, mode = "out") == 1]$name
        samples_to_filter <- intersect(samples_with_one_incoming, samples_with_one_outgoing)
        relabels_graph <- subgraph(relabels_graph, samples_to_filter)
    }
    
    ## 1. Find relabel cycles without using ghosts or unknowns
    new_relabels <- .find_all_relabel_cycles(relabels_graph)
    relabels_graph <- relabels_graph - new_relabels$relabel_from
    relabels <- rbind(relabels, new_relabels)
    
    ## 2. Find relabel cycles allowing for ghosts
    if (allow_ghosts) {
        all_excess_labels <- names(V(relabels_graph)[degree(relabels_graph, mode = "in") == 0])
        V(relabels_graph)$relabel_component_id <- components(relabels_graph)$membership
        all_relabel_component_ids <- unique(V(relabels_graph)$relabel_component_id)
        ## For each component, connect ghost labels with labels that have no more incoming edges, 
        ## indicating that there are no genotyped samples that can take their label
        for (curr_relabel_component_id in all_relabel_component_ids) {
            component_labels <- names(V(relabels_graph)[V(relabels_graph)$relabel_component_id == curr_relabel_component_id])
            component_ghost_labels <- intersect(all_ghost_labels, component_labels)
            component_excess_labels <- intersect(all_excess_labels, component_labels)
            new_edges <- expand.grid(relabel_from=component_ghost_labels, relabel_to=component_excess_labels)
            if (nrow(new_edges) > 0) {
                new_relabels_graph <- graph_from_data_frame(new_edges)
                relabels_graph <- union(relabels_graph, new_relabels_graph)   
            }
        }
        new_relabels <- .find_all_relabel_cycles(relabels_graph)
        relabels_graph <- relabels_graph - new_relabels$relabel_from
        relabels <- rbind(relabels, new_relabels)
    }
    
    ## 3. Find relabel cycles allowing for unknowns
    if (allow_unknowns) {
        all_excess_labels <- names(V(relabels_graph)[degree(relabels_graph, mode = "in") == 0])
        V(relabels_graph)$relabel_component_id <- components(relabels_graph)$membership
        all_relabel_component_ids <- unique(V(relabels_graph)$relabel_component_id)
        ## For each component, connect unknown labels with labels that have no more incoming edges, 
        ## indicating that there are no genotyped samples that can take their label
        for (curr_relabel_component_id in all_relabel_component_ids) {
            component_labels <- names(V(relabels_graph)[V(relabels_graph)$relabel_component_id == curr_relabel_component_id])
            component_unknown_labels <- intersect(all_unknown_labels, component_labels)
            component_excess_labels <- intersect(all_excess_labels, component_labels)
            new_edges <- expand.grid(relabel_from=component_unknown_labels, relabel_to=component_excess_labels)
            if (nrow(new_edges) > 0) {
                new_relabels_graph <- graph_from_data_frame(new_edges)
                relabels_graph <- union(relabels_graph, new_relabels_graph)    
            }
        }
        new_relabels <- .find_all_relabel_cycles(relabels_graph)
        relabels_graph <- relabels_graph - new_relabels$relabel_from
        relabels <- rbind(relabels, new_relabels)
    }
    
    return(relabels)
}

################################################################################
##################            VALIDATION HELPERS              ##################
################################################################################

.validate_sample_genotype_data <- function(sample_genotype_data) {
    required_columns <- c("Sample_ID", "Subject_ID", "Genotype_Group_ID")
    missing_columns <- setdiff(required_columns, colnames(sample_genotype_data))
    assert_that(
        length(missing_columns) == 0,
        msg = glue("'sample_genotype_data' is missing required column(s) {paste(missing_columns, collapse=\", \")}")
    )
    
    duplicated_samples <- sample_genotype_data$Sample_ID[duplicated(sample_genotype_data$Sample_ID)]
    assert_that(
        length(duplicated_samples) == 0,
        msg = glue("'sample_genotype_data' has non-unique 'Sample_ID'(s) {paste(duplicated_samples, collapse=\", \")}")
    )
    
    assert_that(
        all(!is.na(sample_genotype_data$Sample_ID)),
        msg = "'sample_genotype_data' missing value(s) in 'Sample_ID' column"
    )
    
    assert_that(
        all(!is.na(sample_genotype_data$Subject_ID)),
        msg = "'sample_genotype_data' missing value(s) in 'Subject_ID' column"
    )
}

.validate_swap_cats <- function(sample_genotype_data, swap_cats) {
    assert_that(
        colnames(swap_cats)[[1]] == "Sample_ID",
        msg = "'swap_cats' first column must have name 'Sample_ID'"
    )
    
    assert_that(
        colnames(swap_cats)[[2]] == "SwapCat_ID",
        msg = "'swap_cats' second column must have name 'SwapCat_ID'"
    )
    
    missing_values <- swap_cats[is.na(swap_cats$SwapCat_ID), "Sample_ID"]
    assert_that(
        length(missing_values) == 0,
        msg = glue("'swap_cats' is missing SwapCat_ID for Sample_ID(s) {paste(missing_values, collapse=\", \")}")
    )
    
    duplicated_samples <- swap_cats$Sample_ID[duplicated(swap_cats$Sample_ID)]
    assert_that(
        length(duplicated_samples) == 0,
        msg = glue("'swap_cats' has non-unique Sample_ID(s) {paste(duplicated_samples, collapse=\", \")}")
    )
    
    missing_samples <- setdiff(sample_genotype_data$Sample_ID, swap_cats$Sample_ID)
    assert_that(
        length(missing_samples) == 0,
        msg = glue("'swap_cats' is missing Sample_ID(s) that exist in 'sample_genotype_data', check Sample_ID(s) {paste(missing_samples, collapse=\", \")}")
    )
}

.validate_putative_subjects <- function(sample_genotype_data, putative_subjects) {
    required_columns <- c("Genotype_Group_ID", "Subject_ID")
    missing_columns <- setdiff(required_columns, colnames(putative_subjects))
    assert_that(
        length(missing_columns) == 0,
        msg = glue("'putative_subjects' is missing required column(s) {paste(missing_columns, collapse=\", \")}")
    )
    
    extra_genotype_groups <- setdiff(na.omit(putative_subjects$Genotype_Group_ID), sample_genotype_data$Genotype_Group_ID)
    assert_that(
        length(extra_genotype_groups) == 0,
        msg = glue("'putative_subjects' has 'Genotype_Group_ID'(s) not found in 'sample_genotype_data', check {paste(extra_genotype_groups, collapse=\", \")}")
    )
    extra_subjects <- setdiff(na.omit(putative_subjects$Subject_ID), sample_genotype_data$Subject_ID)
    assert_that(
        length(extra_subjects) == 0,
        msg = glue("'putative_subjects' has 'Subject_ID'(s) not found in 'sample_genotype_data', check {paste(extra_subjects, collapse=\", \")}")
    )
    
    duplicated_genotype_groups <- putative_subjects$Genotype_Group_ID[duplicated(na.omit(putative_subjects$Genotype_Group_ID))]
    assert_that(
        length(duplicated_genotype_groups) == 0,
        msg = glue("'putative_subjects' does not map 'Subject_ID' to 'Genotype_Group_ID' one-to-one, check 'Genotype_Group'(s) {paste(duplicated_genotype_groups, collapse=\", \")}")
    )
    duplicated_subjects <- putative_subjects$Subject_ID[duplicated(na.omit(putative_subjects$Subject_ID))]
    assert_that(
        length(duplicated_subjects) == 0,
        msg = glue("'putative_subjects' does not map 'Subject_ID' to 'Genotype_Group_ID' one-to-one, 'check Subject_ID'(s) {paste(duplicated_samples, collapse=\", \")}")
    )
}

.validate_relabels <- function(sample_genotype_data, relabels) {
    required_columns <- c("relabel_from", "relabel_to")
    missing_columns <- setdiff(required_columns, names(relabels))
    assert_that(
        length(missing_columns) == 0,
        msg = glue("'relabels' is missing required column(s) {paste(missing_columns, collapse=\", \")}")
    )
    
    missing_relabel_from <- setdiff(relabels$relabel_from, sample_genotype_data$Sample_ID)
    assert_that(
        length(missing_relabel_from) == 0,
        msg = glue("'relabels$relabel_from' has Sample_ID(s) not in 'sample_genotype_data', check {paste(missing_relabel_from, collapse=\", \")}")
    )
    missing_relabel_to <- setdiff(relabels$relabel_to, sample_genotype_data$Sample_ID)
    assert_that(
        length(missing_relabel_to) == 0,
        msg = glue("'relabels$relabel_to' has Sample_ID(s) not in 'sample_genotype_data', check {paste(missing_relabel_to, collapse=\", \")}")
    )
    
    duplicated_relabel_from <- relabels$relabel_from[duplicated(relabels$relabel_from)]
    assert_that(
        length(duplicated_relabel_from) == 0,
        msg = glue("'relabels$relabel_from' has non-unique Sample_ID(s), check Sample_ID(s) {paste(duplicated_relabel_from, collapse=\", \")}")
    )
    duplicated_relabel_to <- relabels$relabel_to[duplicated(relabels$relabel_to)]
    assert_that(
        length(duplicated_relabel_to) == 0,
        msg = glue("'relabels$relabel_to' has non-unique Sample_ID(s), check Sample_ID(s) {paste(duplicated_genotype_groups, collapse=\", \")}")
    )
    
    creating_duplicates <- setdiff(relabels$relabel_to, relabels$relabel_from)
    assert_that(
        length(creating_duplicates) == 0,
        msg = glue("'relabels$relabel_to' contains Sample_ID(s) that don't exist in 'relabels$relabel_from', check Sample_ID(s) {paste(creating_duplicates, collapse=\", \")}")
    )
}

.validate_anchor_samples <- function(sample_genotype_data, anchor_samples) {
    extra_samples <- setdiff(anchor_samples, sample_genotype_data$Sample_ID)
    assert_that(
        length(extra_samples) == 0,
        msg = glue("'anchor_samples' contains Sample_ID(s) not in 'sample_genotype_data', check {paste(extra_samples, collapse=\", \")}")
    )
    
    ## Check that there are no cases where either
    ## 1. Two or more anchor samples with the different Subject_ID(s) have same Genotype_Group_ID
    ## 2. Two or more anchor samples with the same Subject_ID have different Genotype_Group_ID(s)
    anchor_samples_consistency <- data.frame(Sample_ID = anchor_samples) %>% 
        left_join(
            sample_genotype_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")],
            by="Sample_ID"
        ) %>% 
        filter(!is.na(Genotype_Group_ID)) %>% 
        group_by(Genotype_Group_ID) %>% 
        mutate(
            n_Subject_ID = length(unique(Subject_ID))
        ) %>% 
        ungroup() %>% 
        group_by(Subject_ID) %>% 
        mutate(
            n_Genotype_Group_ID = length(unique(Genotype_Group_ID))
        ) %>% 
        ungroup()
    subject_inconsistent_samples <- anchor_samples_consistency %>% 
        filter(n_Subject_ID != 1) %>% 
        arrange(Sample_ID) %>% 
        pull(Sample_ID)
    genotype_inconsistent_samples <- anchor_samples_consistency %>% 
        filter(n_Genotype_Group_ID != 1) %>% 
        arrange(Sample_ID) %>% 
        pull(Sample_ID)
    assert_that(
        length(subject_inconsistent_samples) == 0,
        msg = glue("'anchor_samples' contains Sample_ID(s) in the same Genotype_Group_ID
                   but different Subject_ID(s), check {paste(subject_inconsistent_samples, collapse=\", \")}")
    )
    assert_that(
        length(genotype_inconsistent_samples) == 0,
        msg = glue("'anchor_samples' contains Sample_ID(s) with the same Subject_ID
                   but in different Genotype_Group_ID(s), check {paste(genotype_inconsistent_samples, collapse=\", \")}")
    )
    
    ## TODO: output a warning if an anchor sample contradicts a majority sample
}

################################################################################
##################               GRAPH HELPERS                ##################
################################################################################

# .plot_graph <- function(graph, to_write=FALSE) {
#     edge_colors <- if (!is.null(E(graph)$edge_colors)) E(graph)$edge_colors else "grey"
#     color <- if (!is.null(V(graph)$color)) V(graph)$color else "orange"
#     layout_custom <- with_seed(layout_nicely(graph), seed=1987)
#     vertex.label.cex <- 0.6
#     vertex.size <- 5
#     edge.arrow.size <- 0.8
#     edge.width <- 3
#     if (to_write) {
#         vertex.label.cex <- 2
#         vertex.size <- 5
#         edge.arrow.size <- 3
#         edge.width <- 8
#     }
#     
#     my_plot <- plot(graph, vertex.size=vertex.size, vertex.label.cex=vertex.label.cex, edge.arrow.size=edge.arrow.size, 
#                     edge.width=edge.width, edge.color=edge_colors, vertex.color=color, vertex.label.color="black", 
#                     vertex.frame.color="transparent", layout=layout_custom)
#     return(my_plot)
# }

.generate_graph <- function(
        relabel_data, 
        graph_type=c("label", "genotype", "combined"), 
        ghost_relabel_data=NULL, 
        anchor_samples=character(0),
        populate_plotting_attributes=FALSE
) {
    graph_type_mapping <- list(
        label = "Subject_ID",
        genotype = "Genotype_Group_ID",
        combined = NA_character_
    )
    
    graph_type <- as.character(graph_type)
    graph_type <- match.arg(graph_type)
    
    ghost_samples <- character(0)
    if (!is.null(ghost_relabel_data)) {
        ghost_samples <- ghost_relabel_data[, "Sample_ID", drop=TRUE]
    }
    
    if (graph_type == "combined") {
        genotype_graph <- .generate_graph(relabel_data, "genotype")
        E(genotype_graph)$genotypes <- TRUE
        label_graph <- .generate_graph(relabel_data, "label", ghost_relabel_data)
        E(label_graph)$labels <- TRUE 
        graph <- graph.union(genotype_graph, label_graph, byname=TRUE)
        E(graph)[is.na(E(graph)$genotypes)]$genotypes <- FALSE
        E(graph)[is.na(E(graph)$labels)]$labels <- FALSE
        E(graph)$concordant <- E(graph)$genotypes & E(graph)$labels
    } else {
        if(graph_type == "label" & !is.null(ghost_relabel_data)) {
            common_cols <- intersect(colnames(relabel_data), colnames(ghost_relabel_data))
            relabel_data <- rbind(relabel_data[, common_cols], ghost_relabel_data[, common_cols])
        }
        group_col <- graph_type_mapping[[graph_type]]
        edges <- relabel_data %>%
            group_by_at(group_col) %>%
            mutate(
                sample_a = Sample_ID,
                sample_b = list(Sample_ID)
            ) %>%
            ungroup() %>% 
            unnest(sample_b) %>%
            transmute(
                sample1 = pmin(sample_a, sample_b),
                sample2 = pmax(sample_a, sample_b)
            ) %>%
            filter(sample1 != sample2) %>%
            distinct()
        vertices <- relabel_data[, "Sample_ID", drop=FALSE]
        graph <- graph_from_data_frame(edges, vertices=vertices, directed=FALSE)
    }
    
    if (!populate_plotting_attributes) {return(graph)}
    ## Specify vertex and edge attributes for plotting
    sample_shapes <- data.frame(Sample_ID=names(V(graph))) %>% 
        left_join(unsolved_all_data, by="Sample_ID") %>% 
        pull(SwapCat_Shape)
    V(graph)$shape <- sample_shapes
    
    anchor_samples <- intersect(anchor_samples, V(graph)$name)
    label_not_found_samples <- V(graph)$name[grepl(LABEL_NOT_FOUND, V(graph)$name)]
    V(graph)$color <- "orange"
    V(graph)[anchor_samples]$color <- "forestgreen"
    V(graph)[label_not_found_samples]$color <- "firebrick"
    V(graph)$size <- 12
    if (graph_type == "combined") {
        V(graph)[ghost_samples]$color <- "lightgrey"
        E(graph)$color <- ifelse(E(graph)$concordant, "forestgreen", ifelse(E(graph)$genotypes, "orange", "cornflowerblue"))
        E(graph)[.from(ghost_samples)]$color <- "lightgrey"
    } else if (graph_type == "label") {
        V(graph)[ghost_samples]$color <- "lightgrey"
        E(graph)$color <- "cornflowerblue"
        E(graph)[.from(ghost_samples)]$color <- "lightgrey"
    } else {
        E(graph)$color <- "orange"
    }
    E(graph)$width <- 6
    
    return(graph)
}

 .generate_corrections_graph <- function(relabel_data) {
    applied_relabels <- relabel_data %>% 
        filter(Init_Sample_ID != Sample_ID) %>% 
        select(
            Init_Sample_ID,
            Sample_ID
        )
    corrections_graph <- graph_from_data_frame(applied_relabels, directed=TRUE)
    ghost_samples <- intersect(V(corrections_graph)$name, relabel_data %>% filter(is.na(Genotype_Group_ID)) %>% pull(Sample_ID))
    V(corrections_graph)$color <- "orange"
    V(corrections_graph)[ghost_samples]$color <- "lightgrey"
    V(corrections_graph)$size <- 24
    E(corrections_graph)$width <- 9
    return(corrections_graph)
}

## TODO: drop duplicate cycles
.find_directed_cycles <- function(graph, cutoff=1) {
    assert_that(is_directed(graph), 
                msg="param 'graph' must be directed")
    cycles <- list()
    for (vertex in V(graph)) {
        in_neighbors <- names(neighbors(graph, vertex, mode="in"))
        for (in_neighbor in in_neighbors) {
            simple_paths <- all_simple_paths(graph, vertex, in_neighbor, mode="out", cutoff=cutoff)
            cycles <- append(cycles, simple_paths)
        }
    }
    cycles <- lapply(cycles, names)
    return(cycles)
}

.find_all_relabel_cycles <- function(relabels_graph) {
    ## TODO: add warning if there are multiple deficiencies
    relabels <- data.frame(relabel_from=character(0), relabel_to=character(0))
    all_relabeled_samples <- NULL
    all_cycles <- list()
    cutoff <- 1
    while (vcount(relabels_graph) > 0 && cutoff < max(table(components(relabels_graph)$membership))) {
        curr_cycles <- .find_directed_cycles(relabels_graph, cutoff=cutoff)
        for (curr_cycle in curr_cycles) {
            if (!any(curr_cycle %in% all_relabeled_samples)) {
                all_cycles <- append(all_cycles, list(curr_cycle))
                all_relabeled_samples <- c(all_relabeled_samples, curr_cycle)
                relabels_graph <- delete_vertices(relabels_graph, curr_cycle)
            }
        }
        cutoff <- cutoff + 1
    }
    
    ## From the found cycles, construct relabels dataframe
    for (curr_cycle in all_cycles) {
        n <- length(curr_cycle)
        curr_relabels <- data.frame(
            relabel_from = curr_cycle,
            relabel_to = c(curr_cycle[2:n], curr_cycle[1])
        )
        relabels <- rbind(relabels, curr_relabels)
    }
    
    return(relabels)
}

.swap_cats_to_graph <- function(swap_cats) {
    n <- nrow(swap_cats)
    if (ncol(swap_cats) < 2) {
        swap_cats_adjacency_mat <- matrix(TRUE, nrow=n, ncol=n)
    } else {
        swap_cats_adjacency_mat <- matrix(FALSE, nrow=n, ncol=n)
        swap_cats_int <- interaction(swap_cats[, 2:ncol(swap_cats)], sep=":")
        for (level in levels(swap_cats_int)) {
            idx <- swap_cats_int == level
            swap_cats_adjacency_mat[idx, idx] <- TRUE
        }
    }
    rownames(swap_cats_adjacency_mat) <- colnames(swap_cats_adjacency_mat) <- swap_cats$Sample_ID
    swap_cats_graph <- graph_from_adjacency_matrix(swap_cats_adjacency_mat, mode="directed")
    return(swap_cats_graph)
} 

################################################################################
##################                  TESTING                   ##################
################################################################################

load_test_case <- function(test_name) {
    library(openxlsx)
    setwd("/Users/charlesdeng/Workspace/mislabeling")
    sample_genotype_data <- read.xlsx("simulations/tests/sample_genotype_data_tests.xlsx", sheet=test_name)
    swap_cats <- NULL
    try({
        swap_cats <- read.xlsx("simulations/tests/swap_cats_tests.xlsx", sheet=test_name)
        for (i in 2:ncol(swap_cats)) {
            swap_cats[[i]] <- as.factor(swap_cats[[i]])
        }
    }, silent=TRUE)
    my_mislabel_solver <- new("MislabelSolver", sample_genotype_data, swap_cats)
    return(my_mislabel_solver)
}

#my_mislabel_solver <- new("MislabelSolver", sample_genotype_data, swap_cats, anchor_samples); object <- my_mislabel_solver
#my_mislabel_solver <- load_test_case("Example"); object <- my_mislabel_solver
#my_mislabel_solver <- load_test_case("1-3C-2G-unswapped"); object <- my_mislabel_solver
# my_mislabel_solver <- load_test_case("comprehensive_search_puzzle"); object <- my_mislabel_solver
# my_mislabel_solver <- load_test_case("1-4C-1-3C-1-2C-2G"); object <- my_mislabel_solver
#my_mislabel_solver <- load_test_case("1-4C-2G-incycle"); object <- my_mislabel_solver
# my_mislabel_solver <- load_test_case("1-4C-1G"); object <- my_mislabel_solver
# my_mislabel_solver <- load_test_case("example"); object <- my_mislabel_solver
# my_mislabel_solver <- load_test_case("2G"); object <- my_mislabel_solver


# setMethod("solve_local_search_bulk", "MislabelSolver",
#           function(object, objective=c("genotype_entropy", "hamming_distance"), n_iter=1, frac_per_iter=0.20, include_ghost=FALSE) {
#               if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
#                   return(object)
#               }
# 
#               search_objectives <- list(
#                   genotype_entropy = .objective_genotype_entropy,
#                   hamming_distance = .objective_hamming_distance
#               )
#               objective_options <- names(search_objectives)
#               objective <- as.character(objective)
#               objective <- match.arg(objective, names(search_objectives))
# 
#               ## Bulk local search only works for genotype_entropy
#               assert_that(objective == "genotype_entropy")
# 
#               if (objective %in% objective_options) {
#                   objective_function <- search_objectives[[objective]]
#               } else {
#                   stop(glue("unrecognized value for 'objective': {objective}"))
#               }
# 
#               for (i in 1:n_iter) {
#                   print(glue("Iteration {i}:: local search with 'objective' \"{objective}\""))
#                   neighbors <- .find_neighbors(object, include_ghost) %>%
#                       left_join(
#                           rbind(object@.solve_state$unsolved_relabel_data[, c("Sample_ID", "Component_ID")],
#                                 object@.solve_state$unsolved_ghost_data[, c("Sample_ID", "Component_ID")]),
#                           by=c("Sample_A"="Sample_ID")
#                       ) %>%
#                       dplyr::rename(relabel_from=Sample_A, relabel_to=Sample_B)
#                   print(glue("{nrow(neighbors)} candidate swaps being evaluated..."))
#                   all_component_ids <- sort(unique(object@.solve_state$unsolved_relabel_data$Component_ID))
#                   relabels_list <- list()
#                   curr_idx <- 1
#                   for (curr_component_id in all_component_ids) {
#                       cc_relabel_data <- rbind(object@.solve_state$unsolved_relabel_data %>% filter(Component_ID == curr_component_id),
#                                                object@.solve_state$unsolved_ghost_data %>% filter(Component_ID == curr_component_id))
#                       cc_neighbors <- neighbors %>% filter(Component_ID == curr_component_id)
#                       calc_swapped_objective <- function(swap_from, swap_to) {
#                           swap_from_index <- which(cc_relabel_data$Sample_ID == swap_from)[[1]]
#                           swap_to_index <- which(cc_relabel_data$Sample_ID == swap_to)[[1]]
#                           swap_from_subject <- cc_relabel_data[swap_from_index, "Subject_ID"][[1]]
#                           swap_to_subject <- cc_relabel_data[swap_to_index, "Subject_ID"][[1]]
#                           cc_relabel_data[swap_from_index, "Sample_ID"] <- swap_to
#                           cc_relabel_data[swap_from_index, "Subject_ID"] <- swap_to_subject
#                           cc_relabel_data[swap_to_index, "Sample_ID"] <- swap_from
#                           cc_relabel_data[swap_to_index, "Subject_ID"] <- swap_from_subject
#                           new_objective <- objective_function(cc_relabel_data[!is.na(cc_relabel_data$Genotype_Group_ID), ])
#                           return(new_objective)
#                       }
# 
#                       if (nrow(cc_neighbors) == 0) {next}
#                       cc_neighbor_objectives <- cc_neighbors %>%
#                           mutate(
#                               base = objective_function(cc_relabel_data),
#                               objective = mapply(calc_swapped_objective, swap_from=relabel_from, swap_to=relabel_to),
#                               delta = objective - base
#                           )
#                       cc_candidate_relabels <- cc_neighbor_objectives %>%
#                           filter(delta < 0) %>%
#                           left_join(cc_relabel_data[, c("Sample_ID", "Genotype_Group_ID")], by=c("relabel_from"="Sample_ID")) %>%
#                           left_join(cc_relabel_data[, c("Sample_ID", "Genotype_Group_ID")], by=c("relabel_to"="Sample_ID")) %>%
#                           arrange(delta)
# 
#                       if (nrow(cc_candidate_relabels) == 0) {next}
# 
#                       ## Find minimum number of mislabels for the component
#                       n_min_mislabels <- sum(cc_relabel_data %>%
#                                                  group_by(Genotype_Group_ID) %>%
#                                                  summarize(
#                                                      count = n(),
#                                                      max_count = max(table(Subject_ID)),
#                                                      min_mislabels = count - max_count
#                                                  ) %>% pull(min_mislabels))
#                       ## Determine number of swaps to look for
#                       n_swaps_limit <- max(1, floor(n_min_mislabels * frac_per_iter * 0.5))
# 
#                       cc_relabels <- data.frame(matrix(data=NA, nrow=n_swaps_limit, ncol=2, dimnames=list(c(), c("relabel_from", "relabel_to"))))
#                       for (j in 1:n_swaps_limit) {
#                           top_relabel <- cc_candidate_relabels %>%
#                               head(1)
#                           cc_relabels[j, ] <- top_relabel[, c("relabel_from", "relabel_to")]
#                           from_genotype_id <- "Genotype_Group_ID.x"
#                           to_genotype_id <- "Genotype_Group_ID.x"
#                           cc_candidate_relabels <- cc_candidate_relabels %>%
#                               filter(!(Genotype_Group_ID.x %in% c(from_genotype_id, to_genotype_id)),
#                                      !(Genotype_Group_ID.y %in% c(from_genotype_id, to_genotype_id)))
#                           if (nrow(cc_candidate_relabels) == 0) {
#                               break
#                           }
#                       }
# 
#                       cc_relabels
#                       relabels_list[[curr_idx]] <- cc_relabels
#                       ##relabels[curr_idx, c("relabel_from", "relabel_to")] <- cc_relabels
#                       curr_idx <- curr_idx + 1
#                   }
# 
#                   relabels <- do.call(rbind, relabels_list)
#                   relabels <- rbind(relabels, data.frame(relabel_from=relabels$relabel_to, relabel_to=relabels$relabel_from))
#                   object <- .relabel_samples(object, relabels)
#               }
# 
#               return(object)
#           }
# )

# .genotype_group_vote <- function(object, unsolved=TRUE) {
#     if (unsolved) {
#         relabel_data <- object@.solve_state$unsolved_relabel_data
#     } else {
#         relabel_data <- object@.solve_state$relabel_data
#     }
#     votes <- relabel_data %>% 
#         group_by(Subject_ID, Genotype_Group_ID) %>% 
#         summarize(n=n(), .groups='drop') %>% 
#         pivot_wider(names_from=Subject_ID, values_from=n) %>% 
#         mutate_all(~ifelse(is.na(.), 0, .)) %>% 
#         tibble::column_to_rownames("Genotype_Group_ID") %>% 
#         as.matrix()
#     return(votes)
# }
# setMethod("solve_local_search", "MislabelSolver", 
#           function(object, objective=c("genotype_entropy", "hamming_distance"), 
#                    n_iter=1, 
#                    include_ghost=FALSE, 
#                    filter_concordant_vertices=FALSE) {
#               print("Starting local search")
#               if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
#                   return(object)
#               }
#               
#               search_objectives <- list(
#                   genotype_entropy = .objective_genotype_entropy,
#                   hamming_distance = .objective_hamming_distance
#               )
#               objective_options <- names(search_objectives)
#               objective <- as.character(objective)
#               objective <- match.arg(objective, names(search_objectives))
#               if (objective %in% objective_options) {
#                   objective_function <- search_objectives[[objective]]
#               } else {
#                   stop(glue("unrecognized value for 'objective': {objective}"))
#               }
#               
#               for (i in 1:n_iter) {
#                   print(glue("Iteration {i}:: local search with 'objective' \"{objective}\", 'include_ghost'={include_ghost},
#                              'filter_concordant_vertices'={filter_concordant_vertices}"))
#                   neighbors <- .find_neighbors(object, include_ghost, filter_concordant_vertices) %>% 
#                       left_join(
#                           rbind(object@.solve_state$unsolved_relabel_data[, c("Sample_ID", "Component_ID")],
#                                 object@.solve_state$unsolved_ghost_data[, c("Sample_ID", "Component_ID")]), 
#                           by=c("Sample_A"="Sample_ID")
#                       )
#                   print(glue("{nrow(neighbors)} candidate swaps being evaluated..."))
#                   all_component_ids <- sort(unique(object@.solve_state$unsolved_relabel_data$Component_ID))
#                   relabels <- data.frame(matrix(data=NA, nrow=length(all_component_ids), ncol=2, dimnames=list(c(), c("relabel_from", "relabel_to"))))
#                   curr_idx <- 1
#                   for (curr_component_id in all_component_ids) {
#                       cc_relabel_data <- rbind(object@.solve_state$unsolved_relabel_data %>% filter(Component_ID == curr_component_id),
#                                                object@.solve_state$unsolved_ghost_data %>% filter(Component_ID == curr_component_id))
#                       cc_neighbors <- neighbors %>% filter(Component_ID == curr_component_id)
#                       calc_swapped_objective <- function(swap_from, swap_to) {
#                           swap_from_index <- which(cc_relabel_data$Sample_ID == swap_from)[[1]]
#                           swap_to_index <- which(cc_relabel_data$Sample_ID == swap_to)[[1]]
#                           swap_from_subject <- cc_relabel_data[swap_from_index, "Subject_ID"][[1]]
#                           swap_to_subject <- cc_relabel_data[swap_to_index, "Subject_ID"][[1]]
#                           cc_relabel_data[swap_from_index, "Sample_ID"] <- swap_to
#                           cc_relabel_data[swap_from_index, "Subject_ID"] <- swap_to_subject
#                           cc_relabel_data[swap_to_index, "Sample_ID"] <- swap_from
#                           cc_relabel_data[swap_to_index, "Subject_ID"] <- swap_from_subject
#                           new_objective <- objective_function(cc_relabel_data[!is.na(cc_relabel_data$Genotype_Group_ID), ])
#                           return(new_objective)
#                       }
#                       
#                       if (nrow(cc_neighbors) == 0) {next}
#                       cc_neighbor_objectives <- cc_neighbors %>% 
#                           mutate(
#                               base = objective_function(cc_relabel_data),
#                               objective = mapply(calc_swapped_objective, swap_from=Sample_A, swap_to=Sample_B),
#                               delta = objective - base
#                           )
#                       cc_relabels <- cc_neighbor_objectives %>%
#                           filter(delta < 0, delta == min(delta)) 
#                       if (nrow(cc_relabels) == 0) {next}
#                       cc_relabels <- cc_relabels %>% 
#                           sample_n(1) %>%
#                           transmute(
#                               relabel_from=Sample_A,
#                               relabel_to=Sample_B
#                           )
#                       relabels[curr_idx, c("relabel_from", "relabel_to")] <- cc_relabels
#                       
#                       curr_idx <- curr_idx + 1
#                   }
#                   
#                   relabels <- relabels[!is.na(relabels[, 1]) & !is.na(relabels[, 2]), ]
#                   relabels <- rbind(relabels, data.frame(relabel_from=relabels$relabel_to, relabel_to=relabels$relabel_from))
#                   object <- .relabel_samples(object, relabels)
#               }
#               
#               return(object)
#           }
# )

# .generate_corrections_graph <- function(relabel_data, swap_cats) {
#     applied_relabels <- relabel_data %>% 
#         filter(Init_Sample_ID != Sample_ID) %>% 
#         select(
#             Init_Sample_ID,
#             Sample_ID
#         ) %>% 
#         left_join(swap_cats, by="Sample_ID") %>% 
#             
#             mutate(
#             SwapCat_ID = ifelse(
#                 is.na(Subject_ID),
#                 sapply(Sample_ID, \(x) unlist(strsplit(x, split=":"))[2]),
#                 Subject_ID
#             )
#         )
#     corrections_graph <- graph_from_data_frame(applied_relabels, directed=TRUE)
#     ghost_samples <- intersect(V(corrections_graph)$name, relabel_data %>% filter(is.na(Genotype_Group_ID)) %>% pull(Sample_ID))
#     V(corrections_graph)$color <- "orange"
#     V(corrections_graph)[ghost_samples]$color <- "lightgrey"
#     V(corrections_graph)$size <- 12
#     E(corrections_graph)$width <- 6
#     E(corrections_graph)$color <- "forestgreen"
#     
#     return(corrections_graph)
# }


# .objective_genotype_entropy = function(relabel_data) {
#     p <- table(relabel_data$Subject_ID, relabel_data$Genotype_Group_ID)
#     genotype_counts <- colSums(p)
#     p <- scale(p, scale=genotype_counts, center=FALSE)
#     return(sum(genotype_counts * colSums(-(p * log(p)), na.rm=TRUE)))
# }
# 
# .objective_hamming_distance <- function(relabel_data) {
#     labels_graph <- .generate_graph(relabel_data, "label")
#     genotype_graph <- .generate_graph(relabel_data, "genotype")
#     int_graph <- graph.intersection(labels_graph, genotype_graph)
#     return(ecount(labels_graph) + ecount(genotype_graph) - 2*ecount(int_graph))
# }



