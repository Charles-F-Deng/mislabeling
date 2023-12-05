setClass("MislabelSolver",
         representation(
             sample_genotype_data = "data.frame",
             swap_cats = "data.frame",
             anchor_samples = "character",
             .solve_state = "list"
         ),
         prototype(
             swap_cats = NULL,
             anchor_samples = character(0),
             .solve_state = list()
         ),
         
)

setMethod("initialize", "MislabelSolver",
          function(.Object, sample_genotype_data, swap_cats=NULL, anchor_samples=character(0)) {
              .validate_sample_genotype_data(sample_genotype_data)
              
              if (is.null(swap_cats)) {
                  swap_cats <- sample_genotype_data[, "Sample_ID", drop=FALSE]
              }
              .validate_swap_cats(sample_genotype_data, swap_cats)
              .validate_anchor_samples(sample_genotype_data, anchor_samples)
              
              .Object@sample_genotype_data <- sample_genotype_data
              .Object@swap_cats <- swap_cats
              .Object@anchor_samples <- anchor_samples
              
              swap_cats_graph <- .swap_cats_to_graph(.Object@swap_cats)
              relabel_data <- .Object@sample_genotype_data %>% 
                  mutate(
                      Init_Sample_ID = Sample_ID,
                      Init_Subject_ID = Subject_ID,
                      Solved = FALSE
                  )
              unsolved_relabel_data <- relabel_data
              putative_subjects <- data.frame(Genotype_Group_ID = character(0),
                                              Subject_ID = character(0))
              ambiguous_subjects <- list()
              .Object@.solve_state <- list(
                  swap_cats_graph = swap_cats_graph,
                  relabel_data = relabel_data,
                  unsolved_relabel_data = unsolved_relabel_data,
                  putative_subjects = putative_subjects,
                  ambiguous_subjects = ambiguous_subjects
              )
              .Object <- .update_solve_state(.Object, initialize=TRUE)
              return(.Object)
          }
)

setMethod("plot", "MislabelSolver",
          function(x, y=NULL, filter_solved=TRUE, corrections=FALSE) {
              if (corrections) {
                  .plot_graph(.generate_corrections_graph(x@.solve_state$relabel_data))
                  
              } else {
                  if (filter_solved) {
                      relabel_data <- x@.solve_state$unsolved_relabel_data
                  } else {
                      relabel_data <- x@.solve_state$relabel_data
                  }
                  .plot_graph(.generate_graph(relabel_data, graph_type = "combined"))
              }
          }
)

setGeneric("solve_majority_search", function(object) {
    standardGeneric("solve_majority_search")
})
setGeneric("solve_comprehensive_search", function(object) {
    standardGeneric("solve_comprehensive_search")
})
setGeneric("solve_local_search", function(object, ...) {
    standardGeneric("solve_local_search")
})
setGeneric("solve", function(object, ...) {
    standardGeneric("solve")
})
setGeneric("write_corrections", function(object) {
    standardGeneric("write_corrections")
})

setMethod("solve_majority_search", "MislabelSolver",
          function(object) {
              if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                  return(object)
              }
              
              ## 1. Update putative subjects
              votes <- .genotype_group_vote(object, unsolved=TRUE)
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
              relabels <- .find_relabel_cycles_from_putative_subjects(object)
              
              ## 3. Relabel samples and update solve state
              object <- .relabel_samples(object, relabels)
              return(object)
          }
)

## TODO: fix the case that I showed Ryan, 4 samples and 2 groups but it's actually unambiguous
setMethod("solve_comprehensive_search", "MislabelSolver",
          function(object) {
              if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                  return(object)
              }
              
              ## 1. Update putative subjects
              MAX_GENOTYPES <- 9
              component_ids <- unique(object@.solve_state$unsolved_relabel_data$Component_ID)
              possible_neighbors <- .find_neighbors(object)
              for (component_id in component_ids) {
                  connected_component <- object@.solve_state$unsolved_relabel_data %>% filter(Component_ID == component_id)
                  
                  ## For now, pass out of components where number of Genotype_Group(s) doesn't match number of Subject_ID(s)
                  if (length(unique(connected_component$Genotype_Group_ID)) != length(unique(connected_component$Subject_ID))) {
                      next
                  }
                  
                  locked_genotypes <- object@.solve_state$putative_subjects %>% 
                      filter(Genotype_Group_ID %in% unique(connected_component$Genotype_Group_ID))
                  
                  ## Hotfix for the case that I showed Ryan, 4 samples and 2 groups but it's actually unambiguous
                  ## If it's a plurality within it's own genotype group, and it can't be swapped, lock it
                  votes <- .genotype_group_vote(object)
                  for (genotype_group in rownames(votes)) {
                      if (genotype_group %in% locked_genotypes$Genotype_Group_ID) {next}
                      max_votes <- max(votes[genotype_group, ])
                      plurality_subjects <- colnames(votes)[votes[genotype_group, ] == max_votes]
                      genotype_locked <- FALSE
                      for (subject_id in plurality_subjects) {
                          sample_ids <- connected_component %>%
                              filter(Subject_ID == subject_id,
                                     Genotype_Group_ID == genotype_group) %>%
                              pull(Sample_ID)
                          # If any of the sample ids can't be swapped, lock this genotype/subject
                          for (sample_id in sample_ids) {
                              if (!(sample_id %in% possible_neighbors$Sample_A | sample_id %in% possible_neighbors$Sample_B)) {
                                  locked_genotypes %<>%
                                      rbind(data.frame(
                                          Genotype_Group_ID = genotype_group,
                                          Subject_ID = subject_id))
                                  genotype_locked <- TRUE
                                  break
                              }
                          }
                          if (genotype_locked) {break}
                      }
                  }

                  free_genotypes <- setdiff(connected_component$Genotype_Group_ID, locked_genotypes$Genotype_Group_ID)
                  free_subjects <- setdiff(connected_component$Subject_ID, locked_genotypes$Subject_ID)
                  if (length(free_genotypes) > MAX_GENOTYPES) {
                      next
                  }
                  if (length(free_genotypes) > 0) {
                      perm_genotypes <- as.data.frame(t(simplify2array(permn(free_subjects))))
                      colnames(perm_genotypes) <- sort(free_genotypes)
                      for (genotype_id in locked_genotypes$Genotype_Group_ID) {
                          subject_id <- locked_genotypes[locked_genotypes$Genotype_Group_ID == genotype_id, "Subject_ID"][[1]]
                          perm_genotypes[[genotype_id]] <- subject_id
                      }
                  } else {
                      perm_genotypes <- as.data.frame(t(locked_genotypes$Subject_ID))
                      colnames(perm_genotypes) <- locked_genotypes$Genotype_Group_ID
                  }
                  perm_genotypes <- as.matrix(perm_genotypes)
                  
                  ## Only including cases where the number of samples in Genotype_Group_ID and Subject_ID align
                  # genotype_counts <- connected_component %>%
                  #     count(Genotype_Group_ID) %>%
                  #     column_to_rownames(var = "Genotype_Group_ID")
                  # genotype_counts <- genotype_counts[colnames(perm_genotypes), ]
                  # subject_counts <- connected_component %>%
                  #     count(Subject_ID) %>%
                  #     column_to_rownames(var = "Subject_ID")
                  # count_concordance <- sapply(1:nrow(perm_genotypes), \(x) all(subject_counts[perm_genotypes[x, ], ] == genotype_counts))
                  # perm_genotypes <- perm_genotypes[count_concordance, , drop=FALSE]
                  # if (nrow(perm_genotypes) == 0) {
                  #     next
                  # }
                  
                  ## For each Genotype_Group_ID/Subject_ID permutation, determine the implied Sample_ID/Subject_ID permutation
                  perm_samples <- matrix(nrow=nrow(perm_genotypes), ncol=nrow(connected_component), dimnames=list(NULL, connected_component$Sample_ID))
                  for (i in 1:nrow(perm_samples)) {
                      perm_samples[i, ] <- perm_genotypes[i, ][connected_component$Genotype_Group_ID]
                  }
                  
                  ## For each Sample_ID/Subject_ID permutation, determine number of implied mislabels
                  curr_subjects <- t(replicate(n=nrow(perm_samples), connected_component$Subject_ID))
                  colnames(curr_subjects) <- colnames(perm_samples)
                  min_mislabels <- connected_component %>% 
                      group_by(Genotype_Group_ID) %>% 
                      summarize(
                          n = n(),
                          n_max_count = max(table(Subject_ID)),
                          n_min_mislabels = n - n_max_count
                      )
                  n_min_mislabels <- sum(min_mislabels$n_min_mislabels)
                  perm_genotypes <- as.data.frame(perm_genotypes) %>% 
                      mutate(
                          n_min_mislabels = n_min_mislabels,
                          n_implied_mislabels = rowSums(perm_samples != curr_subjects)
                      ) %>% 
                      arrange(n_implied_mislabels)
                  
                  ## If there were ties, track equal alternatives in .solve_state$ambiguous_subjects
                  tied_solutions <- perm_genotypes %>% 
                      filter(n_implied_mislabels == min(perm_genotypes$n_min_mislabels)) %>% 
                      select(-n_min_mislabels, -n_implied_mislabels)
                  if (nrow(tied_solutions > 0)) {
                      for (genotype_id in colnames(tied_solutions)) {
                          unique_subject_ids <- unique(tied_solutions[[genotype_id]])
                          if (length(unique_subject_ids) > 1) {
                              object@.solve_state$ambiguous_subjects[[genotype_id]] <- unique_subject_ids
                          }
                      }
                  }
                  
                  ## To find a single solution, take top row (fewest mislabels)
                  new_putative_subjects <- perm_genotypes[1, ] %>% 
                      select(-n_min_mislabels, -n_implied_mislabels) %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      rownames_to_column()
                  colnames(new_putative_subjects) <- c("Genotype_Group_ID", "Subject_ID")
                  new_putative_subjects %<>% 
                      anti_join(object@.solve_state$putative_subjects, by=c("Genotype_Group_ID", "Subject_ID"))
                  object <- .update_putative_subjects(object, new_putative_subjects)
                  
                  ## Find relabel cycles
                  relabels <- .find_relabel_cycles_from_putative_subjects(object)
                  
                  ## Relabel samples and update solve state
                  object <- .relabel_samples(object, relabels) 
              }
              
              return(object)
          }
)

setMethod("solve_local_search", "MislabelSolver",
          function(object, objective=c("genotype_entropy", "hamming_distance"), n_iter=1) {
              search_objectives <- list(
                  genotype_entropy = .objective_genotype_entropy,
                  hamming_distance = .objective_hamming_distance
              )
              objective_options <- names(search_objectives)
              objective <- as.character(objective)
              objective <- match.arg(objective, names(search_objectives))
              if (objective %in% objective_options) {
                  objective_function <- search_objectives[[objective]]
              } else {
                  stop(glue("unrecognized value for 'objective': {objective}"))
              }
              
              calc_swapped_objective <- function(swap_from, swap_to) {
                  relabel_data <- object@.solve_state$unsolved_relabel_data
                  swap_from_index <- which(relabel_data$Sample_ID == swap_from)[[1]]
                  swap_to_index <- which(relabel_data$Sample_ID == swap_to)[[1]]
                  swap_from_subject <- relabel_data[swap_from_index, "Subject_ID"][[1]]
                  swap_to_subject <- relabel_data[swap_to_index, "Subject_ID"][[1]]
                  relabel_data[swap_from_index, c("Sample_ID", "Subject_ID")] <- c(swap_to, swap_to_subject)
                  relabel_data[swap_to_index, c("Sample_ID", "Subject_ID")] <- c(swap_from, swap_from_subject)
                  return(objective_function(relabel_data))
              }
              
              for (curr_iter in 1:n_iter) {
                  if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                      return(object)
                  }
                  
                  neighbors <- .find_neighbors(object) %>% 
                      left_join(
                          object@.solve_state$unsolved_relabel_data[, c("Sample_ID", "Component_ID")], 
                          by=c("Sample_A"="Sample_ID")
                      )
                  if (nrow(neighbors) == 0) { break }
                  neighbor_objectives <- neighbors %>% 
                      mutate(
                          base = objective_function(object@.solve_state$unsolved_relabel_data),
                          objective = mapply(calc_swapped_objective, swap_from=Sample_A, swap_to=Sample_B),
                          delta = objective - base
                      )
                  relabels <- neighbor_objectives %>%
                      group_by(Component_ID) %>% 
                      filter(delta < 0, delta == min(delta))
                  if (nrow(relabels) == 0) { break }
                  relabels <- relabels %>%
                      sample_n(1) %>%
                      ungroup(Component_ID) %>% 
                      transmute(
                          relabel_from=Sample_A,
                          relabel_to=Sample_B
                      )
                  relabels <- rbind(relabels, data.frame(relabel_from=relabels$relabel_to, relabel_to=relabels$relabel_from))
                  object <- .relabel_samples(object, relabels)
              }
              
              return(object)
          }
)

setMethod("solve", "MislabelSolver",
          function(object) {
              while (TRUE) {
                  if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
                      break
                  }
                  prev_solve_state <- object@.solve_state$unsolved_relabel_data
                  
                  object <- solve_comprehensive_search(object)
                  object <- solve_majority_search(object)
                  object <- solve_comprehensive_search(object)
                  object <- solve_local_search(object, n_iter=1)
                  
                  ## If current solve state is the same as previous, break
                  if (nrow(prev_solve_state) == nrow(object@.solve_state$unsolved_relabel_data)) {
                      if (all.equal(prev_solve_state, object@.solve_state$unsolved_relabel_data)) {
                          break
                      }
                  }
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
                  transmute(
                      Component_ID = Init_Component_ID,
                      Genotype_Group_ID,
                      Init_Subject_ID,
                      Init_Sample_ID,
                      Final_Subject_ID = Subject_ID,
                      Final_Sample_ID = Sample_ID,
                      n_Genotype_Group,
                      Init_Agreement = paste0(n_agree-1, " out of ", n_Genotype_Group - 1),
                      Status = case_when(
                          is.na(Final_Sample_ID) ~ "removed",
                          Init_Sample_ID != Final_Sample_ID ~ "corrected",
                          n_Genotype_Group == 1 ~ "ignored_single-sample",
                          n_agree < 2 ~ "ignored",
                          TRUE ~ "validated"
                      )
                  ) %>% 
                  arrange(Component_ID, Genotype_Group_ID, Final_Subject_ID, Final_Sample_ID)
              corrections_graph <- generate_corrections_graph(object@.solve_state$relabel_data)
              cluster_info <- components(corrections_graph)
              cluster_sizes <- cluster_info$csize[cluster_info$membership]
              names(cluster_sizes) <- names(cluster_info$membership)
              sample_summary$Mislabel_Type <- NA_character_
              sample_summary[sample_summary$Init_Sample_ID %in% names(cluster_sizes), "Mislabel_Type"] = as.character(cluster_sizes)
              sample_summary %<>% mutate(
                  Mislabel_Type = ifelse(
                      is.na(Mislabel_Type),
                      Mislabel_Type,
                      ifelse(
                          Mislabel_Type == "2", 
                          "Swap",
                          paste0(Mislabel_Type, "-cycle")
                      ) 
                  )
              )
              ## TODO: fix this
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
              
              init_solved_components <- object@.solve_state$relabel_data %>% 
                  group_by(Init_Component_ID) %>% 
                  summarize(
                      n_Subject = length(unique(Init_Subject_ID)),
                      n_Genotype_Group = length(unique(Genotype_Group_ID))
                  ) %>% 
                  filter(n_Subject == 1 & n_Genotype_Group == 1) %>% 
                  pull(Init_Component_ID)
              unsolved_components <- sort(setdiff(object@.solve_state$relabel_data$Init_Component_ID, init_solved_components))
              
              for (component in unsolved_components) {
                  save_file <- file.path(dir_name, glue(component, ".png"))
                  png(filename=save_file, width=1600, height=1600)
                  layout_matrix <- matrix(1:4, ncol=2)
                  layout(layout_matrix)
                  .plot_graph(.generate_graph(object@.solve_state$relabel_data %>% 
                                                  filter(Init_Component_ID == component) %>% 
                                                  mutate(Sample_ID = Init_Sample_ID, Subject_ID = Init_Subject_ID), graph_type="combined"), to_write=TRUE)
                  title(main=glue(component, ": Before corrections"), cex.main=2)
                  .plot_graph(.generate_corrections_graph(object@.solve_state$relabel_data %>% filter(Init_Component_ID == component)), to_write=TRUE)
                  title(main=glue(component, ": Applied corrections"), cex.main=2)
                  .plot_graph(.generate_graph(object@.solve_state$relabel_data %>% filter(Init_Component_ID == component), graph_type="combined"), to_write=TRUE)
                  title(main=glue(component, ": After corrections"), cex.main=2)
                  dev.off()
              }
          }
)

################################################################################
##################               SOLVE HELPERS                ##################
################################################################################

.update_solve_state <- function(object, initialize=FALSE) {
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        return(object)
    }
    
    ## 1. Assign a Component_ID for each Sample_ID in object@.solve_state$unsolved_relabel_data
    combined_graph <- .generate_graph(object@.solve_state$unsolved_relabel_data, graph_type="combined")
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data %>% 
        mutate(
            Component_ID = components(combined_graph)$membership[Sample_ID]
        )
    
    ## 2. Determine which Sample_ID(s) are newly solved
    component_data <- unsolved_relabel_data %>% 
        group_by(Component_ID) %>% 
        summarize(
            n_Genotype_Group_ID = length(unique(Genotype_Group_ID)),
            n_Subject_ID = length(unique(Subject_ID)),
            n_Sample_ID = length(Sample_ID)
        ) %>%
        mutate(
            Solved = ifelse(
                n_Genotype_Group_ID == 1 & n_Subject_ID == 1,
                TRUE,
                FALSE
            )
        )
    solved_components <- component_data %>% filter(Solved) %>% pull(Component_ID)
    unsolved_relabel_data %<>% mutate(
        Solved = ifelse(
            Component_ID %in% solved_components,
            TRUE,
            FALSE
        )
    )
    
    ## 3. Re-rank Component_ID(s) so that Component01 is the largest unsolved component
    n_components <- nrow(component_data)
    component_data %<>% arrange(Solved, desc(n_Sample_ID)) %>% 
        mutate(
            new_Component_ID = if (n_components != 0) 1:n_components else character(0),
            new_Component_ID = paste0("Component", formatC(new_Component_ID, width=nchar(n_components), format="d", flag="0"))
        )
    unsolved_relabel_data %<>% 
        left_join(component_data[, c("Component_ID", "new_Component_ID")], by="Component_ID") %>% 
        mutate(
            Component_ID = new_Component_ID
        ) %>% 
        select(-new_Component_ID)
    
    ## 4. Update putative_subjects
    new_putative_subjects <- unsolved_relabel_data %>% 
        filter(Solved) %>% 
        select(Genotype_Group_ID, Subject_ID) %>% 
        distinct() %>% 
        anti_join(object@.solve_state$putative_subjects, by=c("Genotype_Group_ID", "Subject_ID"))
    putative_subjects <- rbind(object@.solve_state$putative_subjects, new_putative_subjects)
    .validate_putative_subjects(object@sample_genotype_data, putative_subjects)
    
    ## 5. Update relabel_data and unsolved_relabel_data
    if (initialize) {
        unsolved_relabel_data %<>%
            select(Init_Sample_ID, Init_Subject_ID, Genotype_Group_ID, Component_ID, Sample_ID, Subject_ID, Solved)
        relabel_data <- unsolved_relabel_data %>% rename(
            Init_Component_ID = Component_ID
        )
    } else {
        relabel_data <- object@.solve_state$relabel_data %>% 
            left_join(unsolved_relabel_data[, c("Sample_ID", "Subject_ID", "Init_Sample_ID", "Solved")], 
                      by="Init_Sample_ID", suffix = c(".x", ".y")) %>% 
            mutate(
                Sample_ID = ifelse(!is.na(Sample_ID.y), Sample_ID.y, Sample_ID.x),
                Subject_ID = ifelse(!is.na(Subject_ID.y), Subject_ID.y, Subject_ID.x),
                Solved = ifelse(!is.na(Solved.y), Solved.y, Solved.x),
            ) %>% 
            select(-ends_with(".x"), -ends_with("y")) %>% 
            select(Init_Sample_ID, Init_Subject_ID, Genotype_Group_ID, Init_Component_ID, Sample_ID, Subject_ID, Solved)
    }
    unsolved_relabel_data <- unsolved_relabel_data %>% filter(!Solved)
    
    ## 6. Overwrite .solve_state
    object@.solve_state$relabel_data <- relabel_data
    object@.solve_state$unsolved_relabel_data <- unsolved_relabel_data
    object@.solve_state$putative_subjects <- putative_subjects
    
    return(object)
}

.update_putative_subjects <- function(object, proposed_putative_subjects) {
    if (nrow(proposed_putative_subjects) == 0) return(object)
    ## Only add Genotype_Group_ID/Subject_ID combinations if neither the 
    ## Genotype_Group_ID nor the Subject_ID are already in putative_subjects
    existing_genotypes <- object@.solve_state$putative_subjects$Genotype_Group_ID
    existing_subjects <- object@.solve_state$putative_subjects$Subject_ID
    proposed_putative_subjects %<>%
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
    .validate_relabels(object@.solve_state$unsolved_relabel_data, relabels)
    relabels %<>% 
        rename("Sample_ID" = "relabel_to") %>% 
        left_join(
            object@sample_genotype_data[, c("Sample_ID", "Subject_ID")],
            by="Sample_ID"
        )
    ## Call it relabeled sample ID instead
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data %>% 
        left_join(
            relabels, by=c("Sample_ID"="relabel_from"), suffix=c(".x", ".y")
        ) %>% 
        mutate(
            Sample_ID = ifelse(!is.na(Sample_ID.y), Sample_ID.y, Sample_ID),
            Subject_ID = ifelse(!is.na(Subject_ID.y), Subject_ID.y, Subject_ID.x)
        ) %>% 
        select(-ends_with(".x"), -ends_with(".y"))
    
    object@.solve_state$unsolved_relabel_data <- unsolved_relabel_data
    object <- .update_solve_state(object)
    return(object)
}

.find_relabel_cycles_from_putative_subjects <- function(object, putative_subjects=NULL) {
    if (is.null(putative_subjects)) {
        putative_subjects <- object@.solve_state$putative_subjects
    } else {
        .validate_putative_subjects(object@sample_genotype_data, putative_subjects)
    }
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    
    putative_subjects_list <- as.list(setNames(putative_subjects$Subject_ID, putative_subjects$Genotype_Group_ID))
    
    potential_mislabel_univ <- unsolved_relabel_data %>%
        mutate(
            Putative_Subject_ID = sapply(Genotype_Group_ID, \(x) putative_subjects_list[[x]]),
            Inferred_Correctly_Labeled = Subject_ID == Putative_Subject_ID
        ) %>% 
        filter(!Inferred_Correctly_Labeled) %>% 
        pull(Sample_ID)
    n <- length(potential_mislabel_univ)
    
    if (n == 0) {
        return(data.frame(relabel_from=character(0), relabel_to=character(0)))
    }
    
    relabels_adjacency_mat <- matrix(FALSE, nrow=n, ncol=n)
    rownames(relabels_adjacency_mat) <- colnames(relabels_adjacency_mat) <- potential_mislabel_univ
    for (sample_id in potential_mislabel_univ) {
        genotype_group_id <- unsolved_relabel_data[unsolved_relabel_data$Sample_ID == sample_id, ]$Genotype_Group_ID[[1]]
        inferred_subject_id <- ifelse(genotype_group_id %in% names(putative_subjects_list), putative_subjects_list[[genotype_group_id]], NA_character_)
        if (!is.na(inferred_subject_id)) {
            potential_relabels <- unsolved_relabel_data %>% 
                filter(
                    Sample_ID %in% potential_mislabel_univ,
                    Subject_ID == inferred_subject_id,
                    Genotype_Group_ID != genotype_group_id
                ) %>% 
                pull(Sample_ID)
            relabels_adjacency_mat[sample_id, potential_relabels] <- TRUE
        }
    }
    relabels_graph <- graph_from_adjacency_matrix(relabels_adjacency_mat, mode="directed")
    relabels_graph <- graph.intersection(relabels_graph, object@.solve_state$swap_cats_graph)
    
    ## Filter out samples that have no edges
    relabels_graph <- subgraph(relabels_graph, V(relabels_graph)[degree(relabels_graph, mode="all") > 0])
    
    relabels <- .find_all_relabel_cycles(relabels_graph)
    return(relabels)
}

.genotype_group_vote <- function(object, unsolved=TRUE) {
    if (unsolved) {
        relabel_data <- object@.solve_state$unsolved_relabel_data
    } else {
        relabel_data <- object@.solve_state$relabel_data
    }
    votes <- relabel_data %>% 
        group_by(Subject_ID, Genotype_Group_ID) %>% 
        summarize(n=n(), .groups='drop') %>% 
        pivot_wider(names_from=Subject_ID, values_from=n) %>% 
        mutate_all(~ifelse(is.na(.), 0, .)) %>% 
        column_to_rownames("Genotype_Group_ID") %>% 
        as.matrix()
    return(votes)
}

.find_neighbors <- function(object, filter_concordant_vertices=FALSE) {
    combined_graph <- generate_graph(object@.solve_state$unsolved_relabel_data, graph_type="combined")
    swap_cats_graph <- object@.solve_state$swap_cats_graph
    v_filtered <- V(combined_graph)
    
    ## Criteria 1: filter out vertices that have at least 1 concordant edge
    if (filter_concordant_vertices) {
        v_filtered <- NULL
        for (v in V(combined_graph)) {
            v_neighbors <- neighbors(combined_graph, v)
            concordant_edge_count <- sum(E(combined_graph)[v %--% v_neighbors]$concordant)
            if (concordant_edge_count == 0) {
                v_filtered <- c(v_filtered, v)
            }
        }
    }
    
    ## Criteria 2: filter only pairs of vertices that are within at exactly 2 edges of each other
    graph_distances <- distances(combined_graph, v=v_filtered, to=v_filtered)
    dist_equals_2 <- graph_distances == 2
    row_indices <- row(dist_equals_2)[dist_equals_2]
    col_indices <- col(dist_equals_2)[dist_equals_2]
    unique_pairs <- data.frame(
        Row = rownames(graph_distances)[row_indices],
        Col = colnames(graph_distances)[col_indices]) %>%
        transmute(
            Sample_A = pmin(Row, Col),
            Sample_B = pmax(Row, Col)
        ) %>% unique()
    
    ## Criteria 3: filter only pairs of vertices that are within the same swap category
    swap_cats_edges <- get.data.frame(swap_cats_graph) %>% 
        rename(Sample_A="from", Sample_B="to")
    unique_pairs <- inner_join(unique_pairs, swap_cats_edges, by=c("Sample_A", "Sample_B"))
    
    ## Criteria 4: filter out pairs of vertices that will violate putative_subjects
    unique_pairs %<>% left_join(
            object@.solve_state$unsolved_relabel_data[, c("Sample_ID", "Genotype_Group_ID")],
            by=c("Sample_A"="Sample_ID")
        ) %>% left_join(
            object@sample_genotype_data[, c("Sample_ID", "Subject_ID")],
            by=c("Sample_B"="Sample_ID")
        ) %>% 
        rename("new_Subject_ID"="Subject_ID") %>% 
        left_join(
            object@.solve_state$putative_subjects,
            by="Genotype_Group_ID"
        ) %>% 
        mutate(
            Valid_Swap = is.na(Subject_ID) | new_Subject_ID == Subject_ID
        ) %>% 
        filter(Valid_Swap) %>% 
        select(Sample_A, Sample_B)

    return(unique_pairs)
}

.objective_genotype_entropy = function(relabel_data) {
    relabel_data %>% 
        group_by(Genotype_Group_ID) %>%
        summarize(genotype_entropy = n()*entropy(table(Subject_ID)))  %>%
        select(genotype_entropy) %>%
        sum()
}

.objective_hamming_distance <- function(relabel_data) {
    labels_graph <- generate_graph(relabel_data, "label")
    genotype_graph <- generate_graph(relabel_data, "genotype")
    int_graph <- graph.intersection(labels_graph, genotype_graph)
    return(ecount(labels_graph) + ecount(genotype_graph) - 2*ecount(int_graph))
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
    
    duplicated_samples <- swap_cats$Sample_ID[duplicated(swap_cats$Sample_ID)]
    assert_that(
        length(duplicated_samples) == 0,
        msg = glue("'swap_cats' has non-unique Sample_ID(s) {paste(duplicated_samples, collapse=\", \")}")
    )
    
    if (ncol(swap_cats) > 1) {
        swap_cats_factor_filter <- sapply(swap_cats[, 2:ncol(swap_cats)], is.factor)
        non_factor_cols <- colnames(swap_cats[, 2:ncol(swap_cats)])[!swap_cats_factor_filter]
        assert_that(
            all(swap_cats_factor_filter),
            msg = glue("'swap_cats' columns 2:n must be of type factor, check column(s) {paste(non_factor_cols, collapse=\", \")}")
        )
    }

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
    
    extra_genotype_groups <- setdiff(putative_subjects$Genotype_Group_ID, sample_genotype_data$Genotype_Group_ID)
    assert_that(
        length(extra_genotype_groups) == 0,
        msg = glue("'putative_subjects' has 'Genotype_Group_ID'(s) not found in 'sample_genotype_data', check {paste(extra_genotype_groups, collapse=\", \")}")
    )
    extra_subjects <- setdiff(putative_subjects$Subject_ID, sample_genotype_data$Subject_ID)
    assert_that(
        length(extra_subjects) == 0,
        msg = glue("'putative_subjects' has 'Subject_ID'(s) not found in 'sample_genotype_data', check {paste(extra_subjects, collapse=\", \")}")
    )
    
    duplicated_genotype_groups <- putative_subjects$Genotype_Group_ID[duplicated(putative_subjects$Genotype_Group_ID)]
    assert_that(
        length(duplicated_genotype_groups) == 0,
        msg = glue("'putative_subjects' does not map 'Subject_ID' to 'Genotype_Group_ID' one-to-one, check 'Genotype_Group'(s) {paste(duplicated_genotype_groups, collapse=\", \")}")
    )
    duplicated_subjects <- putative_subjects$Subject_ID[duplicated(putative_subjects$Subject_ID)]
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
        pull(Sample_ID)
    genotype_inconsistent_samples <- anchor_samples_consistency %>% 
        filter(n_Genotype_Group_ID != 1) %>% 
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
    
    ## TODO: compare this with putative subjects
}

################################################################################
##################               GRAPH HELPERS                ##################
################################################################################

.plot_graph <- function(graph, to_write=FALSE) {
    edge_colors <- if (!is.null(E(graph)$edge_colors)) E(graph)$edge_colors else "grey"
    vertex_colors <- if (!is.null(V(graph)$vertex_colors)) V(graph)$vertex_colors else "orange"
    layout_custom <- with_seed(layout_nicely(graph), seed=1987)
    vertex.label.cex <- 0.6
    vertex.size <- 5
    edge.arrow.size <- 0.8
    edge.width <- 3
    if (to_write) {
        vertex.label.cex <- 2
        vertex.size <- 2
        edge.arrow.size <- 2
        edge.width <- 8
    }
    my_plot <- plot(graph, vertex.size=vertex.size, vertex.label.cex=vertex.label.cex, edge.arrow.size=edge.arrow.size, 
                    edge.width=edge.width, edge.color=edge_colors, vertex.color=vertex_colors, vertex.label.color="black", 
                    vertex.frame.color="transparent", layout=layout_custom)
    return(my_plot)
}

.generate_graph <- function(relabel_data, graph_type=c("label", "genotype", "combined")) {
    graph_type_mapping <- list(
        label = "Subject_ID",
        genotype = "Genotype_Group_ID",
        combined = NA_character_
    )
    
    graph_type <- as.character(graph_type)
    graph_type <- match.arg(graph_type)
    
    if (graph_type == "combined") {
        genotype_graph <- .generate_graph(relabel_data, "genotype")
        E(genotype_graph)$genotypes <- TRUE
        label_graph <- .generate_graph(relabel_data, "label")
        E(label_graph)$labels <- TRUE 
        graph <- graph.union(genotype_graph, label_graph, byname=TRUE)
        E(graph)[is.na(E(graph)$genotypes)]$genotypes <- FALSE
        E(graph)[is.na(E(graph)$labels)]$labels <- FALSE
        E(graph)$concordant <- E(graph)$genotypes & E(graph)$labels
    } else {
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
    
    if (graph_type == "combined") {
        E(graph)$edge_colors <- ifelse(E(graph)$concordant, "forestgreen", ifelse(E(graph)$genotypes, "orange", "cornflowerblue"))
    } else if (graph_type == "label") {
        E(graph)$edge_colors <- "cornflowerblue"
    } else {
        E(graph)$edge_colors <- "orange"
    }
    
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
    # if (ecount(relabels_graph) == 0) return(relabels)
    # print(table(components(relabels_graph)$membership))
    relabels <- data.frame(relabel_from=character(0), relabel_to=character(0))
    all_relabeled_samples <- NULL
    all_cycles <- list()
    cutoff <- 1
    while (vcount(relabels_graph) > 0 && cutoff < max(table(components(relabels_graph)$membership))) {
        curr_cycles <- find_directed_cycles(relabels_graph, cutoff=cutoff)
        for (curr_cycle in curr_cycles) {
            if (!any(curr_cycle %in% all_relabeled_samples)) {
                all_cycles <- append(all_cycles, list(curr_cycle))
                all_relabeled_samples <- c(all_relabeled_samples, curr_cycle)
                relabels_graph <- delete_vertices(relabels_graph, curr_cycle)
            }
        }
        cutoff <- cutoff + 1
    }
    
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

# my_mislabel_solver <- new("MislabelSolver", sample_genotype_data, swap_cats)
my_mislabel_solver1 <- load_test_case("1-4C-1-2C-D")
my_mislabel_solver2 <- load_test_case("1-4C")
my_mislabel_solver1 <- solve(my_mislabel_solver1)
plot(my_mislabel_solver, filter_solved=FALSE)
# plot(my_mislabel_solver, )
my_mislabel_solver <- solve_comprehensive_search(my_mislabel_solver)
plot(my_mislabel_solver, filter_solved=FALSE)
my_mislabel_solver <- solve_majority_search(my_mislabel_solver)
plot(my_mislabel_solver, filter_solved=FALSE)
my_mislabel_solver <- solve_comprehensive_search(my_mislabel_solver)
plot(my_mislabel_solver, filter_solved=FALSE)
my_mislabel_solver <- solve_majority_search(my_mislabel_solver)
plot(my_mislabel_solver, filter_solved=FALSE)
my_mislabel_solver <- solve_local_search(my_mislabel_solver, n_iter=3)
plot(my_mislabel_solver, filter_solved=TRUE)

plot(my_mislabel_solver, filter_solved=FALSE)
my_mislabel_solver <- solve_comprehensive_search(my_mislabel_solver)
plot(my_mislabel_solver, filter_solved=FALSE)
plot(my_mislabel_solver, filter_solved=TRUE)
my_mislabel_solver <- solve_comprehensive_search(my_mislabel_solver)
plot(my_mislabel_solver, filter_solved=FALSE)

