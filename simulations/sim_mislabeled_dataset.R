digits <- function(x, base = 10) {
    1 + floor(log(x, base = base))
}

generate_ids <- function(prefix, n, num_digits=digits(max(n))) {
    if (length(n) == 1) {
        if (n == 0) {
            return(NULL)
        }
        n <- seq_len(n)
    }
    digit_format <- str_c("%0", num_digits, "i")
    str_c(prefix, sprintf(digit_format, n))
}

reorder_elist_samples <- function(elist, sample_ids) {
    if (all(sample_ids == elist$targets$E_Sample_ID)) {
        return(elist)
    }
    
    assert_that(length(sample_ids) == dim(elist)[[2]], 
                msg = glue("sample_ids param length {length(sample_ids)} does not match elist sample dimension {dim(elist)[[2]]}"))
    
    reordered_elist <- elist
    reordered_elist$E <- elist$origs$E_true[, sample_ids]
    colnames(reordered_elist) <- colnames(elist)
    reordered_elist$targets %<>% mutate(
        E_Sample_ID = sample_ids,
        E_Subject_ID = str_extract(sample_ids, rex(start, not("_"))),
        ## this is hacky, will fix this when we swap labels instead of whatever it is I'm doing now
        Genotype_Group_ID = str_c("Genotype_Group", as.character(as.numeric(factor(E_Subject_ID)))),
        Mislabeled = !(Sample_ID == E_Sample_ID),
        Mislabeled_Subject = !(Subject_ID == E_Subject_ID)
    )
    return(reordered_elist)
}

swap_elist_samples <- function(elist, swaps) {
    if (nrow(swaps) == 0) {return(elist)}
    
    unexpected_samples <- setdiff(unique(unlist(swaps)), elist$targets$Sample_ID)
    assert_that(length(unexpected_samples) == 0, 
                msg = glue("Cannot swap samples {unexpected_samples} that aren't part of Sample_ID set"))
    
    sample_ids <- elist$targets$E_Sample_ID
    for (row in 1:nrow(swaps)) {
        index1 <- match(swaps[row, 1], elist$targets$Sample_ID)
        index2 <- match(swaps[row, 2], elist$targets$Sample_ID)
        sample1 <- sample_ids[[index1]]
        sample2 <- sample_ids[[index2]]
        sample_ids[[index1]] = sample2
        sample_ids[[index2]] = sample1
    }
    return(reorder_elist_samples(elist, sample_ids))
}

genswaps_random_with_replacement = function(elist, n_mislabels, swap_cats) {
    ## Figure out the maximum number of new mislabels you can make
    ## Samples that are singletons in their connected component cannot be mislabeled
    swap_cats_graph <- swap_cats_to_graph(swap_cats)
    components <- components(swap_cats_graph)
    single_vertex_components <- which(components$csize == 1)
    single_vertices <- NULL
    for (component_index in single_vertex_components) {
        single_vertices <- c(single_vertices, names(which(components$membership == component_index)))
    }
    n_correct_labels <- nrow(elist$targets %>% filter(!Mislabeled))
    n_max_new_mislabels <- n_correct_labels - length(single_vertices)
    assert_that(n_mislabels <= n_correct_labels, 
                msg = glue("Can't generate {n_mislabels} new mislabels, max that can be generated is {n_max_new_mislabels}"))
    n_curr_mislabels <- nrow(elist$targets %>% filter(Mislabeled))
    n_goal_mislabels <- n_curr_mislabels + n_mislabels
    n_rounds <- 0
    while (n_curr_mislabels < n_goal_mislabels) {
        n_new_mislabels <- n_goal_mislabels - n_curr_mislabels
        edges <- get.data.frame(swap_cats_graph) %>% 
            left_join(elist$targets[, c("Sample_ID", "Mislabeled")], by=c("from"="Sample_ID")) %>% 
            left_join(elist$targets[, c("Sample_ID", "Mislabeled")], by=c("to"="Sample_ID")) %>% 
            mutate(
                n_mislabels = Mislabeled.x + Mislabeled.y
            )
        both_correct <- edges %>% filter(n_mislabels == 0) %>% select(from, to)
        one_correct <- edges %>% filter(n_mislabels == 1) %>% select(from, to)
        all_edges <- edges %>% select(from, to)
        if (n_rounds > 8) {
            break
        } else if (n_new_mislabels == 1) {
            if (nrow(one_correct) < 1) break
            swaps <- one_correct[sample(nrow(one_correct), 1), ]
        } else if (n_new_mislabels == 2) {
            if (nrow(both_correct) == 0) {
                if (nrow(one_correct) < 2) break
                swaps <- one_correct[sample(nrow(one_correct), 2), ]
            }
            swaps <- both_correct[sample(nrow(both_correct), 1), ]
        } else if (n_rounds > 5) {
            swaps <- both_correct[sample(nrow(both_correct), floor(n_new_mislabels/2), replace=TRUE), ]
        } else {
            swaps <- all_edges[sample(nrow(all_edges), floor(n_new_mislabels/2), replace=TRUE), ]
        }
        elist <- swap_elist_samples(elist, swaps)
        n_rounds <- n_rounds + 1
        n_curr_mislabels <- sum(elist$targets$Mislabeled)
    }
    return(elist)
}

swap_cats_to_graph <- function(swap_cats) {
    n <- nrow(swap_cats)
    if (ncol(swap_cats) < 2) {
        swap_cats_matrix <- matrix(TRUE, nrow=n, ncol=n)
    } else {
        swap_cats_matrix <- matrix(FALSE, nrow=n, ncol=n)
        swap_cats_int <- interaction(swap_cats[, 2:ncol(swap_cats)], sep=":")
        for (level in levels(swap_cats_int)) {
            idx <- swap_cats_int == level
            swap_cats_matrix[idx, idx] <- TRUE
        }
    }
    rownames(swap_cats_matrix) <- colnames(swap_cats_matrix) <- swap_cats$Sample_ID
    return(graph_from_adjacency_matrix(swap_cats_matrix, mode="directed"))
}

local({
    n_subjects_per_group = 30
    n_samples_per_group = 150
    n_swap_cats = 2
    n_features = 2
    n_mislabels = 100
    n_sv = 2
    fraction_de_case = 0.5
    case_sd = 0.1
    subject_sd = 0.1
    sv_sd = 0.1
    resid_sd = 0.1
    sex_sd = 0.1
    age_sd = 0.1
    seed = 1986
})

sim_mislabeled_data <- function (
        n_subjects_per_group,
        n_samples_per_group,
        n_swap_cats,
        n_features,
        n_mislabels,
        n_sv,
        fraction_de_case,
        case_sd,
        sv_sd,
        subject_sd,
        resid_sd,
        sex_sd = 0,
        age_sd = 0,
        seed
) {
    if (!missing(seed)) {
        sim_params <- mget(c(
            "n_subjects_per_group",
            "n_samples_per_group",
            "n_swap_cats",
            "n_features",
            "n_mislabels",
            "n_sv",
            "fraction_de_case",
            "case_sd",
            "subject_sd",
            "sv_sd",
            "resid_sd",
            "sex_sd",
            "age_sd"
        ))
        res <- with_seed(seed, do.call(sim_mislabeled_data, sim_params))
        res$params$seed <- seed
        return(res)
    }
    
    n_samples <- 2 * n_samples_per_group
    
    assert_that(
        n_subjects_per_group > 1,
        n_samples_per_group >= n_subjects_per_group,
        n_swap_cats >= 0,
        n_features > 0,
        n_mislabels >= 0 & n_mislabels != 1,
        n_mislabels <= n_samples,
        n_sv >= 0,
        fraction_de_case >= 0,
        fraction_de_case <= 1,
        case_sd >= 0,
        subject_sd >= 0,
        sv_sd >= 0,
        resid_sd >= 0,
        sex_sd >= 0,
        age_sd >= 0
    )
    if (n_sv < 1 && sv_sd > 0) {
        warning("Setting SV signal to zero since there are not SVs")
        n_sv <- 1
        sv_sd <- 0
    }
    
    ## Generate the subjects
    ctrl_subjects <- generate_ids("Ctrl", n_subjects_per_group)
    case_subjects <- generate_ids("Case", n_subjects_per_group)
    subject_meta_table <- tibble(
        Subject_ID = c(ctrl_subjects, case_subjects),
        Status = c(rep_along(ctrl_subjects, "Control"), rep_along(case_subjects, "Case")) %>%
            factor(levels = c("Control", "Case"))
    ) %>% 
        mutate(
            Sex = sample(c("Male", "Female"), size=n(), replace=TRUE),
            Age = rbinom(n=n(), size=60, prob=0.2) + 20,
            n_samples = 1
        )
    ## Determine number of samples per subject
    addl_samples <- c(
        table(sample(ctrl_subjects, n_samples_per_group - n_subjects_per_group, replace = TRUE)),
        table(sample(case_subjects, n_samples_per_group - n_subjects_per_group, replace = TRUE))) %>% 
        as.data.frame() %>% 
        rename(c("addl_samples" = "."))
    addl_samples$Subject_ID <- rownames(addl_samples)
    subject_meta_table %<>%
        left_join(addl_samples, by="Subject_ID") %>% 
        replace_na(list(addl_samples=0)) %>% 
        mutate(
            n_samples = n_samples + addl_samples
        ) %>% select(-addl_samples)
    
    ## Generate the samples
    genotype_digit_format <- str_c("%0", digits(2*n_subjects_per_group), "i")
    sample_meta_table <- tibble(
        Subject_ID = rep.int(subject_meta_table$Subject_ID, times=subject_meta_table$n_samples)
    ) %>% 
        group_by( Subject_ID ) %>% 
        mutate(
            Sample_ID = generate_ids(str_c(Subject_ID, "_Sample"), n(), num_digits=digits(max(subject_meta_table$n_samples)))
        ) %>% 
        ungroup() %>% 
        left_join(
            subject_meta_table, by="Subject_ID"
        ) %>% 
        mutate(
            E_Sample_ID = Sample_ID,
            E_Subject_ID = Subject_ID,
            Genotype_Group_ID = str_c("Genotype_Group", sprintf(genotype_digit_format, (as.numeric(factor(E_Subject_ID))))),
            Mislabeled = FALSE,
            Mislabeled_Subject = FALSE
        ) %>% 
        select(-n_samples) %>% 
        relocate(Sample_ID)
    
    assert_that(length(unique(sample_meta_table$Sex)) > 1, msg = "Only one sex simulated")
    
    swap_cats = data.frame(Sample_ID=sample_meta_table$Sample_ID)
    if (n_swap_cats > 0) {
        swap_cat_groups <- generate_ids("SwapCat_Group", n_swap_cats)
        swap_cats[["SwapCat1"]] <- as.factor(sample(swap_cat_groups, nrow(sample_meta_table), replace=TRUE))
    }
    
    sv_names <- generate_ids("SV", n_sv)
    for (svi in sv_names) {
        sample_meta_table[[svi]] <- rnorm(nrow(sample_meta_table))
    }
    
    n_de_features <- round(n_features * fraction_de_case)
    if (n_de_features == 0 && fraction_de_case != 0) {
        warning("fraction_de_case is too small, no features are DE")
    }
    
    ## Generate the features
    feature_meta_table <- tibble(
        Feature_Number = seq_len(n_features),
        DE = Feature_Number <= n_de_features,
        Feature_ID = generate_ids("Feature", n_features) %>%
            str_c(if_else(DE, "_DE", ""))
    ) %>%
        select(Feature_ID, Feature_Number, DE, everything())
    
    ## Generate the coefs
    design_case <- model.matrix(~1 + Status, sample_meta_table)[, 2, drop = FALSE] %>%
        ## Convert to sum-to-zero
        subtract(mean(unique(.)))
    design_sex <- model.matrix(~1 + Sex, sample_meta_table)[, 2, drop = FALSE] %>%
        subtract(mean(unique(.)))
    design_age <- as.matrix(sample_meta_table["Age"]) %>%
        ## Make sure that age design has same scale as other SVs (mean=0, sd=1)
        scale()
    design_subject <- model.matrix(~0 + Subject_ID, sample_meta_table)
    design_sv <- as.matrix(sample_meta_table[sv_names])
    ## Use abs so that all true logFC are positive
    eta_case <- design_case %*% matrix(ncol = n_features, nrow = ncol(design_case), abs(rnorm(n_features * ncol(design_case)))) %>%
        ## Set all eta to 0 for non-DE genes
        scale(center = FALSE, scale = if_else(feature_meta_table$DE, 1, Inf))
    eta_sex <- design_sex %*% matrix(ncol = n_features, nrow = ncol(design_sex), rnorm(n_features * ncol(design_sex)))
    eta_age <- design_age %*% matrix(ncol = n_features, nrow = ncol(design_age), rnorm(n_features * ncol(design_age)))
    eta_subject <- design_subject %*% matrix(ncol = n_features, nrow = ncol(design_subject), rnorm(n_features * ncol(design_subject)))
    eta_sv <- design_sv %*% matrix(ncol = n_features, nrow = ncol(design_sv), rnorm(n_features * ncol(design_sv)))
    eta_resid <- matrix(ncol = n_features, nrow = nrow(sample_meta_table), rnorm(n_features * nrow(sample_meta_table)))
    
    ## Put all the expression effects together
    expr_mat <- t(
        eta_case * case_sd +
            eta_sex * sex_sd + 
            eta_age * age_sd +
            eta_subject * subject_sd + 
            eta_sv * sv_sd +
            eta_resid * resid_sd
    ) %>%
        set_colnames(sample_meta_table$Sample_ID) %>%
        set_rownames(feature_meta_table$Feature_ID)
    attr(expr_mat, "scaled:scale") <- NULL
    
    ## Build the EList
    elist <- new("EList", list(
        E = expr_mat,
        genes = as.data.frame(feature_meta_table),
        targets = as.data.frame(sample_meta_table),
        swap_cats = as.data.frame(swap_cats),
        other = list(
            eta_sv = t(eta_sv),
            eta_case = t(eta_case),
            eta_sex = t(eta_sex),
            eta_age = t(eta_age),
            eta_resid = t(eta_resid)
        ),
        params = mget(c(
            "n_subjects_per_group",
            "n_samples_per_group",
            "n_swap_cats",
            "n_features",
            "n_mislabels",
            "n_sv",
            "fraction_de_case",
            "case_sd",
            "subject_sd",
            "sv_sd",
            "resid_sd",
            "sex_sd",
            "age_sd"
        )),
        ## Store true expressions for reference when relabeling samples
        origs = list(
            E_true = expr_mat
        )
    ))
    
    ## Discard unused eta matrices
    if (case_sd == 0 || n_de_features == 0) {
        elist$other$eta_case <- NULL
    }
    sd_eta_pairs = list(c("sex_sd", "eta_sex"),
                        c("age_sd", "eta_age"),
                        c("sv_sd", "eta_sv"),
                        c("resid_sd", "eta_resid"))
    for (pair in sd_eta_pairs) {
        var_sd <- pair[[1]]
        eta_var <- pair[[2]]
        if (var_sd == 0) {
            elist$other[[eta_var]] <- NULL
        }
    }
    
    ## Apply dimnames to all relevant sub-elements in the EList
    dimnames(elist) <- dimnames(elist)
    
    ## Create mislabeled samples
    message(glue("Swapping samples to create mislabels"))
    elist <- genswaps_random_with_replacement(elist, n_mislabels, swap_cats)
    
    ## Store original targets so we can access initial Sample_ID/E_Data 
    ## pairs to assess the results of the label correction algorithm
    elist$origs$targets_orig <- elist$targets
    
    return(elist)
}

