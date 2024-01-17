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

swap_samples <- function(sample_meta_data, swaps) {
    if (nrow(swaps) == 0) {return(sample_meta_data)}
    
    unexpected_samples <- setdiff(unique(unlist(swaps)), sample_meta_data$Sample_ID)
    assert_that(length(unexpected_samples) == 0, 
                msg = glue("Cannot swap samples {unexpected_samples} that aren't part of Sample_ID set"))
    
    for (row in 1:nrow(swaps)) {
        index1 <- match(swaps[row, 1], sample_meta_data$Sample_ID)
        index2 <- match(swaps[row, 2], sample_meta_data$Sample_ID)
        sample1 <- sample_meta_data[index1, "Sample_ID"][[1]]
        subject1 <- sample_meta_data[index1, "Subject_ID"][[1]]
        sample2 <- sample_meta_data[index2, "Sample_ID"][[1]]
        subject2 <- sample_meta_data[index2, "Subject_ID"][[1]]
        sample_meta_data[index2, c("Sample_ID", "Subject_ID")] <- c(sample1, subject1)
        sample_meta_data[index1, c("Sample_ID", "Subject_ID")] <- c(sample2, subject2)
    }
    sample_meta_data <- sample_meta_data %>% 
        mutate(
            Mislabeled = Sample_ID != E_Sample_ID,
            Mislabeled_Subject = Subject_ID != E_Subject_ID
        )
    return(sample_meta_data)
}

genswaps_random_with_replacement = function(sample_meta_data, n_mislabels, swap_cats) {
    ## Figure out the maximum number of new mislabels you can make
    ## since samples that are singletons in their swap_cats cannot be mislabeled
    swap_cats_counts <- swap_cats %>%
        group_by(SwapCat1) %>% 
        summarize(n_samples = n(), prob = n_samples/nrow(swap_cats))
    all_swap_cats <- swap_cats_counts$SwapCat1
    single_vertex_swap_cats <- swap_cats_counts[swap_cats_counts$n_samples == 1, "SwapCat1"]
    eligible_swap_cats <- all_swap_cats[!(all_swap_cats %in% single_vertex_swap_cats)]
    single_vertices <- swap_cats[swap_cats$SwapCat1 %in% single_vertex_swap_cats, "Sample_ID"]
    n_correct_labels <- nrow(sample_meta_data %>% filter(!Mislabeled))
    n_max_new_mislabels <- n_correct_labels - length(single_vertices)
    assert_that(n_mislabels <= n_correct_labels, 
                msg = glue("Can't generate {n_mislabels} new mislabels, max that can be generated is {n_max_new_mislabels}"))
    n_curr_mislabels <- nrow(sample_meta_data %>% filter(Mislabeled))
    n_goal_mislabels <- n_curr_mislabels + n_mislabels
    n_rounds <- 0
    while (n_curr_mislabels < n_goal_mislabels) {
        sample_mislabels <- sample_meta_data[sample_meta_data$Mislabeled, "Sample_ID"]
        n_new_mislabels <- n_goal_mislabels - n_curr_mislabels
        sampled_swap_cats <- sample(swap_cats_counts$SwapCat1, size=n_new_mislabels, replace=TRUE, prob=swap_cats_counts$prob)
        sample_counts <- 5*2*table(sampled_swap_cats)
        samples_to_swap <- c()
        for (swap_cat_name in names(sample_counts)) {
            num_to_sample <- sample_counts[swap_cat_name]
            new_samples_to_swap <- sample(swap_cats[swap_cats$SwapCat1 == swap_cat_name, "Sample_ID"], num_to_sample, replace=TRUE)
            samples_to_swap <- c(samples_to_swap, new_samples_to_swap)
        }
        edges <- matrix(samples_to_swap, ncol = 2, byrow = TRUE) %>% 
            as.data.frame()
        colnames(edges) <- c("from", "to")
        edges <- edges %>% 
            mutate(
                from_mislabeled = from %in% sample_mislabels,
                to_mislabeled = to %in% sample_mislabels,
                n_mislabels = from_mislabeled + to_mislabeled
            )
        both_correct <- edges %>% filter(n_mislabels == 0) %>% select(from, to)
        one_correct <- edges %>% filter(n_mislabels == 1) %>% select(from, to)
        all_edges <- edges %>% select(from, to)
        if (n_rounds > 8) {
            break
        } else if (n_new_mislabels == 1) {
            if (nrow(one_correct) < 1) {
                n_rounds <-  n_rounds + 1
                next
            }
            swaps <- one_correct[sample(nrow(one_correct), 1), ]
        } else if (n_new_mislabels == 2) {
            if (nrow(both_correct) == 0) {
                if (nrow(one_correct) < 2) {
                    n_rounds <- n_rounds + 1
                }
                swaps <- one_correct[sample(nrow(one_correct), 2), ]
            }
            swaps <- both_correct[sample(nrow(both_correct), 1), ]
        } else if (n_rounds > 5) {
            swaps <- both_correct[sample(nrow(both_correct), floor(n_new_mislabels/2), replace=FALSE), ]
        } else {
            swaps <- all_edges[sample(nrow(all_edges), floor(n_new_mislabels/2), replace=FALSE), ]
        }
        sample_meta_data <- swap_samples(sample_meta_data, swaps)
        n_rounds <- n_rounds + 1
        n_curr_mislabels <- sum(sample_meta_data$Mislabeled)
    }
    return(sample_meta_data)
}

local({
    n_subjects_per_group = 30
    n_samples_per_group = 90
    n_swap_cats = 2
    n_mislabels = 15
    seed = 1986
})

sim_mislabeled_sample_meta_data <- function (
        n_subjects_per_group,
        n_samples_per_group,
        n_swap_cats,
        n_mislabels,
        seed
) {
    if (!missing(seed)) {
        sim_params <- mget(c(
            "n_subjects_per_group",
            "n_samples_per_group",
            "n_swap_cats",
            "n_mislabels"
        ))
        res <- with_seed(seed, do.call(sim_mislabeled_sample_meta_data, sim_params))
        res$params$seed <- seed
        return(res)
    }
    
    n_samples <- 2 * n_samples_per_group
    
    assert_that(
        n_subjects_per_group > 1,
        n_samples_per_group >= n_subjects_per_group,
        n_swap_cats >= 0,
        n_mislabels >= 0 & n_mislabels != 1,
        n_mislabels <= n_samples
    )
    
    ## Generate the subjects
    ctrl_subjects <- generate_ids("Ctrl", n_subjects_per_group)
    case_subjects <- generate_ids("Case", n_subjects_per_group)
    subject_meta_table <- tibble(
        Subject_ID = c(ctrl_subjects, case_subjects),
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
    sample_meta_data <- tibble(
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
    
    swap_cats = data.frame(Sample_ID=sample_meta_data$Sample_ID)
    if (n_swap_cats > 0) {
        swap_cat_groups <- generate_ids("SwapCat_Group", n_swap_cats)
        swap_cats[["SwapCat1"]] <- as.factor(sample(swap_cat_groups, nrow(sample_meta_data), replace=TRUE))
    }
    
    sample_meta_data <- as.data.frame(sample_meta_data)
    swap_cats <- as.data.frame(swap_cats)
    
    ## Create mislabeled samples
    message(glue("Swapping samples to create mislabels"))
    sample_meta_data <- genswaps_random_with_replacement(sample_meta_data, n_mislabels, swap_cats)
    
    return(list(sample_meta_data = sample_meta_data, swap_cats = swap_cats))
}

