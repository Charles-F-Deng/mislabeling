random = function(elist, n_mislabels) {
    correct_ctrl_samples <- elist$targets %>% filter(Status == "Control" & !Mislabeled) %>% pull(Sample_ID)
    correct_case_samples <- elist$targets %>% filter(Status == "Case" & !Mislabeled) %>% pull(Sample_ID)
    n_swaps <- floor(n_mislabels / 2)
    swaps <- cbind(
        sample(correct_ctrl_samples, n_swaps),
        sample(correct_case_samples, n_swaps)
    )
    return(swap_elist_samples(elist, swaps))
}

random_with_replacement = function(elist, n_mislabels) {
    n_correct_labels <- nrow(elist$targets %>% filter(!Mislabeled))
    assert_that(n_mislabels <= n_correct_labels, 
                msg = glue("Can't generate {n_mislabels} new mislabels when there are only {n_correct_samples} correctly labeled samples"))
    n_curr_mislabels <- nrow(elist$targets %>% filter(Mislabeled))
    n_goal_mislabels <- n_curr_mislabels + n_mislabels
    n_rounds <- 0
    while (n_curr_mislabels < n_goal_mislabels) {
        n_new_mislabels <- n_goal_mislabels - n_curr_mislabels
        mislabeled_samples <- elist$targets %>% filter(Mislabeled) %>% pull(Sample_ID)
        correct_samples <- elist$targets %>% filter(!Mislabeled) %>% pull(Sample_ID)
        all_samples <- c(mislabeled_samples, correct_samples)
        if (n_new_mislabels == 1) {
            swaps <- cbind(
                sample(mislabeled_samples, 1),
                sample(correct_samples, 1)
            )
        } else if (n_new_mislabels == 2) {
            swaps <- matrix(sample(correct_samples, 2), ncol = 2)
        } else if (n_rounds > 5) {
            swaps <- cbind(
                sample(mislabeled_samples, n_new_mislabels),
                sample(correct_samples, n_new_mislabels)
            )
        } else {
            n_swaps <- floor(n_new_mislabels / 2)
            swaps <- matrix(sample(all_samples, 2 * n_swaps, replace = TRUE), ncol = 2)
        }
        elist <- swap_elist_samples(elist, swaps)
        n_rounds <- n_rounds + 1
        n_curr_mislabels <- nrow(elist$targets %>% filter(Mislabeled))
    }
    return(elist)
}

systematic = function(elist, n_mislabels) {
    correct_samples <- elist$targets %>% filter(!Mislabeled) %>% pull(Sample_ID)
    n_swaps <- floor(n_mislabels / 2)
    swaps <- cbind(
        swaps <- matrix(sample(correct_samples, n_swaps * 2), ncol = 2)
    )
    return(swap_elist_samples(elist, swaps))
}

handle_swaps_options <- ls()
