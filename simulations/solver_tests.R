n_subjects = 10
n_samples_per_subject = 10
n_swap_cats = 3
fraction_mislabel = 0.2
fraction_anchor = 0.2
fraction_ghost = 0.2
seed = 1986
output_dir = "/Users/charlesdeng/Workspace/mislabeling/test_output"
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
with_seed(seed, ghost_samples <- sample(sample_meta_data %>% pull(Sample_ID), n_ghost_samples, replace=FALSE))
sample_genotype_data[sample_genotype_data$Sample_ID %in% ghost_samples, "Genotype_Group_ID"] <- NA_character_
with_seed(seed, anchor_samples <- sample(sample_meta_data %>% filter(!Mislabeled) %>% pull(Sample_ID), n_anchor_samples, replace=FALSE))

mislabel_solver <- new("MislabelSolver", sample_genotype_data, swap_cats, anchor_samples)

object <- mislabel_solver
object <- solve_comprehensive_search(object)

object <- mislabel_solver
object <- solve_majority_search(object)

object <- mislabel_solver
object <- solve_majority_search(object, unambiguous_only=TRUE)

object <- mislabel_solver
object <- solve_local_search(object)

object <- mislabel_solver
object <- solve_local_search(object, n_iter=3)

object <- mislabel_solver
object <- solve_local_search(object, n_iter=3, include_ghost=TRUE)

object <- mislabel_solver
object <- solve_local_search(object, n_iter=3, include_ghost=TRUE, filter_concordant_vertices=TRUE)

object <- mislabel_solver
object <- solve_local_search(object)
object <- solve_local_search(object); plot(object)
object <- solve_local_search(object, n_iter=8); plot(object)
object <- solve(mislabel_solver)
plot(object, unsolved=FALSE)

