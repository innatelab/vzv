# Virginie Girault VZV project - FP MMP8 kO


Sys.setenv(MKL_NUM_THREADS = 1)

#job.args <- c("vgirault_vzvkofp", "vgirault_vzvkofp_msglm", "20230207", "6685")

if (!exists('job.args')) {
  job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message("Project ID=", project_id)
job_name <- as.character(job.args[[2]])
job_version <- job.args[[3]]
fit_version <- job_version
job_chunk <- as.integer(job.args[[4]])

message('Job ', job_name, '(id=',job_chunk,
        " fit_version=", fit_version,
        " running on ", Sys.info()["nodename"], ")")

#source('~/R/config.R')#Use this line for local calculation
source("/projects/R/config.R")#Use this line for cluster
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

rdata_filepath <- file.path(scratch_path, paste0(project_id, "_msglm_data_", fit_version, ".RData"))
message("Loading data from ", rdata_filepath)
load(rdata_filepath)

if (Sys.getenv('SLURM_CPUS_PER_TASK') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) # SLURM way
} else if (Sys.getenv('NSLOTS') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('NSLOTS')) # SGE way
} else {
  mcmc_nchains <- 8
}

require(rlang)
require(dplyr)
require(msglm)
require(tidyr)
require(stringr)

sel_object_ids <- msdata$objects$object_id[[job_chunk]]
job_chunk <- dplyr::filter(msdata$objects, object_id == sel_object_ids) %>% dplyr::pull(chunk) %>% .[[1]]
message(paste0(sel_object_ids, collapse=" "), " ", msdata$msentities[['object']], " ID(s): ",
        paste0(sort(unique(dplyr::filter(msdata$objects, object_id %in% sel_object_ids) %>% dplyr::pull(object_label))), collapse=' '))
model_data <- msglm_data(msglm_def, msdata, sel_object_ids,
                         max_quantobjects = 50L,
                         eager_msprotocols = FALSE,
                         #specificity_msexp_group_cols = "strain",
                         #cooccurrence_msexp_group_cols = "condition",
                         quantobject_fdr = 0.01)

# remove unneeded data to free some memory
msdata <- NULL
gc()

message('Running STAN in NUTS mode...')
options(mc.cores=mcmc_nchains)

msglm.stan_fit <- fit_model(model_data, adapt_delta=0.9,
                            iter=4000L, chains=mcmc_nchains, #change back to 4000
                            stanmodel_options = list(
                               object_labu_min_scale = 3,
                               batch_effect_sigma = 1.0,
                               object_msprobe_shift_tau = 0.1, object_msprobe_shift_df = 2,
                               effect_slab_scale = 1, effect_slab_df = 16, #Decrease scale and increase df to make the regularisation more stringent!!!
                               quant_batch_tau = 1.0, quant_batch_df = 2,
                               quantobject_shift_sigma = 1.0, quantobject_shift_df = 4
                            ))

#require(shinystan)
#require(rstan)
#launch_shinystan(rstan::read_stan_csv(msglm.stan_fit$output_files()))
msglm_results <- process.stan_fit(msglm.stan_fit, model_data)

res_prefix <- paste0(project_id, "_msglm")
if (!dir.exists(file.path(scratch_path, res_prefix))) {
  dir.create(file.path(scratch_path, res_prefix))
}
rfit_filepath <- file.path(scratch_path, res_prefix, paste0(res_prefix, '_', fit_version, '_', job_chunk, '.RData'))
message('Saving STAN results to ', rfit_filepath, '...')
results_info <- list(project_id = project_id, msfolder = data_info$msfolder,
                     fit_version = fit_version,
                     job_name = job_name, job_chunk = job_chunk)
save(data_info, results_info, model_data, msglm_results,
     file = rfit_filepath)
message('Done.')
on.exit(unlink(tempdir(), force = TRUE), add=TRUE)
