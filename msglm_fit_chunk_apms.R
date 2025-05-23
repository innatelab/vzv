# Msglm fit chunk for the apms dataset - VZV project - Virginie Girault

#job.args <- c("vgirault_vzvapms", "vgirault_vzvapms_msglm", "20180301", "0", "3506") 

if ( !exists( 'job.args' ) ) {
  job.args <- commandArgs( trailingOnly = TRUE )
}

project_id <- job.args[[1]]
message( 'Project ID=', project_id )

job_name <- as.character( job.args[[2]] )
message( 'Job name is ', job_name )

job_version <- job.args[[3]]
data_version <- job_version
analysis_version <- job_version
message( 'Job dataset version is ', data_version, " analysis version is ", analysis_version )

job_id <- as.integer( job.args[[4]] )
message( 'Job ID is ', job_id )

job_chunk <- as.integer(job.args[[5]])
message( 'Job chunk is ', job_chunk )

source("~/R/config.R")
source( file.path( base_scripts_path, 'misc/setup_base_paths.R' ) )
source( file.path( base_scripts_path, 'misc/setup_project_paths.R' ) )
source( file.path( base_scripts_path, 'msglm/R/stan_process_utils.R') )
source( file.path( base_scripts_path, 'msglm/R/stan_models.R' ) )
source( file.path( base_scripts_path, 'msglm/R/model_data.R' ) )
source( file.path( base_scripts_path, 'msglm/R/msglm_results.R' ) )

# Read the msglm_data 
data_path <- file.path(base_data_path, project_id)
rdata_filepath <- file.path(scratch_path, paste0( project_id, '_msglm_data_', data_version, '.RData' ) )
message('Loading data from ', rdata_filepath)
load( rdata_filepath )

mop.max_nprocesses <- 8

require(dplyr)
require(rstan)

# Extract data of a given protein 
sel_protgroup_ids <- ms_data$protgroups$protgroup_id[[job_chunk]]
message(sel_protgroup_ids, " protgroup ID(s): ",
        paste0(sort(unique(dplyr::filter(ms_data$protgroups, protgroup_id %in% sel_protgroup_ids) %>%
                             .$gene_names)), collapse=' '))
ms_data.df <- dplyr::filter(ms_data$protgroup_intensities, protgroup_id %in% sel_protgroup_ids) %>%
  dplyr::inner_join(ms_data$mschannels %>% dplyr::select(msrun, condition) %>% dplyr::distinct())# %>%

message('Preparing MS GLM data...')
model_data <- list()
model_data$mschannels <- dplyr::select(dplyr::filter(ms_data$mschannels, quant_type != "aggregate"),
                                       bait_orf, bio_batch, contamination, condition, msrun_ix, msrun) %>%
  dplyr::inner_join(total_msrun_shifts.df) %>%
  dplyr::arrange(msrun_ix) %>% dplyr::distinct() %>%
  dplyr::mutate(msrun_ix = row_number())
experiment_shift_col <- 'total_msrun_shift'
model_data$mschannels$model_mschannel_shift <- model_data$mschannels[[experiment_shift_col]]
model_data$conditions <- dplyr::select(model_data$mschannels, condition, bait_orf) %>%
  dplyr::distinct() %>%
  dplyr::arrange(bait_orf) %>%
  dplyr::mutate(condition_ix = as.integer(factor(condition, levels = rownames(conditionXeffect.mtx)))) %>%
  dplyr::arrange(condition_ix)
model_data$mschannels <- dplyr::inner_join(model_data$mschannels, model_data$conditions) %>%
  dplyr::arrange(msrun_ix)

model_data$interactions <- expand.grid(protgroup_id = sel_protgroup_ids,
                                       condition_ix = model_data$conditions$condition_ix) %>%
    dplyr::inner_join(dplyr::select(model_data$conditions, condition_ix, condition) %>% dplyr::distinct()) %>%
    dplyr::left_join(ms_data.df %>% dplyr::select(condition, protgroup_id) %>% dplyr::distinct() %>%
                       dplyr::mutate(is_virtual = FALSE)) %>%
    dplyr::mutate(is_virtual = is.na(is_virtual),
                  iaction_id = paste(condition, protgroup_id))

model_data$interactions <- dplyr::arrange(model_data$interactions, condition_ix, protgroup_id) %>%
  dplyr::mutate(glm_iaction_ix = row_number(),
                glm_object_ix = as.integer(factor(protgroup_id)))

model_data$objects <- dplyr::inner_join(ms_data$protgroups %>% dplyr::select(protgroup_id, majority_protein_acs, protein_acs, contains("is_"), gene_names, protein_names),
                                        dplyr::select(model_data$interactions, glm_object_ix, protgroup_id)) %>%
  dplyr::distinct() %>%
  dplyr::arrange(glm_object_ix)

# entries for an interaction in all replicate experiments
model_data$ms_data <- dplyr::inner_join(model_data$interactions, model_data$mschannels) %>%
  dplyr::left_join(ms_data.df) %>%
  dplyr::arrange(msrun_ix, glm_object_ix) %>%
  dplyr::mutate(glm_observation_ix = seq_len(n()),
                qdata_ix = if_else(!is.na(intensity), cumsum(!is.na(intensity)), NA_integer_),
                mdata_ix = if_else(is.na(intensity), cumsum(is.na(intensity)), NA_integer_))

model_data <- prepare_effects(model_data, underdefined_iactions=TRUE)

msglm.stan_data <- stan.prepare_data(protgroup_instr_calib, model_data, batch_tau=0.8, repl_tau=0.4)

message( 'Running STAN in HMC mode...' )
rstan_options(auto_write=TRUE)
options(mc.cores=mop.max_nprocesses)
msglm.stan_fit <- stan.sampling(msglm.stan_data, adapt_delta=0.9)

min.iteration <- as.integer(1.5 * msglm.stan_fit@sim$warmup)
dims_info <- msglm.prepare_dims_info(model_data, object_cols=c('protgroup_id', "majority_protein_acs", "gene_names"))

background_contrasts <- unique(as.character(filter(contrastXmetacondition.df,
                                                   str_detect(metacondition, "allMinus") & weight < 0)$contrast))
background_contrasts.quantiles_rhs <- lapply(background_contrasts, function(contr) c(0.25, 0.9))
names(background_contrasts.quantiles_rhs) <- background_contrasts

msglm_results <- process.stan_fit(msglm.stan_fit, dims_info,
                                  condition.quantiles_rhs = background_contrasts.quantiles_rhs)

res_prefix <- "vzvapms_msglm"
if (!dir.exists(file.path(scratch_path, res_prefix))) {
  dir.create(file.path(scratch_path, res_prefix))
}
rfit_filepath <- file.path(scratch_path, res_prefix, paste0(res_prefix, '_', analysis_version, '_', job_chunk, '.RData'))
message('Saving STAN results to ', rfit_filepath, '...')
results_info <- list(project_id = project_id, data_version = data_version, analysis_version = analysis_version,
                     job_name = job_name, job_chunk = job_chunk)
save(data_info, results_info,
     model_data, msglm.stan_data, msglm_results,
     dims_info, file = rfit_filepath)
message( 'Done.' )
on.exit(unlink(tempdir(), force = TRUE), add=TRUE)

