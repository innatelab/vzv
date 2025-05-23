# Fit assembly
#

###############################################################################

project_id <- "vgirault_vzvapms"
message('Project ID=', project_id)
data_version <- "20180301"
analysis_version <- "20180301"
message('Dataset version is ', data_version)
message('Analysis version is ', analysis_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

mop.max_nprocesses <- 8
mop.nprocesses <- 8
source(file.path(pipeline_scripts_path, 'init_cluster.R'))

require(rstan)
require(dplyr)
require(stringr)
require(msglm)
require(msimportr)

message('Loading data...')
load(file.path(scratch_path, paste0(project_id, '_msglm_data_', data_version, '.RData')))
load(file.path(scratch_path, paste0(project_id, '_msfull_', data_version, '.RData')))

ms_data$protgroups <- dplyr::mutate(ms_data$protgroups,
                                    protgroup_label = paste0(str_trunc(if_else(!is.na(gene_names), gene_names, "<noname>"), 20),
                                                             ": ", str_trunc(majority_protein_acs, 20), " (", protgroup_id, ")"))

message('Loading GLM models...')
strip_samples <- TRUE

fit_path <- file.path(scratch_path, "vzvapms_msglm")
fit_files <- list.files(fit_path, paste0('vzvapms_msglm_', analysis_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- data.frame(filename = as.character(fit_files), stringsAsFactors = FALSE) %>%
  bind_cols(
    bind_rows(lapply(strsplit(as.character(fit_files), '[_.]', fixed=FALSE), function(file_chunks) {
    data.frame(task_id = as.integer(file_chunks[[4]]),
               stringsAsFactors = FALSE)}))) %>%
  dplyr::arrange(task_id) %>%
  dplyr::mutate(protgroup_id = ms_data$protgroups$protgroup_id[task_id])

id_range_breaks <- which(c(fit_files.df$task_id[-1] - 1L, nrow(ms_data$protgroups)) !=
                         c(fit_files.df$task_id[-length(fit_files.df$task_id)], nrow(ms_data$protgroups)))
fit_files.df[sort(c(id_range_breaks, id_range_breaks+1L)), ]

if (!is.null(mop.cluster)) {
  clusterEvalQ(mop.cluster, library(dplyr))
  clusterExport(mop.cluster, varlist=c('fit_path', 'fit_files.df', 'process_msglm_chunk'))
  fit_reports <- clusterApplyLB(mop.cluster, seq_along(fit_files.df$task_id), process_msglm_chunk)
  clusterEvalQ(mop.cluster, gc())
} else {
  fit_reports <- lapply(seq_along(fit_files.df$task_id), process_msglm_chunk)
}
names(fit_reports) <- sapply(fit_reports, function(report) paste0(report$results_info$analysis_version, '_', report$model_data$objects$protgroup_id[1]))

fit_stats <- lapply(names(fit_reports[[1]]$msglm_results), join_msglm_reports, fit_reports, 'stats')
names(fit_stats) <- names(fit_reports[[1]]$msglm_results)

fit_contrasts <- lapply(names(fit_reports[[1]]$msglm_results), join_msglm_reports, fit_reports, 'contrast_stats')
names(fit_contrasts) <- names(fit_reports[[1]]$msglm_results)

source(file.path(misc_scripts_path, 'ggplot_ext.R'))



median_log2_bck_min <- 0.5
median_log2_min <- 0.5
pvalue_max <- 0.01

protgroup_contrasts.df <- fit_contrasts$iactions %>%
  dplyr::ungroup() %>%
  dplyr::filter(var=="iaction_labu_replCI") %>%
  dplyr::left_join(conditionXmetacondition.df) %>%
  dplyr::left_join(dplyr::select(conditions.df, condition, bait_orf)) %>%
  dplyr::inner_join(dplyr::select(dplyr::filter(ms_data$protgroups, !is_contaminant & !is_reverse), protgroup_id, is_viral)) %>%
  dplyr::left_join(dplyr::group_by(dplyr::inner_join(ms_data$protgroup_idents,
                                                     dplyr::select(ms_data$mschannels, msrun, condition)),
                                   condition, protgroup_id) %>%
                   dplyr::summarise(is_identified = any(ident_type == "By MS/MS", na.rm=TRUE)) %>%
                   dplyr::ungroup()) %>%
  dplyr::mutate(is_identified = !is.na(is_identified) & is_identified,
                trunc_mean = pmax(-5, pmin(5, mean)),
                trunc_median_log2 = pmax(-5, pmin(5, median_log2)),
                p_value = pmin(prob_nonpos,prob_nonneg)
  ) %>%
  dplyr::group_by(contrast, protgroup_id) %>%
  dplyr::filter(row_number() == 1L) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(p_value_adj = pmin(p.adjust(c(prob_nonpos, prob_nonneg), method = "BY")[1:n()],
                                   p.adjust(c(prob_nonneg, prob_nonpos), method = "BY")[1:n()])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(is_signif = p_value <= pvalue_max & abs(median_log2) >= median_log2_min,
                is_hit = is_identified & is_signif & median_log2 > median_log2_min,
                change = if_else(is_hit, if_else(median_log2 < 0, "-", "+"), "0"),
                contrast_type = case_when(str_detect(.$contrast, "_vs_101c") ~ "VS_GFP",
                                 str_detect(.$contrast, "_vs_allMinus") ~ "VS_ALL",
                                 TRUE ~ NA_character_)) %>%
  dplyr::arrange(bait_orf, contrast_type, p_value)

require(openxlsx)

protgroup_batch_effects.df <- dplyr::filter(fit_stats$object_batch_effects, batch_effect=="contamination") %>%
  dplyr::mutate(is_hit = p_value <= 1E-5 & median_log2 >= 2.5)
write_tsv(dplyr::select(protgroup_batch_effects.df,
                        -var_name, -index_object_batch_effect, -glm_object_ix, -batch_effect_glm_object_ix, -report_name) %>%
          dplyr::inner_join(dplyr::select(protgroup)),
          path=file.path(data_path, "csv", paste0("vgirault_vzvapms_protgroup_contamination_", analysis_version, ".txt")),
          na = "")
require(openxlsx)

protgroup_contrasts_report.df <- dplyr::select(protgroup_contrasts.df, -condition, -metacondition, -condition_role, -contrast_weight,
                        -var, -ncombn, -nsamples, -bin_width, -bw, -index_contrast, -glm_object_ix, -report_name,
                        -trunc_mean, -trunc_median_log2, -prob_nonneg, -prob_nonpos) %>%
          dplyr::semi_join(dplyr::distinct(dplyr::select(dplyr::inner_join(ms_data$protgroup_intensities, ms_data$mschannels), protgroup_id, bait_orf))) %>%
          dplyr::left_join(dplyr::select(protgroup_batch_effects.df,
                                         protgroup_id, contamination_p_value = p_value,
                                         contamination_median_log2 = median_log2,
                                         is_putative_contaminant = is_hit)) %>%
          dplyr::mutate(is_hit = is_hit & !is_putative_contaminant)
write_tsv(protgroup_contrasts_report.df,
          path=file.path(data_path, "csv", paste0("vgirault_vzvapms_contrasts_", analysis_version, ".txt")),
          na = "")
write.xlsx(protgroup_contrasts_report.df, file=file.path(data_path, "csv", paste0("vgirault_vzvapms_contrasts_", analysis_version, ".xlsx")))


protgroup_effects.df <- fit_stats$object_effects %>%
  dplyr::filter(var %in% c('obj_effect_replCI', 'obj_effect')) %>%
  dplyr::mutate(mode = if_else(var=="obj_effect", "weak", "strong")) %>%
  dplyr::inner_join(dplyr::select(dplyr::filter(ms_data$protgroups, !is_contaminant & !is_reverse), protgroup_id, protein_names, gene_names)) %>%
  dplyr::mutate(protein_names = str_replace(protein_names, ";+", ""),
                protein_names = if_else(protein_names == "", NA_character_, protein_names),
                protgroup_label = if_else(is.na(gene_names),
                                        if_else(is.na(majority_protein_acs), "NA", majority_protein_acs),
                                        gene_names) %>% str_trunc(width=20) %>% paste0(" (", protgroup_id, ")")
  ) %>%
  dplyr::group_by(mode, effect) %>%
  dplyr::mutate(p_value_adj = p.adjust(p_value, method = "BY")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(is_signif = p_value <= 3E-2, is_hit = is_signif & abs(median_log2) >= 0.1,
                change = if_else(is_signif, if_else(mean < 0, "-", "+"), ".")) %>%
  dplyr::select(-effect_glm_object_ix, -index_object_effect, -var_name, -report_name, -glm_object_ix)
protgroup_effects_weak.df <- dplyr::filter(protgroup_effects.df, mode=="weak")
protgroup_effects.df <- dplyr::filter(protgroup_effects.df, mode=="strong")

protgroup_effects.df <- dplyr::group_by(protgroup_effects.df, protgroup_id) %>%
  dplyr::filter(any(change == "+")) %>%
  dplyr::ungroup()

sel_IFN_protgroups.df <- dplyr::filter(ms_data$protgroups, str_detect(gene_names, "(^|;)IFIT\\d|DDX58|RIGI|ISG\\d|STAT\\d|IRF\\d|OAS\\d|MX1|PARP9|DTX3L"))

effects.df <- dplyr::mutate(effects.df,
                            ms_dependent = str_detect(effect, "ms"),
                            bio_dependent = str_detect(effect, "bio"),
                            treatment_dependent = str_detect(effect, "treatment"))

sel_protgroup_effects.df <- dplyr::inner_join(protgroup_effects.df, dplyr::filter(effects.df, treatment != "mock")) %>%
  dplyr::group_by(protgroup_id, cell, treatment) %>%
  # leave only hits that are batch-independent or in the batches the change is in the same direction
  dplyr::filter(any(is_hit & !bio_dependent) & !(any(change == "+") & any(change == "-"))) %>%
  dplyr::ungroup() %>%
  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
  # remove batch-dependent effects
  dplyr::anti_join(dplyr::filter(protgroup_effects.df, (effect %in% dplyr::filter(effects.df, effect %in% c("bioMar", "bioSept"))$effect)) %>%
    dplyr::group_by(protgroup_id) %>%
    dplyr::filter(any(p_value <= 5E-2 & median_log2 < 0) & row_number()==1) %>% dplyr::ungroup() %>%
    dplyr::select(protgroup_id) %>% dplyr::distinct()) %>%
  dplyr::inner_join(protgroup_effects.df) %>%
  dplyr::mutate(gene_col = as.integer(cut(as.integer(factor(protgroup_label)), 3)))

sel_protgroup_effects.df <- dplyr::filter(protgroup_effects.df,
                                          str_detect(gene_names, paste0("(^|;)(", paste0(shRNAs.df$gene_name, collapse="|"), ")"))) %>%
  dplyr::mutate(gene_col = as.integer(cut(as.integer(factor(protgroup_label)), 2)))

protgroup_effects_wide.df <- reshape(dplyr::select(protgroup_effects.df, protgroup_id, protgroup_label, gene_names, protein_names, majority_protein_acs,
                                                   effect, median_log2, p_value, change, is_signif, p_value_adj, mean_log2, sd_log2) %>%
                                         dplyr::mutate(trunc_median_log2 = pmax(-5.0, pmin(5.0, median_log2))),
                                         direction="wide",
                                         idvar="protgroup_id",
                                         timevar="effect", v.names=c("trunc_median_log2", "median_log2", "p_value", "change", "is_signif", "mean_log2", "sd_log2", "p_value_adj"))
protgroup_effects.mtx <- as.matrix(dplyr::select(protgroup_effects_wide.df, starts_with("trunc_median_log2")))
rownames(protgroup_effects.mtx) <- protgroup_effects_wide.df$protgroup_label
colnames(protgroup_effects.mtx) <- str_replace(colnames(protgroup_effects.mtx), "^trunc_median_log2\\.bait_", "")

norm_protgroup_effects.mtx <- exp(protgroup_effects.mtx) %*% diag(1/colMeans(exp(protgroup_effects.mtx)))

pdf(file = file.path(data_path, "plots", paste0("effects_heatmap2_", analysis_version, ".pdf")), width=20, height=80)
heatmap.2(norm_protgroup_effects.mtx, #Colv = NULL, Rowv = NULL,
          #distfun = function(mtx) if (nrow(mtx) == nrow(bait_intensities_orig_wide.mtx)) bait_intensities_wide.dist else bait_intensities_wide.tdist,
          trace = "none", sepwidth=c(0.05, 0), sepcolor=NA, #Colv = NULL,
          na.color = "gray", key=FALSE, cexRow = 0.2, cexCol = 1.5,
          lmat = matrix(c(0,2,0,3,1,0), ncol=2), lhei = c(0.003, 0.95, 0.01), lwid = c(0.08, 0.89),
          col = redgreen, symbreaks = TRUE)
dev.off()

write_tsv(protgroup_effects.df %>% dplyr::arrange(effect, mode, p_value),
          path=file.path(data_path, "csv", paste0("vgirault_vzvapms_effects_", analysis_version, ".txt")),
          na = "")


protgroup_contrasts_wide.df <- reshape(dplyr::select(protgroup_contrasts.df, protgroup_id, protein_names, gene_names, majority_protein_acs,
                                                     contrast, median_log2, mean_log2, sd_log2, p_value, change),
                                           direction="wide",
                                           idvar="protgroup_id",
                                           timevar="contrast", v.names=c("median_log2", "mean_log2", "sd_log2", "p_value", "change"))

require(readr)


source(file.path(misc_scripts_path, 'fasta_utils.R'))


rfit_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_fit_', analysis_version, '.RData'))
message( 'Saving full analysis results to ', rdata_filepath, '...' )
results_info <- list(project_id = project_id, data_version = data_version, analysis_version = analysis_version)
save(results_info, global_protgroup_labu_shift,
     fit_stats, fit_contrasts,
     protgroup_contrasts.df, protgroup_effects.df,
     ms_data,
     file = rfit_filepath)
message( 'Dataset analysis has been saved' )

bait_stats.df <- baits.df %>% left_join(filter(iactions_4graphml.df, src_protgroup_id != dest_protgroup_id)) %>%
    dplyr::group_by(dest_protgroup_id) %>%
    dplyr::mutate(n_baits = n_distinct(src_protgroup_id)) %>%
    dplyr::group_by(bait_orf) %>%
    dplyr::summarise(n_iactors = n_distinct(dest_protgroup_id),
                     n_unique_iactors = n_distinct(dest_protgroup_id[n_baits == 1]),
                     n_shared_iactors = n_distinct(dest_protgroup_id[n_baits > 1])) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n_iactors)) %>%
    dplyr::mutate(bait_orf = factor(bait_orf, levels=as.character(bait_orf)))

pdf(file = file.path(data_path, "plots", paste0(project_id, "_iactors_perbait_", analysis_version, ".pdf")), width=8, height=4)
ggplot(bait_stats.df) +
    geom_bar(aes(x=bait_orf, y=n_iactors), stat = "identity") +
    scale_y_log10() +
    theme_bw_ast(base_family = "", base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

protgroup_batch_shifts.df <- fit_model_stats.joined$object_batch_shifts %>%
  dplyr::filter(var == 'objXexp_batch_shift') %>%
  dplyr::inner_join(ms_data$protgroups %>% dplyr::select(protgroup_id, gene_names)) %>%
  dplyr::mutate(protgroup_label = paste0(str_trunc(if_else(is.na(gene_names), "<NA>", gene_names), 20), " (", protgroup_id, ")"))

ggplot(protgroup_batch_shifts.df) + geom_density(aes(x=median_log2, fill=mode), alpha=0.5) + scale_x_log10()

protgroup_batch_effects_wide.df <- fit_model_stats.joined$object_batch_effects %>%
  dplyr::filter(var == 'obj_batch_effect' & batch_effect == 'transient_contamination') %>%
  dplyr::mutate(protgroup_label = paste0(str_trunc(if_else(is.na(gene_names), "<NA>", gene_names), 20), " (", protgroup_id, ")")) %>%
  dplyr::select(mode, protgroup_id, protgroup_label, gene_names, mean_log2, median_log2, p_value) %>%
  reshape(direction = "wide", timevar = "mode", idvar = c("protgroup_id", "protgroup_label", "gene_names"), v.names = c("median_log2", "mean_log2", "p_value"))
