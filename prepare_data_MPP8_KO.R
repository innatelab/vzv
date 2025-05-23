# Virginie Girault VZV project - FP MMP8 kO
# Reshape the Maxquant output 
# Define the linear model, effect, contrasts, priors
# Calculate normalization shifts per sample 

project_id <- "vgirault_vzvkofp"
data_version <- "20230207"
msfolder <- "mq20230207"
message('Project ID=', project_id, " data version=", data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(msimportr)
require(devtools)
require(Biostrings)
require(msglm)
require(dplyr)
require(jsonlite)
require(stringr)
require(readr)
require(pheatmap)

require(cmdstanr)

data_info <- list(project_id = project_id, data_ver=data_version, fit_ver=data_version,
                  pepmodstate_mscalib_filename = "mscalib_QEP5_intensity_pepmodstate_cov2_20211108.json",
                  msfolder = msfolder, quant_type = "LFQ", quant_col_prefix = "LFQ_Intensity",
                  pep_quant_type = "intensity")
msdata_path <- file.path(data_path, msfolder)

# import calibration file
message('Loading MS instrument protein calibration data from ', data_info$pepmodstate_mscalib_filename, '...')
pepmodstate_mscalib <- read_mscalib_json(file.path(data_path, data_info$pepmodstate_mscalib_filename))

ignored_rawfiles = c("20191125_QEp5_AnPi_SA_293T_MPP8_191206205800")

# rawfile (msrun)  annotation
msruns.df <- read_tsv(file.path(msdata_path, "experimentalDesignTemplate.txt")) %>%
  tidyr::extract("Experiment", c('treatment','replicate'),
                 "([^_]+)_([^_]+)$", remove = FALSE) %>%
  tidyr::extract("Name", c('msrun_date', 'instrument'), "^(\\d+)_(QEp\\d+)_",
                 remove = FALSE, convert=FALSE)%>%
 dplyr::rename("raw_file" = "Name") %>%
  dplyr::rename(msexperiment = Experiment) %>%
  dplyr::mutate(treatment = case_when(treatment == "P"& str_detect(msexperiment, "SKN_MPP8_P") ~ "MPP8_P",
                                      TRUE ~ treatment)) %>%
  dplyr::mutate(msrun_ix = row_number(),
                replicate = parse_integer(replicate),
                instrument = instrument,
                treatment = factor(treatment), 
                msprotocol = "FP", 
                is_skipped = raw_file %in% ignored_rawfiles) %>%
  dplyr::mutate(msexperiment = factor(msexperiment, levels=unique(msexperiment)),
                msrun = msexperiment,
                msrun = factor(msrun, levels=msrun),
                treatment = relevel(treatment, "NTC"),
                condition =  paste0(treatment),
                condition = factor(condition, levels=unique(condition)))

# import fasta files used for quanting
fasta.dfs <- list(
  human = bind_rows(read_innate_uniprot_fasta(file.path(data_path, "fasta", "UP000005640_9606.fasta")), read_innate_uniprot_fasta(file.path(data_path, "fasta", "UP000005640_9606_additional.fasta"))))

# import Proteingroup.txt
msdata.wide <- read.MaxQuant.ProteinGroups(file.path(msdata_path), import_data = c(data_info$quant_type, "ident_type"))

uniprot_fasta_extract_ac <- function(fasta_headers) {
 str_replace_all(fasta_headers, "(^|;)(?:sp|tr)\\|([a-zA-Z0-9_-]+)\\|\\S+?(?=$|;)", "\\1\\2")
}

msdata_colgroups <- attr(msdata.wide, "column_groups")

# import evidence.txt
mqevidence <- read.MaxQuant.Evidence(file.path(msdata_path),
                                     mschannel_annotate.f = function(mschans_df) {
                                       res <- dplyr::inner_join(dplyr::select(mschans_df, mstag, rawfile),
                                                                dplyr::select(msruns.df, raw_file, condition, treatment, batch,
                                                                              msexperiment, replicate, instrument, msprotocol,
                                                                              msrun, is_skipped),
                                                                by=c("rawfile"=  "raw_file"))
                                       attr(res, "column_scopes") <- c(treatment = "msexperiment",
                                                                       condition = "msexperiment",
                                                                       replicate = "msexperiment",
                                                                       instrument = "msexperiment",
                                                                       batch = "msexperiment")
                                       return(res)
                                     })

mqevidence$peptides <- read.MaxQuant.Peptides(file.path(msdata_path), file_name='peptides.txt',
                                              import_data='ident_type')
mqevidence$peaks <- NULL # exclude big data frame

strlist_label <- function(strs) {
  str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
  sapply(str_split(strs, delim), strlist_label)
}

# all ms data
msdata_full <- mqevidence[c('msexperiments', 'msruns')]

msdata_full <- append_protgroups_info(msdata_full, msdata.wide, fix_protein_info = TRUE,
                                      proteins_info = dplyr::bind_rows(
                                        dplyr::mutate(fasta.dfs$human, is_viral = FALSE)) %>%
                                        dplyr::mutate(protein_ac_noiso = str_remove(protein_ac, "-\\d+(?:#.+)*$"),
                                                      protein_isoform_ix = replace_na(as.integer(str_match(protein_ac, "-(\\d+)$")[, 2]), 1L)),
                                      import_columns = c("is_viral", "organism"))

msdata_full$peptides <- mqevidence$peptides %>%
  dplyr::left_join(select(dplyr::bind_rows(fasta.dfs), lead_razor_protein_ac = protein_ac, organism)) %>%
  dplyr::mutate(peptide_rank = 1L)
msdata_full$pepmods <- mqevidence$pepmods %>%
  dplyr::mutate(pepmod_rank = 1L)
msdata_full$pepmodstates <- mqevidence$pepmodstates

# redefine protein groups (protregroups)
peptides.df <- dplyr::select(msdata_full$peptides, peptide_id, protgroup_ids, protein_acs, lead_razor_protein_ac,
                             peptide_seq, is_reverse, peptide_rank)

proteins.df <- msdata_full$proteins
save(file = file.path(msdata_path, str_c(project_id, "_", msfolder, '_', data_version, "_peptides.RData")),
     peptides.df, proteins.df)

# .. run protregroup.jl
msdata_full$protregroups <- read_tsv(file.path(data_path, msfolder,
                                               str_c(project_id, "_", msfolder, '_', data_version, "_protregroups.txt")),
                                     col_types = list(protregroup_id = "i"))


uniprot_fasta_extract_ac <- function(fasta_headers) {
  str_replace_all(fasta_headers, "(^|;)(?:sp|tr)\\|([a-zA-Z0-9_-]+)\\|\\S+?(?=$|;)", "\\1\\2")
}

msdata_full$protregroups <- dplyr::mutate(msdata_full$protregroups, 
                                          protein_acs = uniprot_fasta_extract_ac(protein_acs), 
                                          majority_protein_acs = uniprot_fasta_extract_ac(majority_protein_acs)) 


msdata_full$protein2protregroup <- dplyr::select(msdata_full$protregroups, protregroup_id, protein_ac=majority_protein_acs) %>%
  separate_rows(protein_ac, sep=fixed(";"), convert=TRUE) %>%
  dplyr::mutate(is_majority = TRUE) %>%
  dplyr::group_by(protregroup_id) %>%
  dplyr::mutate(protein_ac_rank = row_number()) %>%
  dplyr::ungroup()

msdata_full$protregroup2peptide <- bind_rows(
  select(msdata_full$protregroups, protregroup_id, peptide_id=spec_peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = TRUE),
  select(msdata_full$protregroups, protregroup_id, peptide_id=peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = FALSE)) %>%
  dplyr::group_by(protregroup_id, peptide_id) %>%
  dplyr::summarise(is_specific = any(is_specific)) %>%
  dplyr::ungroup()
msdata_full$protregroup2pepmod <- dplyr::inner_join(msdata_full$protregroup2peptide,
                                                    dplyr::select(msdata_full$pepmods, pepmod_id, peptide_id)) %>%
  dplyr::select(-peptide_id)

msdata_full$protregroups <- msdata_full$protregroups %>%
  dplyr::mutate(gene_label = strlist_label2(gene_names),
                protein_label = strlist_label2(protein_names),
                protein_description = strlist_label2(protein_descriptions),
                is_viral = replace_na(protein_acs %in% fasta.dfs$mcherry$protein_ac, FALSE),
                protac_label = strlist_label2(majority_protein_acs),
                protregroup_label = case_when(is_viral ~ protein_label,
                                              !is.na(gene_label) ~ gene_label,
                                              !is.na(protac_label) ~ protac_label,
                                              TRUE ~ str_c('#', protregroup_id))) %>%
  dplyr::left_join(dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmods) %>%
                   dplyr::group_by(protregroup_id) %>%
                   dplyr::summarise(npeptides = n_distinct(peptide_id),
                                    npepmods = n_distinct(pepmod_id),
                                    npeptides_unique = n_distinct(peptide_id[is_specific]),
                                    npepmods_unique = n_distinct(pepmod_id[is_specific])) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(npeptides_unique_razor = npeptides_unique,
                                 npeptides_razor = 0L,
                                 npepmods_unique_razor = npepmods_unique,
                                 npepmods_razor = 0L))

msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  dplyr::mutate(is_idented = str_detect(ident_type, "MSMS"))

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod,
                                                    msdata_full$pepmodstate_intensities) %>%
  dplyr::group_by(msexperiment, protregroup_id) %>%
  dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                   nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                   npepmods_idented = sum(is_idented),
                   nspecpepmods_idented = sum(is_idented & is_specific)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))


pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  select(pepmod_id, pepmodstate_id, msrun, intensity, ident_type) %>%
  dplyr::mutate(is_idented = ident_type %in% c("ISO-MSMS", "MULTI-MSMS", "MSMS"))

msdata$protgroup_intensities <- protgroup_intensities

# setup experimental design matrices (KO)
conditions.df <- bind_rows(dplyr::select(msdata_full$msexperiments, condition, treatment, batch) %>%
                           dplyr::distinct() %>%
                           dplyr::mutate(is_virtual = FALSE)) %>%
  dplyr::arrange(is_virtual, condition) %>%
  dplyr::mutate(is_negctrl = (treatment == 'NTC')) %>% 
  dplyr::mutate(is_background = treatment == "NTC")

conditionXeffect.mtx <- model.matrix(~ 1 + treatment, conditions.df)

conditionXeffect.mtx <- conditionXeffect.mtx[, colSums(conditionXeffect.mtx != 0) > 0]
dimnames(conditionXeffect.mtx) <- list(condition = conditions.df$condition,
                                       effect = colnames(conditionXeffect.mtx))
conditionXeffect.mtx <- conditionXeffect.mtx[, setdiff(colnames(conditionXeffect.mtx), '(Intercept)'), drop=FALSE] # handled separately
  
all_conditions <- as.character(conditions.df$condition)

if (!dir.exists(plot_path)) dir.create(plot_path)
pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, na_col = "gray", border_color = "white",
         filename = file.path(plot_path,
                              paste0(project_id, "_conditionXeffect_", data_info$msfolder, "_", data_info$fit_ver, ".pdf")),
         width = 20, height = 12)

# set of priors 
effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  dplyr::mutate(treatment = effect_factor(effect, "treatment", levels(as.factor(conditions.df$treatment)), NA_character_),
                is_positive = FALSE) %>%
  dplyr::mutate(effect_type = "treatment") %>%
  dplyr::mutate(prior_tau =  1.0,
                prior_mean = 0.0,
                prior_df1 =  2.0,
                prior_df2 = prior_df1)


# combine msglm model 
msglm_def <- msglm_model(conditionXeffect.mtx, conditions.df, effects.df, verbose=TRUE) 

# combine msdata 
msdata <- import_msglm_data(msdata_full, msglm_def,
                            object="protregroup", quantobject="pepmodstate",
                            verbose = TRUE)

# normalize experiments:
msdata4norm.df <- msdata_full$pepmodstate_intensities %>%
    dplyr::semi_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id)) %>%
    dplyr::semi_join(dplyr::semi_join(msdata_full$pepmods,
        dplyr::filter(msdata_full$protregroups, !is_reverse & !is_contaminant & !is_viral) %>%
        dplyr::inner_join(msdata_full$protregroup2pepmod) %>% dplyr::select(pepmod_id)) %>%
        dplyr::filter(!is_reverse & !is_contaminant) %>%
        dplyr::select(pepmod_id))

require(cmdstanr)
options(mc.cores=8L)
Sys.setenv(MKL_NUM_THREADS=1L)

# normalize experiments:
msruns_hnorm <- multilevel_normalize_experiments(msdata$pepmodstate_mscalib,
                                                 dplyr::select(msdata$msruns, msrun, condition, treatment, msprotocol),
                                                 msdata4norm.df,
                                                 quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
                                                 mcmc.iter = 1000L,
                                                 verbose=TRUE,
                                                 stan_method = "mcmc",
                                                 quant_ratio.max = 3.0,
                                                 norm_levels = list(msrun = list(cond_col="msrun", nmschan_ratio.min=1.0, max_objs=1000L),
                                                                    condition = list(cond_col="condition", max_objs=1000L)))
                                           

msdata$msrun_shifts <- dplyr::semi_join(msruns_hnorm$mschannel_shifts, dplyr::select(msdata$msruns, msrun))

rmsglmdata_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_data_', data_info$fit_ver, '.RData'))
message('Saving data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata, msglm_def, msruns_hnorm,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', data_info$data_ver, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full, fasta.dfs,
     file = rfulldata_filepath)

message('Done.')

