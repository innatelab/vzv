# Virginie Girault VZV project - APMS data 
# Reshape the Maxquant output 
# Define the linear model, effect, contrasts, priors
# Calculate normalization shifts per sample 

project_id <- 'vgirault_vzvapms'
message('Project ID=', project_id)
data_version <- "20180301"
mq_folder <- "mq20180301"
message('Dataset version is ', data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(misc_scripts_path, 'reshape.R'))
source(file.path(misc_scripts_path, 'ggplot_ext.R'))
source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, 'matrix_utils.R'))

require(msglm)
require(msimportr)
require(dplyr)
require(rstan)
require(rjson)
require(stringr)
require(gdata)
require(multidplyr)
require(gplots)

data_path <- file.path(base_data_path, project_id)
mqdata_path <- file.path(data_path, mq_folder)

# import fasta files used for quanting
v5_fasta.df <- read_innate_uniprot_fasta(file.path(data_path, "fasta/V5.fasta"))
gfp_fasta.df <- read_innate_uniprot_fasta(file.path(data_path, "fasta/GFP.fasta"))
vzv_fasta.df <- read_innate_uniprot_fasta(file.path(data_path, "fasta/VZV_vOka_20171208.fasta"))
uphuman_fasta.df <- dplyr::bind_rows(dplyr::mutate(read_innate_uniprot_fasta(file.path(base_data_path, "mqfasta/UP000005640_9606.fasta")), protset = "primary"),
                                     dplyr::mutate(read_innate_uniprot_fasta(file.path(base_data_path, "mqfasta/UP000005640_9606_additional.fasta")), protset = "additional"))

# import calibration file 
protgroup_instr_calib_json_filename <- "instr_qep5_protgroup_LFQ_calib_vgirault_vzvapms_20180103_borg.json"
message( 'Loading MS instrument calibration data from ', protgroup_instr_calib_json_filename, '...' )
protgroup_instr_calib <- fromJSON(file = file.path(data_path, protgroup_instr_calib_json_filename))$instr_calib

# rawfile (msrun)  annotation
msruns.df <- read_tsv(file.path(data_path, mq_folder, "msruns_info.txt"), col_names = TRUE) %>%
    dplyr::mutate(is_after_cleaning = is_after_cleaning == 1,
                  is_used = sample_type != "proteome_fraction",
                  plate_row = factor(plate_row, ordered = TRUE),
                  ms_date = as.POSIXct(strptime(ms_date, "%m/%d/%Y"))) %>%
    dplyr::rename(msrun_mq = msrun, msrun = msrun_fixed)
rawfiles.df <- read_tsv(file.path(data_path, mq_folder, "rawfiles_info.txt"), col_names = TRUE) %>%
    dplyr::mutate(file_date = as.POSIXct(strptime(file_date, "%m/%d/%Y %H:%M %p")),
                  raw_file = str_replace(raw_file, "\\.raw$", "")) %>%
    dplyr::rename(ms_date = file_date)
msruns.df <- dplyr::left_join(dplyr::select(msruns.df, -ms_date), rawfiles.df)

conditions.df <- dplyr::select(dplyr::filter(msruns.df, sample_type != "proteome_fraction"), bait_orf, sample_type) %>%
    dplyr::distinct() %>%
    dplyr::mutate(condition = if_else(sample_type != "proteome", bait_orf, "SKNproteome"),
                  sample_type = factor(sample_type, levels=c("proteome", "APMS_NEG", "APMS_POS", "APMS")),
                  condition = factor(condition))

baits.df <- dplyr::select(dplyr::filter(msruns.df, str_detect(sample_type, "APMS")), bait_orf, sample_type) %>%
    dplyr::distinct() %>%
    dplyr::mutate(bait_type = case_when(.$sample_type == "APMS" ~ "SA",
                                        .$sample_type == "APMS_POS" ~ "POS",
                                        .$sample_type == "APMS_NEG" ~ "NEG") %>% factor(levels=c("NEG", "POS", "SA"))) %>%
    dplyr::mutate(bait_orf_ix = as.integer(str_replace(bait_orf, "[ac]$", ""))) %>%
    dplyr::arrange(bait_type, bait_orf_ix, bait_orf) %>%
    dplyr::mutate(bait_orf = factor(bait_orf, levels=unique(bait_orf))) %>%
    dplyr::mutate(gene_name = if_else(bait_type == "SA", paste0("ORF", bait_orf), NA_character_)) %>%
    dplyr::mutate(gene_name = case_when(.$bait_orf == "335" ~ "ORF33",
                                        .$bait_orf == "101c" ~ "GFP",
                                        .$bait_orf == "103c" ~ "AIFM1",
                                        .$bait_orf == "60c" ~ "ORF60",
                                        .$bait_orf == "66c" ~ "ORF66",
                                        .$bait_orf == "9a" ~ "ORF9",
                                        TRUE ~ .$gene_name)) %>% # FIX
    dplyr::left_join(dplyr::filter(dplyr::select(bind_rows(vzv_fasta.df, dplyr::filter(uphuman_fasta.df, protset=="primary"), gfp_fasta.df),
                                                 gene_name, protein_ac, protein_name, protein_description),
                                   !is.na(gene_name))) %>%
    dplyr::select(-sample_type)

conditions.df <- dplyr::left_join(dplyr::select(conditions.df, -bait_orf),
                                  dplyr::mutate(baits.df, condition=factor(bait_orf, levels=levels(conditions.df$condition)))) %>%
    dplyr::arrange(sample_type, bait_orf_ix, bait_orf)

# Read Maxquant output
quant_type <- "LFQ"
quant_col_prefix <- "LFQ_Intensity"
ms_data.wide <- read.MaxQuant.ProteinGroups(file.path(mqdata_path, 'combined/txt'), import_data = c(quant_type, "ident_type"))
ms_data_colgroups <- attr(ms_data.wide, "column_groups")

ms_data <- list()

# fix protein ac of viral baits
fix_virus_protgroups <- function(protgroups.df) {
  vir_mask <- !is.na(protgroups.df$fasta_headers) & str_detect(protgroups.df$fasta_headers, "OS=Varicella-zoster")
  protgroups.df$is_viral <- vir_mask
  if (!any(vir_mask)) return(protgroups.df)
  vir_fasta_chunks <- str_match(protgroups.df$fasta_headers[vir_mask],
                                "([\\w_]+)\\s(.+) OS=Varicella-zoster virus \\(strain Oka vaccine\\) GN=(\\w+)")
  print(vir_fasta_chunks)
  protgroups.df[vir_mask, 'protein_names'] <- vir_fasta_chunks[,3]
  protgroups.df[vir_mask, 'gene_names'] <- vir_fasta_chunks[,4]
  protgroups.df$is_viral <- vir_mask
  return(protgroups.df)
}

# dataset reshaping
ms_data$protgroups = dplyr::select_(ms_data.wide, .dots=ms_data_colgroups$protgroup) %>%
  dplyr::distinct() %>% dplyr::arrange(protgroup_id) %>%
  dplyr::mutate(protein_names = if_else(is.na(protein_names), as.character(fasta_headers), protein_names))

ms_data$protein2protgroups <- dplyr::bind_rows(
    dplyr::mutate(expand_protgroup_acs(ms_data$protgroups, acs_col="majority_protein_acs"),
                  is_majority = TRUE) %>%
    dplyr::rename(protein_ac = majority_protein_ac),
    dplyr::mutate(expand_protgroup_acs(ms_data$protgroups, acs_col="protein_acs"),
                  is_majority = FALSE)) %>%
    dplyr::group_by(protgroup_id, protein_ac) %>%
    dplyr::summarise(is_majority = any(is_majority)) %>%
    dplyr::ungroup()

ms_data$proteins <- dplyr::bind_rows(dplyr::mutate(vzv_fasta.df, is_viral = TRUE),
                                     dplyr::mutate(uphuman_fasta.df, is_viral = FALSE),
                                     dplyr::mutate(v5_fasta.df, is_viral = FALSE),
                                     dplyr::mutate(gfp_fasta.df, is_viral = FALSE)) %>%
    dplyr::semi_join(ms_data$protein2protgroups) %>%
    dplyr::left_join(dplyr::select(dplyr::filter(ms_data$protein2protgroups, is_majority), protein_ac, protgroup_id)) %>%
  dplyr::mutate(protein_existence = as.integer(protein_existence),
                protein_ac_noiso = stripUniprotIsoform(protein_ac)) %>%
  dplyr::group_by(protein_ac_noiso) %>%
  dplyr::mutate(protein_existence = if_else(is.na(protein_existence), as.integer(min(protein_existence, na.rm=TRUE)), protein_existence)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-protein_ac_noiso) %>%
  dplyr::mutate(ac_penalty = if_else(!is.na(gene_name), 0L, 10L) +
                             if_else(!str_detect(protein_ac, "-\\d+$"), 0L, 1L) +
                             if_else(src_db == "sp", 0L, 3L) +
                             if_else(!is.na(protein_existence), protein_existence - 1L, 6L))

# fix proein ac 
fix_protgroups.df <- dplyr::left_join(dplyr::select(ms_data$protgroups, protgroup_id), ms_data$proteins) %>%
    dplyr::mutate(alt_protein_ac_noiso = str_match(protein_description, "^Isoform of ([^ ,]+)")[,2],
                  protein_ac_noiso = stripUniprotIsoform(protein_ac),
                  has_iso = protein_ac != protein_ac_noiso) %>%
    dplyr::arrange(protgroup_id, ac_penalty, gene_name, protein_ac) %>%
    dplyr::group_by(protgroup_id) %>%
    dplyr::mutate(has_canonical = (protein_ac_noiso %in% unique(protein_ac)) | (alt_protein_ac_noiso %in% unique(protein_ac)),
                  is_canonical = (protein_ac == protein_ac_noiso) & is.na(alt_protein_ac_noiso),
                  use4prot_name = !any(is_canonical | !has_canonical, na.rm = TRUE) | is_canonical | !has_canonical,
                  use4gene_name = !any(!is.na(gene_name)) | !is.na(gene_name)) %>%
    dplyr::summarise(protein_names = str_c(unique(protein_description[use4prot_name]), collapse=";"),
                     gene_names = str_c(unique(gene_name[use4gene_name]), collapse=';'),
                     is_viral = any(is_viral)) %>%
    dplyr::ungroup()

all(fix_protgroups.df$protgroup_id == ms_data$protgroups$protgroup_id)

ms_data$protgroups$protein_names[is.na(ms_data$protgroups$protein_names)] <-
    fix_protgroups.df$protein_names[is.na(ms_data$protgroups$protein_names)]
ms_data$protgroups$gene_names[is.na(ms_data$protgroups$gene_names)] <-
    fix_protgroups.df$gene_names[is.na(ms_data$protgroups$gene_names)]
ms_data$protgroups$is_viral <- fix_protgroups.df$is_viral
ms_data$protgroups$is_contaminant[!is.na(ms_data$protgroups$gene_names) & str_detect(ms_data$protgroups$gene_names, '^(GFP|V5tag)$')] <- FALSE

ms_data$mschannels <- dplyr::filter(msruns.df, is_used) %>%
  dplyr::mutate(quant_type = quant_type,
                bait_orf = factor(bait_orf, levels=levels(conditions.df$bait_orf)),
                bio_batch = factor(bio_batch),
                sample_type = factor(sample_type, levels=levels(conditions.df$sample_type))) %>%
  dplyr::left_join(conditions.df) %>%
  dplyr::arrange(bait_orf, condition, bio_replicate, tech_replicate, msrun) %>%
  dplyr::mutate(condXbiobatch = paste0(condition,'_',bio_batch),
                condXbiobatch = factor(condXbiobatch),
                msrun_ix = row_number(),
                msrun = factor(msrun))

protgroup_intensities_all.df <- reshape(ms_data.wide[, c('protgroup_id', ms_data_colgroups$LFQ)],
                                        direction="long", idvar=c("protgroup_id"),
                                        timevar = "msrun_mq", sep=" ", v.names = paste0(quant_col_prefix, ".Sum")) %>%
  dplyr::rename_(.dots = c(intensity = paste0("`", paste0(quant_col_prefix, ".Sum"),"`")))
ms_data$protgroup_intensities <- dplyr::filter(protgroup_intensities_all.df, !is.na(intensity)) %>%
  dplyr::inner_join(dplyr::select(ms_data$mschannels, msrun, msrun_mq)) %>%
  dplyr::select(-msrun_mq)

ms_data$mschannels <- dplyr::mutate(ms_data$mschannels,
                                    has_data = msrun %in% unique(ms_data$protgroup_intensities$msrun),
                                    mstag = factor("Sum", levels = c("L", "M", "H", "Sum")))

ms_data$protgroup_intensities_wide <- ms_data.wide[, c('protgroup_id', ms_data_colgroups$LFQ)]
ms_data$protgroup_idents_wide <- ms_data.wide[, c('protgroup_id', ms_data_colgroups$ident_type)]
ms_data$protgroup_idents <- reshape(ms_data$protgroup_idents_wide, idvar = 'protgroup_id',
                                    direction = "long", timevar="msrun", v.names = "ident_type") %>%
    dplyr::filter(!is.na(ident_type))

# keep only proteins identified in APMS
ms_data$protgroups <- semi_join(ms_data$protgroups,
                                dplyr::select(filter(semi_join(ms_data$protgroup_idents,
                                                               filter(ms_data$mschannels, str_detect(sample_type, "APMS"))),
                                                     ident_type == "By MS/MS"), protgroup_id))
ms_data$protgroup_idents_wide <- dplyr::semi_join(ms_data$protgroup_idents_wide, ms_data$protgroups)
ms_data$protgroup_idents <- dplyr::semi_join(ms_data$protgroup_idents, ms_data$protgroups)
ms_data$protgroup_intensities_wide <- dplyr::semi_join(ms_data$protgroup_intensities_wide, ms_data$protgroups)
ms_data$protgroup_intensities <- dplyr::semi_join(ms_data$protgroup_intensities, ms_data$protgroups)

bait_checks.df <- dplyr::left_join(baits.df, dplyr::select(ms_data$proteins, protein_ac, protgroup_id)) %>%
    dplyr::left_join(dplyr::select(dplyr::filter(ms_data$protgroup_idents, ident_type=="By MS/MS"), protgroup_id, msrun)) %>%
    dplyr::left_join(dplyr::select(ms_data$mschannels, msrun, observed_bait = bait_orf, condition)) %>%
    dplyr::arrange(bait_orf, protgroup_id, msrun) %>%
    dplyr::group_by(bait_orf, protgroup_id) %>%
    dplyr::mutate(msruns = paste0(unique(msrun), collapse=";"),
                  observed_baits = paste0(unique(observed_bait), collapse=";"),
                  conditions = paste0(unique(condition), collapse=";")) %>%
    dplyr::filter(row_number()==1L) %>%
    dplyr::select(-msrun, -observed_bait, -condition) %>%
    dplyr::ungroup()


ms_data$protgroup_stats <- dplyr::semi_join(protgroup_intensities_all.df, ms_data$protgroups) %>%
    dplyr::inner_join(ms_data$mschannels) %>%
    dplyr::filter(mstag == 'Sum') %>%
    dplyr::group_by(protgroup_id) %>%
    summarize(nmsruns = n(),
              nconditions = n_distinct(condition),
              nmsrun_quanted = sum(!is.na(intensity)),
              ncondition_quanted = n_distinct(condition[!is.na(intensity)]),
              median_intensity = median(intensity, na.rm=TRUE)) %>%
    dplyr::ungroup()

ms_data$protgroups <- dplyr::mutate(ms_data$protgroups,
    is_full_quant = protgroup_id %in% (dplyr::filter(ms_data$protgroup_stats, ncondition_quanted==nconditions) %>% .$protgroup_id),
    is_top_quant = protgroup_id %in% (dplyr::filter(ms_data$protgroup_stats, percent_rank(-median_intensity) <= 0.1) %>% .$protgroup_id))

# setup experimental design matrices (APMS + APMS:Baiti)
conditionXeffect.mtx <- model.matrix(~ 1 + is_APMS + bait_orf,
                                     dplyr::mutate(conditions.df,
                                                   is_APMS = str_detect(sample_type, "APMS"),
                                                   bait_orf = relevel(factor(if_else(sample_type == "proteome", "none", as.character(bait_orf))),
                                                                      "none")))

# intercept handled separately, msMar could not be decoupled from bioMar
conditionXeffect.mtx <- conditionXeffect.mtx[,setdiff(colnames(conditionXeffect.mtx),c('(Intercept)'))]
dimnames(conditionXeffect.mtx) <- list(condition = conditions.df$condition,
                                       effect = colnames(conditionXeffect.mtx))

effects.df <- data.frame(effect=colnames(conditionXeffect.mtx),
                         stringsAsFactors=FALSE) %>%
  dplyr::mutate(bait_orf = effect_factor(effect, "bait_orf", levels(baits.df$bait_orf), NA),
                #batch = effect_factor(effect, "batch", levels(conditions.df$batch), NA),
                is_positive = FALSE) %>%
  dplyr::mutate(tau = case_when(.$effect == "is_APMSTRUE" ~ 5.0,
                                !is.na(.$bait_orf) ~ 1.0,
                                TRUE ~ 1.0))
effects.df$effect_label <- sapply(1:nrow(effects.df), function(i) {
  comps <- c()
  if (effects.df$effect[[i]] == "is_APMSTRUE") comps <- append(comps, "APMS")
  if (!is.na(effects.df$bait_orf[[i]]))  comps <- append(comps, paste0("bait_orf(", effects.df$bait_orf[[i]], ")"))
  paste0(comps, collapse="+")
})

conditionXeffect.df <- conditionXeffect_frame(conditionXeffect.mtx, effects.df)
inv_conditionXeffect.mtx <- frame2matrix(conditionXeffect.df,
                                         "condition", "effect", "cond_w",
                                         rows = rownames(conditionXeffect.mtx),
                                         cols = colnames(conditionXeffect.mtx))
heatmap.2(inv_conditionXeffect.mtx, margins=c(16, 16), Rowv=FALSE, Colv=FALSE, dendrogram = "none")


msrunXreplEffect.mtx <- replicate_effects_matrix(
  dplyr::group_by(ms_data$mschannels, condition, bio_batch) %>%
  dplyr::mutate(bio_replicate = if_else(!is.na(bio_replicate), bio_replicate, 1L)) %>%
  dplyr::ungroup(),
  replicate_col = "bio_replicate", cond_col = "condition")

msrunXreplEffect.df <- as.data.frame(as.table(msrunXreplEffect.mtx)) %>%
  dplyr::filter(Freq != 0) %>% dplyr::select(-Freq)

#Metaconditions required to design contrasts
compound_metaconditions <- paste0("allMinus", levels(conditions.df$bait_orf))
all_metaconditions <- c(levels(conditions.df$condition), compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in levels(conditions.df$condition)) {
    conditionXmetacondition.mtx[cname, cname] <- TRUE
}
for (cname in levels(conditions.df$bait_orf)) {
    allMinusBait <- paste0("allMinus", cname)
    conditionXmetacondition.mtx[, allMinusBait] <- TRUE
    conditionXmetacondition.mtx[c("SKNproteome", cname), allMinusBait] <- FALSE
}
conditionXmetacondition.mtx['9a', 'allMinus9'] <- FALSE
conditionXmetacondition.mtx['9', 'allMinus9a'] <- FALSE
conditionXmetacondition.mtx['66', 'allMinus66c'] <- FALSE
conditionXmetacondition.mtx['66c', 'allMinus66'] <- FALSE
conditionXmetacondition.mtx['60', 'allMinus60c'] <- FALSE
conditionXmetacondition.mtx['60c', 'allMinus60'] <- FALSE


conditionXmetacondition.df <- as.data.frame(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(Freq) %>% dplyr::select(-Freq)

#Contrasts 
bait_vs_gfp_l <- as.character(dplyr::filter(baits.df, bait_type != "NEG")$bait_orf)
bait_vs_gfp_r <- rep(as.character(dplyr::filter(baits.df, bait_type == "NEG")$bait_orf), length(bait_vs_gfp_l))
bait_vs_others_l <- levels(baits.df$bait_orf)
bait_vs_others_r <- paste0("allMinus", bait_vs_others_l)
contrast_l <- c(bait_vs_gfp_l, bait_vs_others_l)
contrast_r <- c(bait_vs_gfp_r, bait_vs_others_r)
all_contrasts <- paste0(contrast_l, "_vs_", contrast_r)

contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
for (i in 1:length(all_contrasts)) {
    contrastXmetacondition.mtx[all_contrasts[[i]], c(contrast_l[[i]], contrast_r[[i]])] <- c(1, -1)
}

contrastXmetacondition.df <- as.data.frame(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(Freq != 0) %>%
  dplyr::rename(weight = Freq) %>%
  dplyr::mutate(contrast_type = 'filter',
                condition_role = if_else(weight > 0, "signal", "background"))

contrastXcondition.df <- as.data.frame(as.table(conditionXmetacondition.mtx)) %>% dplyr::filter(Freq != 0) %>%
  dplyr::select(-Freq) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition)

ms_data$mschannel_stats <- dplyr::left_join(tidyr::expand(dplyr::filter(ms_data$protgroup_intensities, !is.na(protgroup_id)), protgroup_id, msrun),
                                            dplyr::filter(ms_data$protgroup_intensities, !is.na(protgroup_id))) %>%
  dplyr::left_join(ms_data$protgroup_idents) %>%
  dplyr::inner_join(dplyr::select(ms_data$mschannels, msrun, condition) %>% dplyr::distinct()) %>%
  dplyr::group_by(protgroup_id, condition) %>% dplyr::filter(any(!is.na(intensity))) %>%
  dplyr::group_by(msrun) %>%
  summarize(log_intensity.F_mean = mean(log(intensity[!is.na(intensity)])),
            log_intensity.F_sd = sd(log(intensity[!is.na(intensity)])),
            n = n(),
            n_imputed = sum(is.na(intensity)),
            n_matching = sum(ident_type=="By matching", na.rm = TRUE),
            n_msms = sum(ident_type=="By MS/MS", na.rm = TRUE)) %>%
  dplyr::ungroup()

protgroup_intensities.df <- dplyr::inner_join(protgroup_intensities_all.df,
                                              dplyr::select(ms_data$mschannels, msrun, msrun_mq)) %>%
  dplyr::inner_join(ms_data$mschannel_stats) %>%
  dplyr::group_by(protgroup_id) %>%
  dplyr::mutate(has_missing_intensity = any(is.na(intensity))) %>%
  dplyr::group_by(msrun) %>%
  dplyr::mutate(intensity_imputed.F = if_else(is.na(intensity),
                                              exp(rnorm(n(), mean=log_intensity.F_mean-1.8, sd=log_intensity.F_sd*0.3)),
                                              intensity)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(#intensity_norm.F = exp(-total_msrun_shift)*intensity,
                #intensity_imputed_norm.F = exp(-total_msrun_shift)*intensity_imputed.F,
                msrun = factor(msrun))

bait_intensities.df <- dplyr::select(protgroup_intensities.df, protgroup_id,
                                     has_missing_intensity, msrun, intensity.F=intensity, intensity_imputed.F#, intensity_norm.F, intensity_imputed_norm.F
                                     ) %>%
                       dplyr::inner_join(dplyr::select(dplyr::filter(ms_data$protgroups, is_viral | protgroup_id %in% bait_checks.df$protgroup_id), protgroup_id, majority_protein_acs)) %>%
                       dplyr::inner_join(dplyr::select(msruns.df, msrun, ms_bait_orf=bait_orf, bio_replicate, tech_replicate)) %>%
                       dplyr::left_join(dplyr::select(dplyr::inner_join(dplyr::select(baits.df, bait_orf, protein_ac), ms_data$proteins), protgroup_id, bait_orf)) %>%
                       dplyr::left_join(ms_data$protgroup_idents) %>%
                       dplyr::mutate(bait_ix = as.integer(bait_orf),
                                     ms_bait_ix = as.integer(ms_bait_orf)) %>%
                       dplyr::arrange(bait_ix, ms_bait_ix, bio_replicate, tech_replicate)

bait_intensities_wide.df <- reshape(dplyr::mutate(bait_intensities.df,
                                                  intensity = ifelse(ident_type == "By MS/MS", intensity.F, NA_real_)) %>%
                                    dplyr::select(msrun, protgroup_id, bait_orf, has_missing_intensity, intensity),
                                    direction = "wide", idvar = c("protgroup_id", "has_missing_intensity", "bait_orf"),
                                    timevar = "msrun", v.names = "intensity")
bait_intensities_wide.mtx <- dplyr::select(bait_intensities_wide.df,
                                           -has_missing_intensity, -protgroup_id, -bait_orf) %>% log %>% as.matrix
colnames(bait_intensities_wide.mtx) <- str_replace(colnames(bait_intensities_wide.mtx), "intensity.", "")
rownames(bait_intensities_wide.mtx) <- bait_intensities_wide.df$bait_orf
bait_intensities_wide_imp.df <- reshape(dplyr::select(bait_intensities.df, msrun, bait_orf, protgroup_id, intensity = intensity_imputed.F),
                                    direction = "wide", idvar = c("protgroup_id", "bait_orf"),
                                    timevar = "msrun", v.names = "intensity")
bait_intensities_wide_imp.mtx <- bait_intensities_wide_imp.df %>% #dplyr::semi_join(dplyr::filter(bait_intensities_wide_imp.df, TRUE),
                                                  #dplyr::filter(ms_data$protgroups, !is_contaminant & !is_reverse)) %>%
                                 dplyr::select(-protgroup_id, -bait_orf) %>% log %>% as.matrix
colnames(bait_intensities_wide_imp.mtx) <- str_replace(colnames(bait_intensities_wide_imp.mtx), "intensity_imputed_norm.F.", "")
bait_intensities_wide_imp.dist <- dist(as.matrix(bait_intensities_wide_imp.mtx))
bait_intensities_wide_imp.tdist <- dist(t(as.matrix(bait_intensities_wide_imp.mtx)))
pdf(file = file.path(data_path, "plots", paste0("bait_heatmap_MSMS_Intensity_2_h.pdf")), width=20, height=50)
heatmap.2(t(as.matrix(bait_intensities_wide.mtx)), Colv = NULL, Rowv = NULL,
          trace = "none", #Colv = NULL,
          na.color = "gray", key=FALSE, cexRow = 1.5, cexCol = 1.5,
          lmat = matrix(c(0,2,0,3,1,0), ncol=2), lhei = c(0.003, 0.95, 0.01), lwid = c(0.1, 0.9))
dev.off()

protgroup_intensities_wide.df <- reshape(dplyr::select(protgroup_intensities.df, protgroup_id,
                                                       msrun, intensity = intensity) %>%
                                         dplyr::semi_join(filter(ms_data$mschannels, str_detect(sample_type, "APMS"))) %>%
                                         dplyr::semi_join(ms_data$protgroups),
                                         direction = "wide", idvar = c("protgroup_id"),
                                         timevar = "msrun", v.names = "intensity")
protgroup_intensities_wide.mtx <- dplyr::semi_join(dplyr::filter(protgroup_intensities_wide.df, TRUE),
                                                  dplyr::filter(ms_data$protgroups, !is_contaminant & !is_reverse)) %>%
                                 dplyr::select(-protgroup_id) %>% log %>% as.matrix
colnames(protgroup_intensities_wide.mtx) <- str_replace(colnames(protgroup_intensities_wide.mtx), "intensity.", "")
protgroup_intensities_imp_wide.df <- reshape(dplyr::select(protgroup_intensities.df, protgroup_id,
                                                           msrun, intensity = intensity_imputed.F) %>%
                                             dplyr::semi_join(filter(ms_data$mschannels, str_detect(sample_type, "APMS"))) %>%
                                             dplyr::semi_join(ms_data$protgroups),
                                         direction = "wide", idvar = c("protgroup_id"),
                                         timevar = "msrun", v.names = "intensity") %>%
    dplyr::left_join(dplyr::select(dplyr::mutate(ms_data$protgroups, is_used = !is_contaminant & !is_reverse),
                                   protgroup_id, is_used))
protgroup_intensities_imp_wide.mtx <- dplyr::filter(protgroup_intensities_imp_wide.df, is_used) %>%
                                 dplyr::select(-protgroup_id, -is_used) %>% log %>% as.matrix
colnames(protgroup_intensities_imp_wide.mtx) <- str_replace(colnames(protgroup_intensities_imp_wide.mtx), "intensity.", "")
protgroup_intensities_imp_wide.dist <- dist(as.matrix(protgroup_intensities_imp_wide.mtx))
protgroup_intensities_imp_wide.tdist <- dist(t(as.matrix(protgroup_intensities_imp_wide.mtx)))
pdf(file = file.path(data_path, "plots", paste0("heatmap_Intensity_h.pdf")), width=30, height=100)
heatmap.2(as.matrix(protgroup_intensities_wide.mtx),
          distfun = function(mtx) if (nrow(mtx) == nrow(protgroup_intensities_wide.mtx)) protgroup_intensities_imp_wide.dist else protgroup_intensities_imp_wide.tdist,
          trace = "none", #Colv = NULL,
          na.color = "gray", key=FALSE, cexRow = 0.2,
          lmat = matrix(c(0,2,0,3,1,0), ncol=2), lhei = c(0.003, 0.95, 0.01), lwid = c(0.1, 0.9))
dev.off()

require(FactoMineR)
protgroup_intensities_pca <- PCA(protgroup_intensities_imp_wide.mtx, graph = FALSE)
protgroup_intensities_pca.df <- as.data.frame(protgroup_intensities_pca$svd$V)
colnames(protgroup_intensities_pca.df) <- paste0("comp_", 1:ncol(protgroup_intensities_pca.df))
protgroup_intensities_pca.df <- dplyr::mutate(protgroup_intensities_pca.df,
                                              msrun = str_replace(rownames(protgroup_intensities_pca$var$coord), "intensity.", "")) %>%
  dplyr::inner_join(dplyr::select(ms_data$mschannels, msrun, ms_date, lc_column, bait_orf, tech_replicate, bio_replicate, bait_type, bio_batch, lc_column
                                  #, plate_column, plate_row, file_size
                                  )) %>%
  dplyr::left_join(ms_data$mschannel_stats)

protgroup_intensities_pca_wider.df <- reshape(protgroup_intensities_pca.df, direction="wide", idvar=c("bait_orf", "bio_replicate"),
                                              timevar = "tech_replicate",
                                              v.names = c("comp_1", "comp_2", "comp_3"
                                                          #, "file_size"
                                                          ))

protgroup_intensities_pca.df <- dplyr::mutate(protgroup_intensities_pca.df,
                                              contamination = ifelse(comp_2 > median(comp_2) + mad(comp_2), comp_2 - median(comp_2), 0.0),
                                              contamination = contamination/max(contamination))

protgroup_intensities_imp_wide.df[protgroup_intensities_imp_wide.df$is_used, "comp_2"] <- protgroup_intensities_pca$svd$U[, 2]
protgroup_intensities_imp_wide.df <- dplyr::mutate(protgroup_intensities_imp_wide.df,
                                                   is_comp2 = if_else(!is.na(comp_2), comp_2 > median(comp_2, na.rm=TRUE) + mad(comp_2,na.rm=TRUE), FALSE))

ms_data$protgroups <- dplyr::left_join(ms_data$protgroups, dplyr::select(protgroup_intensities_imp_wide.df, protgroup_id, is_comp2))

options(mc.cores=8)

source(file.path(base_scripts_path, "msglm/R/stan_models.R"))
source(file.path(base_scripts_path, "msglm/R/normalize_utils.R"))

ms_data4norm.df <- ms_data$protgroup_intensities %>%
  dplyr::semi_join(dplyr::filter(ms_data$protgroups, !is_reverse & !is_contaminant & !is_viral & !is_comp2 & is_full_quant))

# normalize experiments:
mschannel_hnorm <- multilevel_normalize_experiments(msglm_normalize.stan_model, protgroup_instr_calib,
    dplyr::select(ms_data$mschannels, msrun, sample_type, condition, bio_batch, bio_replicate, tech_replicate) %>%
        dplyr::mutate(condXbiobatch = paste0(condition, '_', bio_batch),
                      sampletype_biobatch = paste0(sample_type, '_', bio_batch)),
    ms_data4norm.df,
    quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
    mcmc.iter = 1000L,
    #mcmc.chains = 6,
    verbose=TRUE,
    norm_levels = list(msrun = list(cond_col = "msrun", max_objs=500L, missing_exp.ratio=0.1),
                       condXbiobatch = list(cond_col = "condXbiobatch", max_objs=200L, missing_exp.ratio=0.1),
                       biobatch = list(cond_col = "sampletype_biobatch", max_objs = 100L),
                       sample_type = list(cond_col = "sample_type", max_objs = 50L)
                       ))
total_msrun_shifts.df <- mschannel_hnorm$mschannel_shifts

global_protgroup_labu_shift <- 0.95*median(log(ms_data$protgroup_intensities$intensity), na.rm=TRUE)

# setup batch effect
msruns.df <- dplyr::left_join(msruns.df, dplyr::select(dplyr::mutate(protgroup_intensities_pca.df, is_contaminated = contamination > 0),
                                                       msrun, is_contaminated, contamination))
ms_data$mschannels <- dplyr::left_join(ms_data$mschannels, dplyr::select(dplyr::mutate(protgroup_intensities_pca.df, is_contaminated = contamination > 0),
                                                                         msrun, is_contaminated, contamination))

#Batch effect: BioChembatchi+MSBatchi+PCAbatch
msrunXbatchEffect.mtx <- model.matrix(~ 1 + bio_batch + ms_batch + contamination,
                                      mutate(ms_data$mschannels,
                                             # 3rd msbatch is the same as the 2nd biobatch
                                             ms_batch = factor(ms_batch),
                                             contamination = if_else(is.na(contamination), 0, contamination)))

# intercept handled separately, msMar could not be decoupled from bioMar
msrunXbatchEffect.mtx <- msrunXbatchEffect.mtx[,setdiff(colnames(msrunXbatchEffect.mtx),c('(Intercept)'))]
dimnames(msrunXbatchEffect.mtx) <- list(msrun = ms_data$mschannels$msrun,
                                        batch_effect = colnames(msrunXbatchEffect.mtx))
pdf(file = file.path(data_path, "plots", paste0("msrunXbatchEffect_", data_version, ".pdf")), width=6, height=20)
heatmap.2(msrunXbatchEffect.mtx, key=FALSE, cexRow = 0.5, cexCol = 1.5,
          lmat = matrix(c(0,2,0,3,1,0), ncol=2), lhei = c(0.003, 0.9, 0.02), lwid = c(0.1, 0.9), trace = "none",
          Rowv=FALSE, Colv=FALSE, dendrogram = "none")
dev.off()

batch_effects.df <- data.frame(batch_effect=colnames(msrunXbatchEffect.mtx),
                               stringsAsFactors=FALSE) %>%
  dplyr::mutate(bio_batch = effect_factor(batch_effect, "bio_batch", levels(ms_data$mschannels$bio_batch), NA),
                ms_batch = effect_factor(batch_effect, "ms_batch", as.character(unique(ms_data$mschannels$ms_batch)), NA),
                is_positive = batch_effect == "contamination",
                batch_effect_label = batch_effect) %>%
  dplyr::mutate(tau = 2.0)

def_norm_data <- protgroup_instr_calib[c('zDetectionFactor', 'zDetectionIntercept',
                                         'detectionMax', 'sigmaScaleHi', 'sigmaScaleLo',
                                         'sigmaOffset', 'sigmaBend', 'sigmaSmooth',
                                         'zShift', 'zScale')]

mqevidence <- read.MaxQuant.Evidence(file.path(data_path, mq_folder, 'combined/txt'), file_name = "evidence_fixed.txt",
                                     evidence.pepobj = "pepmodstate", correct_ratios = FALSE)#, correct_by_ratio.ref_label='M', mode="labeled")

rmq_filepath <- file.path(scratch_path, paste0(project_id, '_mqevidence_', data_version, '.RData'))
message('Saving MQ evidence to ', rmq_filepath, '...')
data_info <- list(project_id = project_id, version = data_version, mq_folder = mq_folder)
save(data_info, mqevidence, file = rmq_filepath)

msdata_full <- mqevidence
msdata_full$protgroups <- ms_data$protgroups
msdata_full$protein2protgroups <- ms_data$protein2protgroups
msdata_full$msruns <- filter(msdata_full$mschannels, mstag == "F") %>%
  dplyr::select(-mstag)
msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  dplyr::inner_join(select(mqevidence$pepmodstates, pepmodstate_id, charge)) %>%
  dplyr::rename_at(vars(one_of("MSMS", "MULTI-MSMS", "MULTI-SECPEP", "MULTI-MATCH")), ~str_c("ident_type.", .)) %>%
  dplyr::select(pepmodstate_id, pepmod_id, msrun, charge,
                starts_with("intensity"), starts_with("ident_type")) %>%
  dplyr::rename_at(vars(matches("^intensity.*\\.F$")), ~str_remove(., "\\.F$")) %>%
  filter(!is.na(intensity)) %>% # FIXME check all intensity channels
  dplyr::mutate(ident_type = factor(case_when(
    #`ident_type.ISO-MSMS` ~ "ISO-MSMS",
    `ident_type.MULTI-MSMS` ~ "MULTI-MSMS",
    `ident_type.MSMS` ~ "MSMS",
    `ident_type.MULTI-MATCH` ~ "MULTI-MATCH",
    `ident_type.MULTI-SECPEP` ~ "MULTI-SECPEP",
    TRUE ~ NA_character_),
    levels = c("ISO-MSMS", "MULTI-MSMS", "MSMS", "MULTI-SECPEP", "MULTI-MATCH"))) %>%
  dplyr::left_join(dplyr::select(msdata_full$msruns, msrun)) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift))


rmsfull_filepath <- file.path(scratch_path, paste0(project_id, '_msfull_', data_version, '.RData'))
message('Saving evidence-based MS data to ', rmsfull_filepath, '...')
data_info <- list(project_id = project_id, version = data_version, mq_folder = mq_folder)
save(data_info, msdata_full, file = rmsfull_filepath)

rdata_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_data_', data_version, '.RData'))
message( 'Saving MS data to ', rdata_filepath, '...' )
data_info <- list(project_id = project_id, version = data_version, mq_folder = mq_folder)
save( data_info,
      ms_data,
      baits.df, bait_checks.df, conditions.df, effects.df,
      conditionXeffect.mtx, inv_conditionXeffect.mtx, conditionXeffect.df,
      conditionXmetacondition.mtx, conditionXmetacondition.df,
      contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
      protgroup_instr_calib, global_protgroup_labu_shift,
      mschannel_hnorm, total_msrun_shifts.df,
      msrunXreplEffect.mtx,
      batch_effects.df, msrunXbatchEffect.mtx,
      file = rdata_filepath )
message( 'Done.' )

