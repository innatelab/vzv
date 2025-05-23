Collection of general scripts for Virginie Girault VZV project using in-house packages

list of other code sources required: 
msglm [Bayesian Random Effects Linear Models for Mass Spectrometry Data]
- v0.5.0 https://github.com/innatelab/msglm/tree/v0.5.0 
- v0.6.0 https://github.com/innatelab/msglm/tree/v0.5.0 
msimportr [R package for importing and initial processing of MaxQuant/Spectronaut output] https://github.com/innatelab/msimportr 
HierarchicalHotNet [Julia implementation of Hierarchical HotNet method] HierarchicalHotNet.jl

#Affinity-purification of V5-tagged VZV proteins in SK-N-BE2 cells (interactomes)
requires: msglm package v0.5.0 and msimportr 
1- prepare_data_apms.R
Reads the maxquant output file: ProteinGroups.txt and evidence_fixed.txt
Outputs:
msfull.Rdata containing the formated and annotated dataset, with protein groups and peptides intensities. 
msglm.Rdata containing the formated and annotated dataset with protein groups intensities, sample associated normalisation shifts, mass spectrometer associated calibration factors, 
the GLM model description (effects, batch effects, metaconditions and contrasts) with corresponding priors. 
2 - msglm_fit_chunk_apms.R
Extracts the data of a single protein group based on its protgroup_id [job_chunk]
Fits the Bayesian model on the protein group.
Outputs:
[job_chunk].Rdata containing the msglm results
3 - msglm_fit-chunks_apms.sge.sh
Bash script instructing the fitting of all the protein groups of the dataset [job_chunk] individually.
4 - assemble_fits_apms.R 
Assembles the [job_chunk].Rdata into one. 
Outputs excel files containing summarizing the the results per effect per protein and per contrast per protein.
The results are compiled into the Supplementary Table S3, Tab 1 - VZV-Host interactions and Tab 2 - VZV ORF baits. 

#Full proteome changes induced by the depletion of the MPP8 gene in SK-N-BE2 cells
requires: msglm package v.0.6.0 and msimportr
1- prepare_data_MPP8_KO.R 
Reads the maxquant output files: peptides.txt and evidence.txt 
Protein groups are refined via the protregroup.jl script - 
Protein groups distinguished by only one specific peptide or with less than 25% different specific peptides are merged to extend the set of peptides used for protein group quantitation and reduce the protein isoform-specific changes.
Outputs:
msfull.Rdata containing the formated and annotated dataset, with protein groups and peptides intensities. 
msglm.Rdata containing the formated and annotated dataset with peptide intensities, sample associated normalisation shifts, mass spectrometer associated calibration factors, 
the GLM model description (effects) with corresponding priors. 
2 - msglm_fit_chunk_MPP8_KO.R
Extracts the data of a single protein group based on its protgroup_id [job_chunk]
Fits the Bayesian model on the protein group.
Outputs:
[job_chunk].Rdata containing the msglm results
3 - msglm_fit-chunks_MPP8_KO.lrz.sh
Bash script instructing the fitting of all the protein groups of the dataset [job_chunk] individually.
4 - assemble_fits_MMP8_KO.R 
Assembles the [job_chunk].Rdata into one. 
Outputs an excel file containing summarizing the the results per effect per protein.
These results are integrated in the compiled "Complementary omics datasets" Table S5

#Integration of the interactome and effectome data by network diffusion analysis 
requires:HierarchicalHotNet
1. hotnet_analysis.jl 
General script for the hotnet analysis, dependent on results from 2. and 4.
Maps the Interactome and Effectome hits on ReactomeFI. Tests the significance of the generated edges of the network against permuted data-based networks. 
2.1 hotnet_treestats_chunks.jl
Calculates Strongly Connected Component Tree statistic per viral bait. 
2.2 hotnet_treestats_chunk.lrz.sh 
Bash script 
3.1 hotnet_perm_chunk.jl
Generate 1000 permuted tree per bait and calculate SCC tree statistics. 
3.2 hotnet_perm_chunk.lrz.sh 
Bash script 
4. hotnet_perm_chunk_assemble.jl
Assembles Hotnet permuted tree results. 
