# Collection of general scripts for Virginie Girault VZV project using in-house packages

list of other code sources required:<br> 
- **msglm** Bayesian Random Effects Linear Models for Mass Spectrometry Data<br> 
v0.5.0 https://github.com/innatelab/msglm/tree/v0.5.0 <br> 
v0.6.0 https://github.com/innatelab/msglm/tree/v0.6.0 <br> 
- **msimportr** R package for importing and initial processing of MaxQuant/Spectronaut output https://github.com/innatelab/msimportr <br> 
- **HierarchicalHotNet**[Julia implementation of Hierarchical HotNet method https://github.com/alyst/HierarchicalHotNet.jl/tree/v1.1.0 <br> 

## Affinity-purification of V5-tagged VZV proteins in SK-N-BE2 cells (interactomes)
requires: **msglm** package v0.5.0 and **msimportr**   <br> 
1. **prepare_data_apms.R** <br> 
Reads the maxquant output file: ProteinGroups.txt and evidence_fixed.txt <br> 
Outputs: <br> 
msfull.Rdata containing the formated and annotated dataset, with protein groups and peptides intensities.  <br> 
msglm.Rdata containing the formated and annotated dataset with protein groups intensities, sample associated normalisation shifts, mass spectrometer associated calibration factors, 
the GLM model description (effects, batch effects, metaconditions and contrasts) with corresponding priors.  <br> 
2. **msglm_fit_chunk_apms.R** <br> 
Extracts the data of a single protein group based on its protgroup_id [job_chunk] <br> 
Fits the Bayesian model on the protein group. <br> 
Outputs: <br> 
[job_chunk].Rdata containing the msglm results <br> 
3. **msglm_fit-chunks_apms.sge.sh** <br> 
Bash script instructing the fitting of all the protein groups of the dataset [job_chunk] individually. <br> 
4. **assemble_fits_apms.R**  <br> 
Assembles the [job_chunk].Rdata into one.  <br> 
Outputs excel files containing summarizing the the results per effect per protein and per contrast per protein. <br> 
The results are compiled into the **Supplementary Table S3, Tab 1 - VZV-Host interactions and Tab 2 - VZV ORF baits**.  <br> 

## Full proteome changes induced by the depletion of the MPP8 gene in SK-N-BE2 cells
requires: **msglm** package v.0.6.0 and **msimportr** <br> 
1. **prepare_data_MPP8_KO.R**  <br> 
Reads the maxquant output files: peptides.txt and evidence.txt  <br> 
Protein groups are refined via the protregroup.jl script -  <br> 
Protein groups distinguished by only one specific peptide or with less than 25% different specific peptides are merged to extend the set of peptides used for protein group quantitation and reduce the protein isoform-specific changes. <br> 
Outputs: <br> 
msfull.Rdata containing the formated and annotated dataset, with protein groups and peptides intensities.  <br> 
msglm.Rdata containing the formated and annotated dataset with peptide intensities, sample associated normalisation shifts, mass spectrometer associated calibration factors,  <br> 
the GLM model description (effects) with corresponding priors.  <br> 
2. **msglm_fit_chunk_MPP8_KO.R** <br> 
Extracts the data of a single protein group based on its protgroup_id [job_chunk] <br> 
Fits the Bayesian model on the protein group. <br> 
Outputs: <br> 
[job_chunk].Rdata containing the msglm results <br> 
3. **msglm_fit-chunks_MPP8_KO.lrz.sh** <br> 
Bash script instructing the fitting of all the protein groups of the dataset [job_chunk] individually. <br> 
4. **assemble_fits_MMP8_KO.R**  <br> 
Assembles the [job_chunk].Rdata into one.  <br> 
Outputs an excel file containing summarizing the the results per effect per protein. <br> 
These results are integrated in the compiled **"Complementary omics datasets" Table S5** <br> 

## Integration of the interactome and effectome data by network diffusion analysis 
requires:**HierarchicalHotNet** <br> 
1. **hotnet_analysis.jl**  <br> 
General script for the hotnet analysis, dependent on results from 2. and 4. <br> 
Maps the Interactome and Effectome hits on ReactomeFI. Tests the significance of the generated edges of the network against permuted data-based networks.  <br> 
2. **hotnet_treestats_chunks.jl** <br> 
Calculates Strongly Connected Component Tree statistic per viral bait.  <br> 
associated bash script: **hotnet_treestats_chunk.lrz.sh**  <br> 
Bash script  <br> 
3. **hotnet_perm_chunk.jl** <br> 
Generate 1000 permuted tree per bait and calculate SCC tree statistics.  <br> 
associated bash script: **hotnet_perm_chunk.lrz.sh**  <br> 
4. **hotnet_perm_chunk_assemble.jl** <br> 
Assembles Hotnet permuted tree results.  <br> 
