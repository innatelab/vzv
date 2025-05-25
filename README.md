# Data analysis scripts for *"Varizella-Zoster Virus proteomic profiling" (Girault et al., 2025)* study

## In-house Bioinformatics Packages

These in-house *R* and *Julia* packages are used by the VZV data analysis scripts:

- [**msglm**](https://github.com/innatelab/msglm), Bayesian Mixed Effects Linear Models for Mass Spectrometry Data (*R*):
  - [v0.5.0](https://github.com/innatelab/msglm/tree/v0.5.0)
  - [v0.6.0](https://github.com/innatelab/msglm/tree/v0.6.0)
- [**msimportr**](https://github.com/innatelab/msimportr), importing and pre-processing of MaxQuant/Spectronaut output (*R*):
  - [v0.3.0](https://github.com/innatelab/msimportr/tree/v0.3.0)
- [**HierarchicalHotNet.jl**](https://github.com/alyst/HierarchicalHotNet.jl),
  Julia implementation of [Hierarchical HotNet method](https://academic.oup.com/bioinformatics/article/34/17/i972/5093236):
  - [v1.1.0](https://github.com/alyst/HierarchicalHotNet.jl/tree/v1.1.0)

## Analysis of VZV Interactome AP-MS Data

The analysis of affinity-purification data of V5-tagged VZV proteins in SK-N-BE2 cells requires
**msglm** (v0.5.0) and **msimportr** (v0.3.0) packages and is performed in the following steps:

1. **prepare_data_apms.R** script
    - *Reads* the MaxQuant files *ProteinGroups.txt* and *evidence.txt* and prepares the data for statistical analysis.
    - *Outputs*:
      - **msfull.RData** containing the formated and annotated dataset, with protein groups and peptides intensities.
      - **msglm.RData** containing the formated and annotated dataset with protein groups intensities,
        per-MS run intensity normalisation factors, MS instrument noise model parameters,
        the GLM model description (matrices of effects, batch effects, and contrasts) including parameter priors.
2. **msglm_fit_chunk_apms.R** script is called for each protein group (*job chunk*) in the data set
    - Gets `protgroup_id` as an input
    - Extracts the data of a given a protein group from **msglm.RData**
    - Fits the Bayesian model defined **msglm.RData** using *msglm* package
    - *Outputs* **<job_chunk>.RData** file with the *msglm* results
3. **msglm_fit-chunks_apms.sge.sh** is a shell script for the parallel execution of **msglm_fit_chunk_apms.R**
 on the compute cluster using [SGE job scheduler](https://computing.sas.upenn.edu/gpc/job/sge).
4. **assemble_fits_apms.R** script assembles all individual **<job_chunk>.RData** files into a single report.
    - Extracts the significance of relevant model contrasts for each protein group.
    - Compiles the **Supplementary Table S3, Tab 1 - VZV-Host interactions and Tab 2 - VZV ORF baits** report.

## Analysis of MPP8-induced proteome changes

The analysis of proteomic changes induced by the depletion of the MPP8 gene in SK-N-BE2 cells
requires **msglm** package (v.0.6.0) and **msimportr** (v0.3.0).
This analysis closely follows the same steps as the analysis of AP-MS data described above:

1. **prepare_data_MPP8_KO.R**
    - Reads the MaxQuant files: *peptides.txt* and *evidence.txt*
    - Corrects protein groups using the **protregroup.jl** script to avoid splitting isoforms of the same
      gene into individual protein groups that differ only by a single specific peptide
      (see *Material and Methods* section for details).
    - *Outputs*:
      - **msfull.RData** containing the formated and annotated dataset
        with protein groups and peptides intensities.
      - **msglm.RData** containing the formated and annotated dataset with peptide intensities,
        per-MS run intensity normalisation factors, MS instrument noise model parameters,
        the GLM model description (matrices of effects, batch effects, and contrasts) including parameter priors.
2. **msglm_fit_chunk_MPP8_KO.R** script is called for each protein group (*job chunk*) in the data set
    - Gets `protgroup_id` as an input
    - Extracts the data of a given a protein group from **msglm.RData**
    - Fits the Bayesian model defined **msglm.RData** using *msglm* package
    - *Outputs* **<job_chunk>.RData** file with the *msglm* results
3. **msglm_fit-chunks_MPP8_KO.lrz.sh** is a shell script for the parallel execution of
   **msglm_fit_chunk_MPP8_KO.R** on the compute cluster using
   [SLURM workload manager](https://slurm.schedmd.com/sbatch.html).
4. **assemble_fits_MMP8_KO.R** script assembles all individual **<job_chunk>.RData**
   files into a single report.
    - Extracts the significance of relevant model contrasts for each protein group.
    - Compiles the **"Complementary omics datasets" Table S5**

## Integration of the interactome and effectome data

The integration of VZV virus-host *interactome* (virus-host protein interaction)
and *effectome* (proteomic changes induced by the expression of individual viral proteins)
is done by diffusing the protein abundance perturbations over the global network of
protein-protein and functional gene interactions within the host cell
([*ReactomeFI* database](https://reactome.org/tools/reactome-fiviz) is used).
It is implemented in [Julia](https://julialang.org/) and uses in-house
[**HierarchicalHotNet.jl**](https://github.com/alyst/HierarchicalHotNet.jl) package
that implements *Hierarchical HotNet* method and provides additional statistics
for the network diffusion process.

1. **hotnet_analysis.jl** is the general script for the HotNet analysis:
    - Reads the interactome and effectome (the results of the *msglm* analysis).
    - Reads the *ReacomeFI* network of protein-protein and functional gene interactions.
    - Conducts the *HotNet* analysis for unperturbed interactome and effectome data
    - Prepares the nodes and edges weights for the *HotNet* network diffusion
      analysis of each viral bait (step 2).
    - Prepares the 1000 random permutations of node and edge weights per viral protein for step 3
      and saves them in **hotnet_perm_input.jlser.zst** file.
    - Reads the network diffusion results from step 2 and
      the combined permutation statistics from step 4 (**hotnet_perm_assembled_<viral_protein>.jlser.zst**).
    - Identifies the significant interactions between the host proteins physically associated
      with viral proteins (interactome) and the proteins effected by virial proteins overexpression
      (effectome)
    - *Outputs* the **Supplementary Table S4 - HotNet analysis results (FIXME)** with the
      significant interactions between the host proteins and viral proteins.
2. **hotnet_treestats_chunks.jl**
    - Performs *HotNet* network diffusion to integrate interactome and effectome for each viral protein.
    - Calculates the *Strongly Connected Component Tree* statistics for each edge weight cutoff threshold.
    - *Outputs* the **hotnet_treestats_<viral_protein>.jlser.zst** file with the
      *HotNet* network diffusion results for each viral protein.
    - **hotnet_treestats_chunk.lrz.sh** is the associated *SLURM* job script for parallel execution on
      the compute cluster.
3. **hotnet_perm_chunk.jl**
    - Reads the random interactome and effectome permutations from **hotnet_perm_input.jlser.zst**.
    - Performs memory and computationally intensive network diffusion for a block of random
      interactome and effectome permutations generated at step 1.
    - Calculates the *Strongly Connected Component Tree* statistics for each edge weight cutoff threshold.
    - *Outputs* the **hotnet_perm_<job_chunk>.jlser.zst** file.
    - **hotnet_perm_chunk.lrz.sh** is the associated *SLURM* shell script for the parallel
      execution on the compute cluster.
4. **hotnet_perm_chunk_assemble.jl**
    - assembles the HotNet permuted tree results (**hotnet_perm_<job_chunk>.jlser.zst**)
      generated at step 3
    - *Outputs* the combined permutation statistics (**hotnet_perm_assembled_<viral_protein>.jlser.zst**).
