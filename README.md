# "Differential responses of the wheat microbiome to Fusarium graminearum and its virulence factor Deoxynivalenol"
### Whitaker, B.K., Rhoades, N., McCormick, S.P., Hao, G.

**zenodo doi**

This repository includes the R code, data files, and small scripts to supplement the manuscript by Whitaker et al. "Differential responses of the wheat microbiome to Fusarium graminearum and its virulence factor Deoxynivalenol".

For the Round 1 (R1) and Round 2 (R2) projects, spring wheat plants (variety Norm) were grown outside in pots and exposed to natural microbiome sources for infection. At maturity, they were brought inside and placed in growth chambers (C1 or C2) where they were inoculated with either a wild type Fusarium graminearum (PH-1 strain), its trichothecene gene 5 knockout unable to produce Deoxynivalenol (PH-1 delta Tri5), or a mock inoculation as control. Disease progress was monitored and wheat heads harvested 7-10 days post infection. Samples were then ground and processed for microbiome determination using 16S-V4 or ITS2 gene regions on an Illumina Miseq v3 (Round1) or NextSeq 2000 P2 (Round2). For the R2 sequencing, the Illumina NextSeq P2 allowed for more project samples to be grouped on one run. Those projects were called BarleyEndos and HackFlex in the metadata and bioinformatic processing files.

The sample metadata are available in these files ("Tri5MB_R1_Dec2023_MBmetadata", "R2_Tri5MB_BarleyEndos_HackFlex_MBMetadata.csv", and "Tri5MB-R2_BarleyEndos_IlluminaPrep_plateLayout.csv").

Disease data are available in these files: "Tri5MB_R1_diseaseHeads.csv", "Tri5MB_R1_diseaseTubes.csv", "Tri5MB_R2_diseaseHeads.csv", "Tri5MB_R2_diseaseTubes.csv", 

The raw sequence data used are available through the NCBI SRA database under accession numbers PRJNA1559029 (R1) and PRJNA1557970 (R2). https://www.ncbi.nlm.nih.gov/sra/PRJNA1559029 and https://www.ncbi.nlm.nih.gov/sra/PRJNA1557970



Bioinformatic processing using DADA2: 
* DADA2_R1Bact.Rmd & DADA2_R1Fung.Rmd
* ./code/DADA2_hpcc_R2bacteria.R & ./code/DADA2_hpcc_R2fungi.R 

Additional Data Cleaning:
* Taxonomy_Curation.Rmd does taxonomic classification of ASVs, removal of plant ASVs, and curation of contaminant ASVs based on analyses of negative controls. 
* ps_Organize.Rmd was used to merge and process the various rounds and mixed projects included after DADA2 processing. 

Core Data Analyses:
* CommunityAnalysis.Rmd analyzes Community Structure and Diversity of the microbial communities, as well as disease data analysis; 
* DESeq2.Rmd performs differential abundance analyses 
* CCA.Rmd performs a direct gradient analysis, aka a Canonical correspondence analysis

Run files in this order:
(various) DADA2**pipelines > Taxonomy_Curation.Rmd > ps_Organize.Rmd > CommunityAnalysis.Rmd > DESeq2.Rmd > CCA.Rmd 

The /code folder contains the batch and R scripts necessary to run small codes in R and on HPC. Data necessary to run the analyses can be found in /data folder.

Please see the manuscript for details and full reference information.
