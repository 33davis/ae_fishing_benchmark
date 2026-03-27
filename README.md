
### methodology and data construction of the project 
1. Identify current bioinformatics standards and methods
2. Create synthetic test data
   - created a database of 20 COI partial sequences for Nototheniideae and Sphenisciformes from NCBI
   - I chose COI sequences based on availability, as well as hoping to utilise the new MARES COI barcode database pipeline that was recently published.
    2a. MARES pipeline: for initial testing, the precompiled database from [https://osf.io/8rdqk/] was used, which was downloaded on 20.02.2024.
   [33davis/MARES_database_pipeline]
    2b. Create a combined MARES + PR2v5 database to simulate a complex eukaryotic database to map published data against
    2c. Understand utility and proof of concept using this database
3. create "fishing" dataset combining empirical and synthetic data
   - utilising gargammel to simulate sedaDNA fragments [add link to code]
   - for COI sequence "synthetic" sequences, utilise the bioconda::gargammel-slim package to create random fragments and simulate damage on partial gene sequences (no gaps)
   - create different "damage level" sequences by subsampling created files [add link to code]
   - run gargammel on targeted files (all deamination files, and subsampled files [add link to code]
   - merge fasta files + combined files back together prior to running through pipelein [add link to code]
   - these different datasets will then be run through target pipelines by themselves, and then also merged with the published dataset to quantify and compare the ability to be "fished" back out. This will help provide reference if simulated ancient Nototheniideae and Sphenisciformes sequences can be identified with the current pipeline and settings, especially in complex metagenomic samples
   - will also highlight the importance of database choice 
4. test data with Linda pipeline (Armbrecht et al., 2020) [https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13162]
   #### It is noted that MALT/ HOPS has not had an update with NCBI since 2022. Most likely, it will not be updated in the future. 
     - beginning with the U1538 site data
     - will begin with U1538 sample data, that is already been collapsed + adapterremoval + intial FastQC/MultiQC qaulity control
       (6a) create conda environment for testing data  available here
       (6b) databases were merged with code 
       (6c) malt-build with  for MARES preload database, and merged MARES + PR2 v5 database
       (6d) run Komplexity - bbtools/dedupe - MALT - blast2RMA  
       (6e) Visualise and quantify runs in MEGAN v. against built databases // previously published data.
  5. Binary Classifier analyses to establish patterns + lowest minimum thresholds based on fragment amounts
   

### Flow chart and code location ### 
#### Pre-procesing ####  --> /sedaDNA_pipeline 
#### synthetic and IODP + synthetic dataset construction #### ---> /synthetic_dataests
#### dataset processing in bioinformatic pipeline #### ---> /sedaDNA_pipeline
#### HOPS + percent damage #### ---> /sedaDNA_pipeline
#### Basic identities + statistics #### ---> /analyses
#### synthetic extraction + visualization #### ---> /analyses
#### Precision, Sensitivity, and F1 Scores #### ---> /analyses
#### Statistical Analyses, Threshold Estimation, and Visualisation #### ---> /analyses
