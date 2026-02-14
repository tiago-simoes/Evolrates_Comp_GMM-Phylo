# Evolutionary rate incongruences in squamates reveal contrasting patterns of evolutionary novelties and innovation

Dataset DOI: [10.5061/dryad.t76hdr8dc](10.5061/dryad.t76hdr8dc)

## GENERAL INFORMATION

This README_Simoes_etal_2022.txt file was generated on 2025-10-016 by Tiago R.Simões

1.  Title of Dataset: Data from: Evolutionary rate incongruences in squamates reveal contrasting patterns of evolutionary novelties and innovation

2.  Author Information Corresponding Investigator Name: Tiago R. Simões Institution: Princeton University, Princeton-NJ, USA Email: [simoes\@princeton.edu](mailto:simoes@princeton.edu){.email}

    Co-investigator 1 Name: Arthur Brum Institution: Universidade do Estado do Rio de Janeiro, Rio de Janeiro-RJ, Brazil

    Co-investigator 2 Name: Stephanie Pierce Institution: Harvard University, Cambridge-MA, USA

3.  Date of data generation: 2025

4.  Geographic location of data collection: N/A

5.  Funding sources that supported the collection of the data:

U.S. National Science Foundation: 2323124, DEB/SBS

National Council for Scientific and Technological Development: 151134/2024-3

Harvard University

6.  Recommended citation for this dataset: Simões et al. (2025), Supplementary material and supplementary data files for: Evolutionary rate incongruences in squamates reveal contrasting patterns of evolutionary novelties and innovation.

## Description of the data and file structure

The different data types included in this repository are listed below, along with brief descriptions for how to work with them.

## S1_Analyses:

### S1.1_Evo rates: Input data, R codes, and outputs for evolutionary rate analyses

#### S1.1.1_Phylo_Clock:

All codes, data, and outputs for extracting evolutionary rates and calculating selection mode from results from phylodynamic (clock) analyses

1.  Number of folders: 6

2.  Number of files: 88

3.  Missing data codes: None

4.  Main folders list and structure

```         
1Mol1Mor_164t: Results from analyses using 1 morphological and 1 molecular clock partition (global morphological clock) & full taxon sampling (164 taxa)

4Mol4Mor_92t: Results from analyses using 4 morphological and 4 molecular clock partitions & matching GMM taxon sampling with 0% fossils (0% taxa as fossils): 92 taxa

4Mol4Mor_106t: Results from analyses using 4 morphological and 4 molecular clock partitions & matching GMM taxon sampling with 50% fossils (14% taxa as fossils): 106 taxa

4Mol4Mor_122t: Results from analyses using 4 morphological and 4 molecular clock partitions & matching GMM taxon sampling with no fossil removal (24% taxa as fossils): 122 taxa

4Mol4Mor_164t: Results from analyses using 4 morphological and 4 molecular clock partitions & full taxon sampling (164 taxa)
```

5.  Abbreviations used: Mol, molecular data; Mor, morphological data; #t, number of taxa in analysis.

6.  Other relevant information:

    Subfolders: each subfolder will include the following classes of results files: summary output time-trees used to infer rates in Phyllip (.t), vectorized images of the rates in trees (.pdf), and tables with rate values per clade (.csv).

    R scripts: R scripts at the root correspond to the scripts used to produce the outputs into folders with corresponding names

#### S1.1.2_BayesTraits:

All codes, data, and outputs for calculating evolutionary rates for GMM shape data using BayesTraits

1.  Number of folders: 17

2.  Number of files: 169

3.  Missing data codes: None

4.  Main folders list and structure:

```         
phyloPCs: outputs using phylogenetic principal components from GMM

StandardPCs: outputs using stands principal components from GMM
```

5.  Abbreviations used: PC, principal component; MCT, maximum compatible tree; 0F, 0 fossils; 50F, 50% of initial fossils; 100F, 100% of fossils; #t, number of taxa in analysis.

6.  Other relevant information:

    Subfolders: each subfolder will include the following classes of files: BayesTraitsV4 execution files, output folders with log (.txt) and input tree files in nexus format (.nex).

    R scripts: R scripts at the root correspond to the scripts used to produce the outputs into folders with corresponding names

#### S1.1.3_RRphylo:

All codes, data, and outputs for calculating evolutionary rates for GMM shape data using RRphylo

1.  Number of folders: 8

2.  Number of files: 111

3.  Missing data codes: None

4.  Main folders list and structure:

```         
phyloPCs: outputs using phylogenetic principal components from GMM  

StandardPCs: outputs using stands principal components from GMM
```

5.  Abbreviations used: PC, principal component; MCT, maximum compatible tree; 0F, 0 fossils; 50F, 50% of initial fossils; 100F, 100% of fossils; #t, number of taxa in analysis.

6.  Other relevant information:

    Subfolders: each subfolder will include the following classes of files: input tree files in nexus format (.nex)., vectorized images of the rates in trees (.pdf), and tables with rate values per clade (.csv).

    R scripts: R scripts at the root correspond to the scripts used to produce the outputs into folders with corresponding names

#### S1.1.4_Comparisons_Among_Methods:

All codes and outputs for comparing evolutionary rates from approaches 1-3 above.

1.  Number of folders: 0

2.  Number of files: 6

3.  Missing data codes: None

4.  Abbreviations used: PC, principal component; MCT, maximum compatible tree; 0F, 0 fossils; 50F, 50% of initial fossils; 100F, 100% of fossils; #t, number of taxa in analysis

5.  Other relevant information:

    Files: input tree files in nexus format (.nex)., tables with rate values per clade (.csv), R scripts: R script for rate comparisons.

### S1.2_GMM:

Input data, R codes, and outputs for geometric morphometrics, morphological disparity and morphospace analyses

All codes, data, and outputs for calculating evolutionary rates for GMM shape data using RRphylo

1.  Number of folders: 2

2.  Number of files: 89

3.  Missing data codes: None

4.  Main folders list and structure:

```         
Input: TPS Landmark files and cruve files used for geometric morphometrics (generalized Procrustes analysis) and phylogenetic trees used as input for statistical analysis: phylogenetic signal, PCA, phyloPCA, morphospace.

Output: Results produced by GMM analyses using GMM_PCA R script, including figures for PCA, phylogenetically corrected PCA (phyPCA), phylomorphospace, phylogenetic signal (PS), and broken sticks graph (BS) for contribution of each principal component to total variance.
```

5.  Abbreviations used: GMM, geometric morphometrics; MCT, maximum compatible tree; PC, principal component; 0F, 0 fossils; 50F, 50% of initial fossils; 100F, 100% of fossils; #t, number of taxa in analysis

6.  Other relevant information:

    Subfolders: each subfolder will include the following classes of files: input tree files in nexus format (.nex)., landmark files (.tps), vectorized images of the rates in trees (.pdf), and tables with rate values per clade (.csv).

### S1.3_Phylogenetics:

Input data, MrBayes codes, and outputs for all phylogenetics/phylodanymic analyses

1.  Number of folders: 2

2.  Number of files: 89

3.  Missing data codes: None

4.  Main folders list and structure:

```         
Data: TPS Landmark files and cruve files used for geometric morphometrics (generalized Procrustes analysis) and phylogenetic trees used as input for statistical analysis: phylogenetic signal, PCA, phyloPCA, morphospace.

--morphological: Morphological character list i(.docx and .pdf), morphological dataset in Nexus format (.nex)

--molecular: Concatenated molecular dataset in Nexus format (.nex), molecular partitions file (.txt)

Output: Results produced by GMM analyses using GMM_PCA R script, including figures for PCA, phylogenetically corrected PCA (phyPCA), phylomorphospace, phylogenetic signal (PS), and broken sticks graph (BS) for contribution of each principal component to total variance.

--Bayesian: input files including MrBayes input file with dataset and MrBayes commands block in Nexus format (.nex), log files in text format (.txt), parameter files (.p), and tree files in Phyllip format (.t).

--BayesianClocks: input files including MrBayes input file with dataset and MrBayes commands block in Nexus format (.nex), log files in text format (.txt), parameter files (.p), and tree files in Phyllip format (.t).

--MolecularTree: input files including MrBayes input file with dataset and MrBayes commands block in Nexus format (.nex), log files in text format (.txt), parameter files (.p), and tree files in Phyllip format (.t).
```

5.  Abbreviations used: Mol, molecular data; Mor, morphological data; #t, number of taxa in analysis.

## S2_R codes:

All R scripts used in this study (they are also provided with the associated input data and output files in their respective directories in the Analyses folder)

1.  Number of folders: 0

2.  Number of files: 13

3.  Missing data codes: None

## S3_Supplementary data tables:

All supplementary data tables (S1-S44) in Excel format

1.  Number of folders: 0

2.  Number of files: 3

3.  Missing data codes: None

### Files and variables

#### File: README.txt

**Description:** Hierarchical README file

#### File: S1_Analyses.zip

**Description:** All analyticial input and output files

#### **File: S2_R_codes.zip**

**Description:** All R scripts used in this study

#### File: S3_Suppl_Data_Tables.zip

**Description:** All supplementary data tables (S1-S44) in Excel format

#### 

## Code/software

R v.4.1.1; MrBayes v.3.2.7; BayesTraits v.4.

## Access information

Other publicly accessible locations of the data:

-   <https://github.com/tiago-simoes/BTprocessR2>
-   <https://github.com/tiago-simoes/EvoPhylo>.
