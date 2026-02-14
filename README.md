# Evolutionary rate incongruences in squamates reveal contrasting patterns of evolutionary novelties and innovation

Dataset DOI: [10.5061/dryad.t76hdr8dc](https://doi.org/10.5061/dryad.t76hdr8dc)

## GENERAL INFORMATION

This README file was generated on 2025-10-16 by Tiago R. Simoes

1.  Title of Dataset: Data from: Evolutionary rate incongruences in squamates reveal contrasting patterns of evolutionary novelties and innovation

2.  Author Information Corresponding Investigator Name: Tiago R. Simoes Institution: Princeton University, Princeton-NJ, USA Email: [simoes\@princeton.edu](mailto:simoes@princeton.edu){.email}

    Co-investigator 1 Name: Arthur Brum Institution: Universidade do Estado do Rio de Janeiro, Rio de Janeiro-RJ, Brazil

    Co-investigator 2 Name: Stephanie Pierce Institution: Harvard University, Cambridge-MA, USA

3.  Date of data generation: 2025

4.  Geographic location of data collection: N/A

5.  Funding sources that supported the collection of the data:

    U.S. National Science Foundation: 2323124, DEB/SBS

    National Council for Scientific and Technological Development: 151134/2024-3

    Harvard University

6.  Recommended citation for this dataset: Simoes et al. (2025), Supplementary material and supplementary data files for: Evolutionary rate incongruences in squamates reveal contrasting patterns of evolutionary novelties and innovation.

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

    Subfolders: each subfolder will include the following classes of results files: summary output time-trees used to infer rates in Phyllip (.t), vectorized images of the rates in trees (.pdf), and tables with rate values per clade (.csv).

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

    Subfolders: each subfolder will include the following classes of files: BayesTraitsV4 execution files, output folders with log (.txt) and input tree files in nexus format (.nex).

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

    Subfolders: each subfolder will include the following classes of files: input tree files in nexus format (.nex)., vectorized images of the rates in trees (.pdf), and tables with rate values per clade (.csv).

    R scripts: R scripts at the root correspond to the scripts used to produce the outputs into folders with corresponding names

#### S1.1.4_Comparisons_Among_Methods:

All codes and outputs for comparing evolutionary rates from approaches 1-3 above.

1.  Number of folders: 0

2.  Number of files: 6

3.  Missing data codes: None

4.  Abbreviations used: PC, principal component; MCT, maximum compatible tree; 0F, 0 fossils; 50F, 50% of initial fossils; 100F, 100% of fossils; #t, number of taxa in analysis

5.  Other relevant information:

    Files: input tree files in nexus format (.nex)., tables with rate values per clade (.csv), R scripts: R script for rate comparisons.

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

All supplementary data tables (S1-S44) in Excel format, distributed across three .xlsx files.

1.  Number of folders: 0

2.  Number of files: 3

3.  General abbreviations used across tables: GMM, geometric morphometrics; PCA, principal component analysis; pPCA, phylogenetic PCA; PACA, phylogenetic-aligned component analysis; PC, principal component; CS, centroid size; spp, species; 3D, three-dimensional; t, number of taxa; sub50, 50% fossil subsampling; sub0, 0% fossil subsampling (extant taxa only).

------------------------------------------------------------------------

### File 1: Suppl_Data_Tables S1-S7_DataSamples.xlsx

**Description:** Taxon sampling, specimen data, landmark definitions, and fossil calibration ages used as inputs across all analyses. Contains 9 sheets (7 data tables + 1 institutional abbreviations sheet + 1 references sheet).

------------------------------------------------------------------------

#### Sheet: S1_Phylo_Morpho sample (Table S1)

**Description:** Morphological taxon sampling for phylogenetic and phylodynamic (clock) analyses. Lists all 164 taxa included in the morphological dataset.

-   **Rows:** 165 (1 header + 164 taxa)
-   **Columns:** 9

| Variable      | Definition                                                                                                                     |
|------------------------------------|------------------------------------|
| Genus         | Genus name of the taxon                                                                                                        |
| Species       | Species epithet of the taxon                                                                                                   |
| genus_species | Concatenated genus and species name                                                                                            |
| Extant/Fossil | Whether the taxon is extant ("E") or fossil ("F")                                                                              |
| "Suborder"    | Higher-level taxonomic group (e.g., Gekkota, Serpentes); in quotes because these ranks are informal                            |
| "Superfamily" | Superfamily-level classification; in quotes because these ranks are informal                                                   |
| Family        | Family-level classification                                                                                                    |
| Subfamily     | Subfamily-level classification. Note: subfamily assignment is not applicable or has not been formally designated for some taxa |
| Accesion \#   | Museum accession/specimen number(s) for the morphological material examined                                                    |

------------------------------------------------------------------------

#### Sheet: S2_Phylo_Mol sample (Table S2)

**Description:** Molecular sequence sampling for phylogenetic analyses. Lists all molecular loci and GenBank accession numbers used.

-   **Rows:** 659 (1 header + 658 sequence entries)
-   **Columns:** 6

| Variable         | Definition                                                        |
|------------------------------------|------------------------------------|
| organism         | Taxon name (genus + species)                                      |
| accession_number | GenBank accession number for the sequence                         |
| locus            | Name of the molecular locus (e.g., 12S_rRNA, 16S_rRNA, CYTB_mDNA) |
| length           | Sequence length in base pairs (without alignment gaps)            |
| length_with_gaps | Sequence length in base pairs (including alignment gaps)          |
| specimen_voucher | Museum voucher specimen number associated with the sequence.      |

------------------------------------------------------------------------

#### Sheet: S3_GMM sample (Table S3)

**Description:** Taxon sampling for geometric morphometric (GMM) analyses. Lists all 122 taxa included in cranial shape analyses.

-   **Rows:** 123 (1 header + 122 taxa)
-   **Columns:** 9

| Variable       | Definition                                                                    |
|------------------------------------|------------------------------------|
| File_Name      | File name identifier for the specimen (genus_species format)                  |
| Genus          | Genus name                                                                    |
| Taxon          | Full taxon name (genus + species)                                             |
| Group          | Major clade assignment (e.g., Serpentes, Scincoidea, Gekkota)                 |
| Group2         | More specific clade assignment (e.g., Caenophidia within Serpentes)           |
| Family         | Family-level classification                                                   |
| Specimen       | Museum specimen number.                                                       |
| Original_media | Type of media used for morphometric data capture (e.g., mesh, photo, CT-scan) |
| Material       | Whether the taxon is extant or fossil                                         |

------------------------------------------------------------------------

#### Sheet: S4_FossilPreserv (Table S4)

**Description:** Fossil specimen preservation data. Documents the number of preserved anatomical elements available for each fossil taxon, used to assess sampling completeness.

-   **Rows:** 91 (1 header + 89 taxa + 1 blank separator row + 1 TOTAL summary row)
-   **Columns:** 12

| Variable          | Definition                                                   |
|------------------------------------|------------------------------------|
| Taxa              | Taxon name (genus + species)                                 |
| Taxon sample      | Total number of specimens examined for that taxon            |
| Skull (lateral)   | Number of specimens preserving the skull in lateral view     |
| Skull (dorsal)    | Number of specimens preserving the skull in dorsal view      |
| Skull (ventral)   | Number of specimens preserving the skull in ventral view     |
| Mandible (lat)    | Number of specimens preserving the mandible in lateral view  |
| Mandible (medial) | Number of specimens preserving the mandible in medial view   |
| Dentary (lat)     | Number of specimens preserving the dentary in lateral view   |
| Dentary (medial)  | Number of specimens preserving the dentary in medial view    |
| SVL               | Number of specimens preserving snout-vent length measurement |
| Humeri            | Number of specimens preserving humeri                        |
| Femora            | Number of specimens preserving femora                        |

**Empty cells (12 cells):** Row 90 is an entirely blank row that serves as a visual separator between the last taxon data row and the TOTAL summary row (row 91). These cells are intentionally empty for formatting purposes.

------------------------------------------------------------------------

#### Sheet: S5_Fossil_Subsampling (Table S5)

**Description:** Fossil taxon subsampling scheme. Lists all 30 fossil species in the full fossil sample and indicates which 15 were retained in the 50% subsample, along with their higher-level clade assignments.

-   **Rows:** 31 (1 header + 30 taxa)
-   **Columns:** 3 effective data columns (A, B, C)

| Variable                                                             | Definition                                                                                                                                                                                                                                     |
|------------------------------------|------------------------------------|
| Column A (header: "100% of fossils = 30 spp (3D: 22/30)")            | Taxon names for all 30 fossil species in the full (100%) fossil sample. "3D: 22/30" indicates 22 of 30 had 3D data available                                                                                                                   |
| Column B (header: "Subsample (50% of fossils) = 15 spp (3D: 12/15)") | Taxon names for the 15 fossil species retained in the 50% subsample. **Empty cells:** cells are empty for the 15 taxa that were excluded from the 50% subsample (i.e., taxa listed only in column A but not selected for the reduced sampling) |
| Column C (no header)                                                 | Higher-level clade assignment for each taxon (e.g., Anguiformes, Sphenodontia, Serpentes, Gekkota, Lacertoidea, Iguania, Mosasauria, Non-lepidosaurs, Early Lepidosauria, Indet)                                                               |

------------------------------------------------------------------------

#### Sheet: S6_Landmarks (Table S6)

**Description:** Landmark definitions for the 2D geometric morphometric analysis of cranial shape.

-   **Rows:** 25 (1 header + 24 landmark definitions)
-   **Columns:** 3

| Variable         | Definition                                                                                                       |
|------------------------------------|------------------------------------|
| Landmark numbers | Landmark number(s) as used in the TPS landmark files                                                             |
| Type             | Landmark type following Bookstein (1991): I = anatomical (Type I) landmark at discrete juxtapositions of tissues |
| Definition       | Anatomical definition of the landmark position on the cranium                                                    |

------------------------------------------------------------------------

#### Sheet: S7_FossilCalAges (Table S7)

**Description:** Fossil calibration ages used for divergence time estimation. Provides stratigraphic and chronostratigraphic information for each fossil taxon, with minimum and maximum ages based on the International Stratigraphic Chart.

-   **Rows:** 75 (1 header + 74 taxa)
-   **Columns:** 7

| Variable            | Definition                                                                      |
|------------------------------------|------------------------------------|
| Genus               | Genus name of the fossil taxon                                                  |
| species             | Species epithet of the fossil taxon                                             |
| Stratigraphic level | Geological formation, stratigraphic unit, and locality information.             |
| Chronostratigraphy  | Chronostratigraphic age assignment (e.g., "Induan-Olenekian, Early Triassic").  |
| MaxAge (ISC 2021)   | Maximum age in millions of years ago (Ma)                                       |
| MinAge (ISC 2021)   | Minimum age in millions of years ago (Ma)                                       |
| Main age references | Literature references for the age assignments                                   |

------------------------------------------------------------------------

#### Sheet: InstAbbrev

**Description:** List of institutional abbreviations referenced throughout the supplementary data tables. Not a numbered supplementary table; provided as a reference key for specimen accession numbers.

-   **Rows:** 63 (1 header + 1 blank separator row + 61 abbreviation entries)
-   **Columns:** 1

------------------------------------------------------------------------

#### Sheet: References

**Description:** List of literature references cited in the supplementary data tables. Not a numbered supplementary table.

-   **Rows:** 67
-   **Columns:** 1

------------------------------------------------------------------------

### File 2: Suppl_Data_Tables S8-37_Results_PCA.xlsx

**Description:** Results from geometric morphometric analyses, including standard PCA, phylogenetic PCA (pPCA), phylogenetic-aligned component analysis (PACA), Procrustes coordinates, and ancestral shape reconstructions. Contains 30 sheets organized into three groups corresponding to three taxon samplings:

-   **Full sample** (122 taxa, 100% fossils = 30 fossil spp): Sheets S8, S11--S19
-   **50% fossil subsample** (106 taxa, 50% fossils = 15 fossil spp): Sheets S9, S20--S28
-   **0% fossils / extant only** (92 taxa, 0 fossil spp): Sheets S10, S29--S37
-   

------------------------------------------------------------------------

#### Sheets: S8_Classifiers, S9_Classifiers_sub50, S10_Classifiers_sub0 (Tables S8, S9, S10)

**Description:** Classifier/metadata tables for each taxon sampling scheme, used to label and color-code taxa in morphospace plots.

-   **Rows:** 123 / 107 / 93 (1 header + 122 / 106 / 92 taxa, respectively)
-   **Columns:** 7

| Variable  | Definition                                                                     |
|------------------------------------|------------------------------------|
| File_Name | File name identifier for the specimen (genus_species format)                   |
| Number    | Sequential numeric identifier for each taxon                                   |
| Genus     | Genus name                                                                     |
| Group.1   | Major clade assignment with subclade detail (e.g., "Serpentes (Caenophidia)")  |
| Group.2   | Subclade or clade assignment used for grouping (e.g., Caenophidia, Scincoidea) |
| Color     | Hex color code assigned to each taxon for plotting                             |
| Material  | Whether the taxon is extant or fossil                                          |

------------------------------------------------------------------------

#### Sheets: S11_Coords_CS, S20_Coords_CS_sub50, S29_Coords_CS_sub0 (Tables S11, S20, S29)

**Description:** Procrustes-aligned landmark coordinates and centroid size values from the generalized Procrustes analysis (GPA) for each taxon sampling scheme.

-   **Rows:** 123 / 107 / 93 (1 header + 122 / 106 / 92 taxa, respectively)
-   **Columns:** 277

| Variable                                              | Definition                                                                                                                                                                                                                                 |
|------------------------------------|------------------------------------|
| File_Name                                             | File name identifier for the specimen                                                                                                                                                                                                      |
| Centroid_size                                         | Centroid size of the landmark configuration (square root of the sum of squared distances of each landmark from the centroid)                                                                                                               |
| Log.CS                                                | Natural logarithm of centroid size                                                                                                                                                                                                         |
| #.X, #.Y (e.g., 1.X, 1.Y, 2.X, 2.Y, ... 137.X, 137.Y) | Procrustes-aligned X and Y coordinates for each of the 137 landmarks (including semilandmarks). Landmark numbering corresponds to the landmark definitions in Table S6 (landmarks 1--24) plus additional semilandmarks placed along curves |

------------------------------------------------------------------------

#### Sheets: S12_PCA, PS21_PCA_sub50, S30_PCA_sub0 (Tables S12, S21, S30)

**Description:** Standard (non-phylogenetic) principal component analysis (PCA) results. Each row represents a taxon with its classifier metadata and scores on all retained principal components.

-   **Rows:** 123 / 107 / 93 (1 header + 122 / 106 / 92 taxa, respectively)
-   **Columns:** 128 / 112 / 98

| Variable                  | Definition                                                                                                                                                                  |
|------------------------------------|------------------------------------|
| File_Name                 | File name identifier for the specimen                                                                                                                                       |
| Number                    | Sequential numeric identifier                                                                                                                                               |
| Genus                     | Genus name                                                                                                                                                                  |
| Group.1                   | Major clade assignment with subclade detail                                                                                                                                 |
| Group.2                   | Subclade or clade assignment                                                                                                                                                |
| Color                     | Hex color code for plotting                                                                                                                                                 |
| Material                  | Extant or fossil                                                                                                                                                            |
| Comp1, Comp2, ... Comp*N* | Scores on each principal component (PC). Number of components equals N-1 where N is the number of taxa: Comp1--Comp121 (full), Comp1--Comp105 (sub50), Comp1--Comp91 (sub0) |

------------------------------------------------------------------------

#### Sheets: S13_PCA_scores, S22_PCA_scores_sub50, S31_PCA_scores_sub0 (Tables S13, S22, S31)

**Description:** Eigenvalues and variance explained by each principal component from the standard PCA.

-   **Rows:** 122 / 106 / 92 (1 header + 121 / 105 / 91 components, respectively)
-   **Columns:** 6

| Variable               | Definition                                                    |
|------------------------------------|------------------------------------|
| Components             | Principal component identifier (e.g., Comp.1, Comp.2, etc.)   |
| d                      | Singular value (square root of eigenvalue) for each component |
| sdv                    | Standard deviation explained by the component                 |
| var                    | Variance (eigenvalue) explained by the component              |
| proportion.of.variance | Proportion of total variance explained by the component       |
| cumulative.proportion  | Cumulative proportion of total variance explained             |

------------------------------------------------------------------------

#### Sheets: S14_PCA_loadings, S23_PCA_loadings_sub50, S32_PCA_loadings_sub0 (Tables S14, S23, S32)

**Description:** Loading values (eigenvectors) for each landmark coordinate on each principal component from the standard PCA.

-   **Rows:** 275 (1 header + 274 landmark coordinate variables)
-   **Columns:** 122 / 106 / 92

| Variable                  | Definition                                                                                                                       |
|------------------------------------|------------------------------------|
| Loadings                  | Landmark coordinate identifier (e.g., 1.X, 1.Y, 2.X, ... 137.X, 137.Y) corresponding to the X and Y coordinates of each landmark |
| Comp1, Comp2, ... Comp*N* | Loading value of each landmark coordinate on the respective principal component                                                  |

------------------------------------------------------------------------

#### Sheets: S15_pPCA, S24_pPCA_sub50, S33_pPCA_sub0 (Tables S15, S24, S33)

**Description:** Phylogenetic principal component analysis (pPCA) results. Same structure as the standard PCA tables (S12/S21/S30) but computed accounting for phylogenetic non-independence among taxa.

-   **Rows:** 123 / 107 / 93 (1 header + 122 / 106 / 92 taxa, respectively)
-   **Columns:** 128 / 112 / 98

| Variable                  | Definition                                      |
|------------------------------------|------------------------------------|
| File_Name                 | File name identifier for the specimen           |
| Number                    | Sequential numeric identifier                   |
| Genus                     | Genus name                                      |
| Group.1                   | Major clade assignment with subclade detail     |
| Group.2                   | Subclade or clade assignment                    |
| Color                     | Hex color code for plotting                     |
| Material                  | Extant or fossil                                |
| Comp1, Comp2, ... Comp*N* | Scores on each phylogenetic principal component |

------------------------------------------------------------------------

#### Sheets: S16_pPCA_scores, S25_pPCA_scores_sub50, S34_pPCA_scores_sub0 (Tables S16, S25, S34)

**Description:** Eigenvalues and variance explained by each component from the phylogenetic PCA. Same structure and variable definitions as Tables S13/S22/S31.

-   **Rows:** 122 / 106 / 92 (1 header + 121 / 105 / 91 components, respectively)
-   **Columns:** 6

| Variable               | Definition                                        |
|------------------------------------|------------------------------------|
| Components             | Phylogenetic principal component identifier       |
| d                      | Singular value for each component                 |
| sdv                    | Standard deviation explained by the component     |
| var                    | Variance (eigenvalue) explained by the component  |
| proportion.of.variance | Proportion of total variance explained            |
| cumulative.proportion  | Cumulative proportion of total variance explained |

------------------------------------------------------------------------

#### Sheets: S17_pPCA_loadings, S27_pPCA_loadings_sub50, S35_pPCA_loadings_sub0 (Tables S17, S27, S35)

**Description:** Loading values (eigenvectors) for each landmark coordinate on each phylogenetic principal component. Same structure and variable definitions as Tables S14/S23/S32.

-   **Rows:** 275 (1 header + 274 landmark coordinate variables)
-   **Columns:** 122 / 106 / 92

| Variable                  | Definition                                                                                   |
|------------------------------------|------------------------------------|
| Loadings                  | Landmark coordinate identifier (e.g., 1.X, 1.Y, ... 137.X, 137.Y)                            |
| Comp1, Comp2, ... Comp*N* | Loading value of each landmark coordinate on the respective phylogenetic principal component |

------------------------------------------------------------------------

#### Sheets: S18_PACA, S26_PACA_sub50, S36_PACA_sub0 (Tables S18, S26, S36)

**Description:** Phylogenetic-aligned component analysis (PACA) results. PACA identifies axes of morphospace that maximize phylogenetic signal.

-   **Rows:** 123 / 107 / 93 (1 header + 122 / 106 / 92 components, respectively)
-   **Columns:** 9

| Variable               | Definition                                                                                             |
|------------------------------------|------------------------------------|
| Component              | Component identifier                                                                                   |
| d                      | Singular value for each component                                                                      |
| sdv                    | Standard deviation                                                                                     |
| var                    | Variance explained by the component                                                                    |
| proportion.of.variance | Proportion of total variance explained                                                                 |
| cumulative.proportion  | Cumulative proportion of total variance explained                                                      |
| RV                     | RV coefficient (multivariate measure of association between the component and the phylogeny)           |
| cumulative.RV          | Cumulative RV coefficient                                                                              |
| K.by.p                 | Blomberg's K statistic divided by the number of variables, measuring phylogenetic signal per component |

------------------------------------------------------------------------

#### Sheets: S19_anc_pPCA, S28_anc_pPCA_sub50, S37_anc_pPCA_sub0 (Tables S19, S28, S37)

**Description:** Ancestral shape reconstructions from phylogenetic PCA. Provides estimated Procrustes-aligned landmark coordinates for each internal node of the phylogeny, reconstructed using the phylogenetic PCA.

-   **Rows:** 122 / 106 / 92 (1 header + 121 / 105 / 91 internal nodes, respectively)
-   **Columns:** 275

| Variable                                              | Definition                                                                                                |
|------------------------------------|------------------------------------|
| Node                                                  | Internal node identifier from the phylogenetic tree                                                       |
| #.X, #.Y (e.g., 1.X, 1.Y, 2.X, 2.Y, ... 137.X, 137.Y) | Reconstructed ancestral Procrustes-aligned X and Y coordinates for each of the 137 landmarks at that node |

------------------------------------------------------------------------

### File 3: Suppl_Data_Tables_S38-44_SumStats_NormRates.xlsx

**Description:** Summary statistics for normalized evolutionary rates across different analytical methods, clock partitions, and taxon samplings. Contains 7 sheets. Summary statistics describe the distribution of normalized (z-transformed) evolutionary rates across clades.

------------------------------------------------------------------------

#### Sheet: S38_Clock_GlobalClock (Table S38)

**Description:** Summary statistics for normalized evolutionary rates from the global (1-partition) morphological clock analysis using the full taxon sampling (164 taxa).

-   **Rows:** 13 (1 header + 12 clade-by-clock entries)
-   **Columns:** 11

| Variable | Definition                                                                     |
|------------------------------------|------------------------------------|
| clade    | Major squamate clade (e.g., Amphisbaenia, Anguiformes, Caenophidia, etc.)      |
| clock    | Clock data type: "Morpho" for morphological clock or "Mol" for molecular clock |
| n        | Number of taxa in that clade                                                   |
| mean     | Mean normalized evolutionary rate                                              |
| sd       | Standard deviation of normalized rates                                         |
| min      | Minimum normalized rate                                                        |
| Q1       | First quartile (25th percentile) of normalized rates                           |
| median   | Median normalized rate                                                         |
| Q3       | Third quartile (75th percentile) of normalized rates                           |
| max      | Maximum normalized rate                                                        |
| range    | Range of normalized rates (max - min)                                          |

------------------------------------------------------------------------

#### Sheet: S39_Clock_4Parts (164t) (Table S39)

**Description:** Summary statistics for normalized evolutionary rates from the 4-partition morphological clock analysis using the full taxon sampling (164 taxa). Rates are broken down by anatomical partition.

-   **Rows:** 41 (1 header + 40 partition-by-clade entries)
-   **Columns:** 12

| Variable  | Definition                                                                                             |
|------------------------------------|------------------------------------|
| partition | Anatomical partition of the morphological data (e.g., Skull_Surface, Mandible, Postcranium, Dentition) |
| clade     | Major squamate clade                                                                                   |
| clock     | Clock partition number identifier                                                                      |
| n         | Number of taxa in that clade                                                                           |
| mean      | Mean normalized evolutionary rate                                                                      |
| sd        | Standard deviation of normalized rates                                                                 |
| min       | Minimum normalized rate                                                                                |
| Q1        | First quartile of normalized rates                                                                     |
| median    | Median normalized rate                                                                                 |
| Q3        | Third quartile of normalized rates                                                                     |
| max       | Maximum normalized rate                                                                                |
| range     | Range of normalized rates                                                                              |

------------------------------------------------------------------------

#### Sheet: S40_Clock_4Parts (GMM_Trees) (Table S40)

**Description:** Summary statistics for normalized evolutionary rates from the 4-partition morphological clock analysis, compared across different GMM tree sizes (taxon samplings with varying fossil representation).

-   **Rows:** 121 (1 header + 120 partition-by-clade-by-tree entries)
-   **Columns:** 12

| Variable  | Definition                                                                                                                                      |
|------------------------------------|------------------------------------|
| partition | Anatomical partition of the morphological data                                                                                                  |
| clade     | Major squamate clade                                                                                                                            |
| n         | Number of taxa in that clade                                                                                                                    |
| mean      | Mean normalized evolutionary rate                                                                                                               |
| sd        | Standard deviation of normalized rates                                                                                                          |
| min       | Minimum normalized rate                                                                                                                         |
| Q1        | First quartile of normalized rates                                                                                                              |
| median    | Median normalized rate                                                                                                                          |
| Q3        | Third quartile of normalized rates                                                                                                              |
| max       | Maximum normalized rate                                                                                                                         |
| range     | Range of normalized rates                                                                                                                       |
| tree_size | Taxon sampling identifier (e.g., "122t" for 122 taxa with 24% fossils, "106t" for 106 taxa with 14% fossils, "92t" for 92 taxa with 0% fossils) |

------------------------------------------------------------------------

#### Sheet: S41_AllMethods_phyloPC (Table S41)

**Description:** Summary statistics for normalized evolutionary rates calculated using phylogenetic principal components across all three rate methods (BayesTraits, RRphylo, and phylodynamic clock), compared across different tree sizes.

-   **Rows:** 91 (1 header + 90 clade-by-tree-by-method entries)
-   **Columns:** 12

| Variable  | Definition                                                                                       |
|------------------------------------|------------------------------------|
| clade     | Major squamate clade                                                                             |
| n         | Number of taxa in that clade                                                                     |
| mean      | Mean normalized evolutionary rate                                                                |
| sd        | Standard deviation of normalized rates                                                           |
| min       | Minimum normalized rate                                                                          |
| Q1        | First quartile of normalized rates                                                               |
| median    | Median normalized rate                                                                           |
| Q3        | Third quartile of normalized rates                                                               |
| max       | Maximum normalized rate                                                                          |
| range     | Range of normalized rates                                                                        |
| tree_size | Taxon sampling identifier (e.g., "106t\_(14%*fossils)","122t*(24%*fossils)","92t*(0%\_fossils)") |
| analysis  | Analytical method used to calculate rates: "BayesTraits", "RRphylo", or "Clock"                  |

------------------------------------------------------------------------

#### Sheet: S42_AllMethods_standPC (Table S42)

**Description:** Summary statistics for normalized evolutionary rates calculated using standard (non-phylogenetic) principal components across all three rate methods, compared across different tree sizes. Same structure and variable definitions as Table S41.

-   **Rows:** 91 (1 header + 90 clade-by-tree-by-method entries)
-   **Columns:** 12

| Variable  | Definition                                                   |
|-----------|--------------------------------------------------------------|
| clade     | Major squamate clade                                         |
| n         | Number of taxa in that clade                                 |
| mean      | Mean normalized evolutionary rate                            |
| sd        | Standard deviation of normalized rates                       |
| min       | Minimum normalized rate                                      |
| Q1        | First quartile of normalized rates                           |
| median    | Median normalized rate                                       |
| Q3        | Third quartile of normalized rates                           |
| max       | Maximum normalized rate                                      |
| range     | Range of normalized rates                                    |
| tree_size | Taxon sampling identifier                                    |
| analysis  | Analytical method used: "BayesTraits", "RRphylo", or "Clock" |

------------------------------------------------------------------------

#### Sheet: S43_AllMethods_Tree_phyloPC (Table S43)

**Description:** Summary statistics of rate metrics aggregated across all clades for each tree size, using phylogenetic PCs. Provides an overview of how evolutionary rate distributions shift across different taxon samplings for all methods combined.

-   **Rows:** 7 (1 header + 6 metric-by-tree entries)
-   **Columns:** 10

| Variable  | Definition                                                                                       |
|------------------------------------|------------------------------------|
| metric    | The summary statistic being aggregated across clades (e.g., "mean", "median")                    |
| tree_size | Taxon sampling identifier (e.g., "106t\_(14%*fossils)","122t*(24%*fossils)","92t*(0%\_fossils)") |
| n         | Number of clade-level observations aggregated                                                    |
| mean      | Mean of the metric across clades                                                                 |
| sd        | Standard deviation of the metric across clades                                                   |
| min       | Minimum value of the metric across clades                                                        |
| Q1        | First quartile of the metric across clades                                                       |
| median    | Median of the metric across clades                                                               |
| Q3        | Third quartile of the metric across clades                                                       |
| max       | Maximum value of the metric across clades                                                        |

------------------------------------------------------------------------

#### Sheet: S44_AllMethods_Tree_standPC (Table S44)

**Description:** Summary statistics of rate metrics aggregated across all clades for each tree size, using standard (non-phylogenetic) PCs. Same structure and variable definitions as Table S43.

-   **Rows:** 7 (1 header + 6 metric-by-tree entries)
-   **Columns:** 10

| Variable  | Definition                                           |
|-----------|------------------------------------------------------|
| metric    | The summary statistic being aggregated across clades |
| tree_size | Taxon sampling identifier                            |
| n         | Number of clade-level observations aggregated        |
| mean      | Mean of the metric across clades                     |
| sd        | Standard deviation of the metric across clades       |
| min       | Minimum value of the metric across clades            |
| Q1        | First quartile of the metric across clades           |
| median    | Median of the metric across clades                   |
| Q3        | Third quartile of the metric across clades           |
| max       | Maximum value of the metric across clades            |

------------------------------------------------------------------------

### Files and variables (other repository contents)

#### File: README.md

**Description:** This README file

#### File: S1_Analyses.zip

**Description:** All analytical input and output files (not included in GitHub repository due to size constraints; available in the Dryad repository)

#### File: S2_R\_codes.zip

**Description:** All R scripts used in this study

#### File: S3_Suppl_Data_Tables.zip

**Description:** All supplementary data tables (S1-S44) in Excel format, as described in detail above

## Code/software 

R v.4.1.1

## Access information

Other publicly accessible locations of the data:

-   <https://github.com/tiago-simoes/BTprocessR2>
-   <https://github.com/tiago-simoes/EvoPhylo>.
