# EpiTrack pipeline

## General information


This suite was developped to facilitate the development of T-Cell-centric SARS-CoV-2 vaccines by enabling the geo-temporal monitoring of the diversification of SARS-CoV-2 T cell epitopes. This is done by levaraging the vast array of SARS-CoV-2-specific genomic data provided by GISAID. Together, the tools available here can:
-	Temporally monitor the diversification of selected T cell epitopes and assess the prevalence of top alternative epitopes 
-	Perform an in-depth analysis of the top lineages/Variants Of Concerns (VOCs) responsible for the diversification of selected epitopes
-	enable the geo-temporal visualization of peptide diversification.


#### External dependencies:
-	pandas 
-	seaborn & matplotlib
-	numpy
-	geopandas
-   sys
-   os

#### Required dataset download:
This toolkit is built around the GISAID SARS-COV-2 genomic alignment file. In order to run the script, users must download the codon-based alignment of sequences from GISAID (named MSAcodonXXXX, where XXXX refers to the build ID) as well as its corresponding metadata, and merge them. Entries with incomplete year/month dates, as well as entries with non-human hosts should be removed. While many informations are available with the metadata, the metadata information crucial to the tools described below are the date, continent, country, and pango-lineage.


## Instructions:
### General instructions:
Prior to running the tools described below, users should modify the indicated lines of the masterscipt, EpiTrack.sh.
Specifically, the following variables should be provided:

-  FOLDER: full path to the directory containing all EpiTrack files
-  GISAID_LATEST_BUILD_FOLDER: full path to the directory where the GISAID files should be stored
-  GISAID_LATEST_BUILD_FILENAME: name of the file featuring the GISAID msa + metadate
-  ORIGINAL_FILE: full path to the directory where the peptide specific files are stored (following processing by the Generating_Peptide_Files.sh script)
-  WORKING_DIRCT: full path to the directory where the EpiTrack suite is being run

### ExtractPeptide_annotated.sh
This function will access the complete GISAID dataset (full MSA and metadata), extract MSA positions corresponding to the provided peptides of interest, and store the peptide-specific MSA + metadata files in folders named after the peptides of interest. peptide-specific folders will be stored within a folder named ExtractedPeptide, found within the same directory where the full GISAID msa + metadata file is stored.

Preset Inputs (set in EpiTrack.sh script): 
-	FOLDER: full path to folder containing EpiTrack scripts
-   GISAID_LATEST_BUILD_FOLDER: full path to directory where latest GISAID build is stored
-   GISAID_LATEST_BUILD_FILENAME: name of file with complete GISAID msa + metadata dataset

commandline inputs:
-   peptide list, in amino acids. Ex: KLPDDFTGC TLNDLNETL NAPRITFGGP VPYNMRVI...

### Alternative_peptide_tracker

This function temporally monitors the diversification of selected T cell epitopes and assesses the prevalence of top alternative epitopes.

Preset Inputs (set in EpiTrack.sh script): 
-	ORIGINAL_FILE: full path to epitope-specific GISAID msa and metadata
-	FOLDER: full path to folder containing EpiTrack scripts
-   WORKING_DIRCT: full path to current working directory

commandline inputs:
-	-z | -zoomed: no to vizualize full breadth of alternative and wild-type peptides; yes to view bottom 10% of alternative epitopes
-	-o | --Output_File: Name of folder where results are to be stored. This folder will be saved in the working directory.

Suggested command: 
./EPITRACK.sh --script Alternative_peptide_tracker.sh -z yes -o RESULTS_ALTERNATIVE_PEPTIDES


### Peptide_Lineage_Tracker
This function performs an in-depth analysis of the top lineages/Variants Of Concerns (VOCs) responsible for the diversification of selected epitopes

Preset Inputs (set in EpiTrack.sh script): 
-	ORIGINAL_FILE: full path to epitope-specific GISAID msa and metadata
-	FOLDER: full path to folder containing EpiTrack scripts
-   WORKING_DIRCT: full path to current working directory

commandline inputs:
-	-o | --Output_File: Name of folder where results are to be stored. This folder will be saved in the working directory.

suggested command: 
./EPITRACK.sh --script Alternative_peptide_tracker.sh -o RESULTS_LINEAGE_TRACKING


### Peptide_Map_Generator
This function enable the geo-temporal visualization of peptide diversification.

Preset Inputs (set in EpiTrack.sh script): 
-	ORIGINAL_FILE: full path to epitope-specific GISAID msa and metadata
-	FOLDER: full path to folder containing EpiTrack scripts
-   WORKING_DIRCT: full path to current working directory

commandline inputs:
-   -g | --Geography: Europe or Worldwide
-   -s | --Peptide_Specific: no to analyse all alternative epitopes for every single epitope given; yes to analyze specific alternative epitopes of interest
-	-o | --Output_File: Name of folder where results are to be stored. This folder will be saved in the working directory.

suggested command:
./EPITRACK.sh --script Peptide_Map_Generator.sh -g Europe -s yes -o RESULTS_MAP


