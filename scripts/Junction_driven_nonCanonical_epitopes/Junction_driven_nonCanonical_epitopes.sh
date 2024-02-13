#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --account=def-carone
#SBATCH --job-name=19_Compute_Allresults_NUC.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-user=david.hamelin.1@umontreal.ca
#SBATCH --mail-type=FAIL

#Translate_FindPeptideCombined_Strands.py

ORIGINAL_FILES=$1 #$ARG1 #$2 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/GISAID/2023_10_24/PeptideExtracted/"

WORKING_DIRCT=$2 #$ARG2 #$3 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_4_Lineages/" #Create working directory 

#mkdir $WORKING_DIRCT$Output_directory #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_4_Lineages/Lineage_summaries_4_TEST/"

Output_directory=$3 #$ARG3 #$4 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_4_Lineages/Lineage_summaries_4_TEST/" #subdirectory of working directory where result files are sted
mkdir $WORKING_DIRCT"/"$Output_directory

Script_Folder=$4 #$ARG4


python3 $Script_Folder/scripts/Junction_driven_nonCanonical_epitopes/Junction_driven_nonCanonical_epitopes.py $Script_Folder/scripts/Junction_driven_nonCanonical_epitopes/All_Data_Multiprocessed_Three_waves_DelAbove1_LongInterval.csv $WORKING_DIRCT $Script_Folder/scripts/Junction_driven_nonCanonical_epitopes 30




python3 $Script_Folder/scripts/Junction_driven_nonCanonical_epitopes/Print_Figures.py NOTRANSLATE_FrameShift_results_First_Three_waves_Combined_Strands_DelAbove3.csv $WORKING_DIRCT"/"$Output_directory

