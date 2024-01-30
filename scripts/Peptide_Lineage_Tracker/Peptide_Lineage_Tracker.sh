#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --account=def-carone
#SBATCH --job-name=19_Compute_Allresults_NUC.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-user=david.hamelin.1@umontreal.ca
#SBATCH --mail-type=FAIL


ORIGINAL_FILES=$1 #$ARG1 #$2 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/GISAID/2023_10_24/PeptideExtracted/"

WORKING_DIRCT=$2 #$ARG2 #$3 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_4_Lineages/" #Create working directory 

#mkdir $WORKING_DIRCT$Output_directory #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_4_Lineages/Lineage_summaries_4_TEST/"

Output_directory=$3 #$ARG3 #$4 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_4_Lineages/Lineage_summaries_4_TEST/" #subdirectory of working directory where result files are sted
mkdir $WORKING_DIRCT"/"$Output_directory

Script_Folder=$4 #$ARG4

cd $ORIGINAL_FILES #Go to raph folder where original files are stored      ####AADLDDFSKQLQ

for FOLDER in *; do #Iterate through folders (for each peptides)
    echo $FOLDER
    mkdir $WORKING_DIRCT"/"$FOLDER #Create a folder in the working directory with the name of the peptide

    cd $FOLDER #enter the original folder
	
    for file in *.data; do  #Iterate through files of the folder
        echo $file
        cat $file | cut -f20 | sort | uniq -c > $WORKING_DIRCT"/"$FOLDER"/"$file".csv" #for each file, extract the proper info and store in working directory of folder as csv    #| uniq -c   | tr " " "," #
    done

    cd $WORKING_DIRCT"/"$FOLDER #Go to the peptide folder in wrking directory

    ls *data.csv > FILES.txt #store all files into FILTE.txt for python program to access


    echo "We're here"
    

    python3 $Script_Folder/scripts/Peptide_Lineage_Tracker/PEPTIDE_LINEAGES_CLEANED_UP.py FILES.txt $WORKING_DIRCT"/"$Output_directory"/"$FOLDER $FOLDER #Right_Format_Events_binding_affinity_immunity_Combine_All_Alternative_peptides.py

    cd $ORIGINAL_FILES

done

#python3 V2f_FullFunction_PAIRWISE_SPARSEMATRIX_SIMPLIFIED_TEST_Multiprocess_With_Network_Make_Figures_With_Frequencies_AllVariants.py 25
