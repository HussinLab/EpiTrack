#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --account=def-carone
#SBATCH --job-name=19_Compute_Allresults_NUC.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-user=david.hamelin.1@umontreal.ca
#SBATCH --mail-type=FAIL


ORIGINAL_FILES="/home/dhamelin/scratch/GISAID_PEPTIDE_DIVERSIFICATION_ANALYSES/RAPH_FILES/MHC_datafiles/"

WORKING_DIRCT="/home/dhamelin/scratch/GISAID_PEPTIDE_DIVERSIFICATION_ANALYSES/ANALYSES/MHCVal_3_MAPS/" #MHCVal_2  Create working directory 

mkdir "/home/dhamelin/scratch/GISAID_PEPTIDE_DIVERSIFICATION_ANALYSES/ANALYSES/MHCVal_3_MAPS/RESULTS_MAPS_PEPTIDE_Specific_GEOPANDAS_CLEANED_UP_TEST/" #Europe

Output_directory="/home/dhamelin/scratch/GISAID_PEPTIDE_DIVERSIFICATION_ANALYSES/ANALYSES/MHCVal_3_MAPS/RESULTS_MAPS_PEPTIDE_Specific_GEOPANDAS_CLEANED_UP_TEST/" #Europe  #subdirectory of working directory where result files are sted

SCRIPT="NOT_Peptide_Specific" #Peptide_Specific

cd $ORIGINAL_FILES #Go to raph folder where original files are stored      ####AADLDDFSKQLQ

for FOLDER in *; do #Iterate through folders (for each peptides)

	#######mkdir $WORKING_DIRCT$FOLDER #Create a folder in the working directory with the name of the peptide

	########cd $FOLDER #enter the original folder
	
	########for file in *.data; do  #Iterate through files of the folder

	    #########cat $file | cut -f3,4,12 | sort | uniq -c > $WORKING_DIRCT$FOLDER"/"$file".csv" #for each file, extract the proper info and store in working directory of folder as csv    #| uniq -c   | tr " " "," #
	##########done

	cd $WORKING_DIRCT$FOLDER #Go to the peptide folder in wrking directory

	#########ls *data.csv > FILES.txt #store all files into FILTE.txt for python program to access


	echo "We're here"
    

    python3 ../Generate_WorldMaps_Peptide_Specific_CLEANED_UP.py FILES.txt $Output_directory$FOLDER $FOLDER $Output_directory $SCRIPT

    cd $ORIGINAL_FILES

done

#python3 V2f_FullFunction_PAIRWISE_SPARSEMATRIX_SIMPLIFIED_TEST_Multiprocess_With_Network_Make_Figures_With_Frequencies_AllVariants.py 25
