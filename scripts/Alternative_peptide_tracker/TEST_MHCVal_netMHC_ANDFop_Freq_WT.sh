#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --account=def-carone
#SBATCH --job-name=19_Compute_Allresults_NUC.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-user=david.hamelin.1@umontreal.ca
#SBATCH --mail-type=FAIL


ORIGINAL_FILES="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/GISAID/2023_10_24/PeptideExtracted/"    ###2023_08_16/MHC_datafiles/ # MHC_nulldatafiles

WORKING_DIRCT="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/" #Create working directory 

mkdir "/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/perMonth_Tables/" #_controls

TEMP_directory="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/perMonth_Tables/" #_controls  subdirectory of working directory where result files are sted


mkdir "/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/NEWVACCINE_NETMHCPAN_AND_TOP_FREQS_WT_Zoomed_REDUCED_CODE_TEST/"

Output_directory="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/NEWVACCINE_NETMHCPAN_AND_TOP_FREQS_WT_Zoomed_REDUCED_CODE_TEST/" #subdirectory of working directory where result files are sted


cd $ORIGINAL_FILES #Go to raph folder where original files are stored      ####AADLDDFSKQLQ

for FOLDER in *; do #Iterate through folders (for each peptides)

	#####mkdir $WORKING_DIRCT$FOLDER #Create a folder in the working directory with the name of the peptide

	###cd $FOLDER #enter the original folder

	cp ./$FOLDER/permonth.tab $TEMP_directory"/"$FOLDER"_permonth.tab"

	######cd $ORIGINAL_FILES

done

rm $Output_directory"AllPeptides_data.csv"

cd $TEMP_directory #Go to the peptide folder in wrking directory

ls *.tab > FILES.txt #store all files into FILTE.txt for python program to access


echo "We're here"


python3 ../GISAID_V4_CLEANED_UP.py FILES.txt $Output_directory $TEMP_directory  #$FOLDER $FOLDER $Output_directory #Right_Format_Events_binding_affinity_immunity_Combine_All_Alternative_peptides.py

    



#python3 V2f_FullFunction_PAIRWISE_SPARSEMATRIX_SIMPLIFIED_TEST_Multiprocess_With_Network_Make_Figures_With_Frequencies_AllVariants.py 25
