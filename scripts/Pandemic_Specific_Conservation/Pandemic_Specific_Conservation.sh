#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --account=def-carone
#SBATCH --job-name=19_Compute_Allresults_NUC.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-user=david.hamelin.1@umontreal.ca
#SBATCH --mail-type=FAIL


#ORIGINAL_FILES="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/GISAID/2023_10_24/PeptideExtracted/"    ###2023_08_16/MHC_datafiles/ # MHC_nulldatafiles

#WORKING_DIRCT="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_5_Raph_Fig/" #Create working directory 

#mkdir "/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_5_Raph_Fig/MSU/" #_controls

#TEMP_directory="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_5_Raph_Fig/MSU/" #_controls  subdirectory of working directory where result files are sted


#mkdir "/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_5_Raph_Fig/RAPH_FIG_V1/"

#Output_directory="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_5_Raph_Fig/RAPH_FIG_V1/" #subdirectory of working directory where result files are sted

ORIGINAL_FILES=$1 #$ARG2  #$2 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/GISAID/2023_10_24/PeptideExtracted/"    ###2023_08_16/MHC_datafiles/ # MHC_nulldatafiles
echo $ORIGINAL_FILES

WORKING_DIRCT=$2 #$ARG3 #$3 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/" #Create working directory 
echo $WORKING_DIRCT
#mkdir  #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/perMonth_Tables/" #_controls

TEMP_directory=$WORKING_DIRCT"/MSU/"  #$4 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/perMonth_Tables/" #_controls  subdirectory of working directory where result files are sted
mkdir $TEMP_directory

Script_Folder=$4 #$ARG5 #$5

#mkdir  #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/NEWVACCINE_NETMHCPAN_AND_TOP_FREQS_WT_Zoomed_REDUCED_CODE_TEST/"

Output_directory=$3 #$ARG4 #$4 #"/home/dhamelin/projects/ctb-hussinju/shared/covid-19/Kovalchik_paper/Peptide_Temporal_analysis/MHCVal_2/NEWVACCINE_NETMHCPAN_AND_TOP_FREQS_WT_Zoomed_REDUCED_CODE_TEST/" #subdirectory of working directory where result files are sted
mkdir $WORKING_DIRCT"/"$Output_directory

cd $ORIGINAL_FILES #Go to raph folder where original files are stored      ####AADLDDFSKQLQ

for FOLDER in *; do #Iterate through folders (for each peptides)

	#####mkdir $WORKING_DIRCT$FOLDER #Create a folder in the working directory with the name of the peptide

	###cd $FOLDER #enter the original folder

    #cat ./$FOLDER/*msu | cut -f1 | sort | uniq -c > $TEMP_directory"/"$FOLDER"_MSU.csv"
    cp ./$FOLDER/*msu $TEMP_directory  #"/"$FOLDER"_MSU.csv"
	#cp ./$FOLDER/permonth.tab $TEMP_directory"/"$FOLDER"_permonth.tab"

	######cd $ORIGINAL_FILES

done

#rm $Output_directory"AllPeptides_data.csv"

cd $TEMP_directory #Go to the peptide folder in wrking directory

ls *.msu > FILES.txt #store all files into FILTE.txt for python program to access


#echo "We're here"


python3 $Script_Folder/scripts/Pandemic_Specific_Conservation/Create_RAPH_Fig_input_Cleaned_up.py FILES.txt $WORKING_DIRCT"/"$Output_directory"/" $TEMP_directory $Script_Folder/scripts/Pandemic_Specific_Conservation/  #$FOLDER $FOLDER $Output_directory #Right_Format_Events_binding_affinity_immunity_Combine_All_Alternative_peptides.py

    



#python3 V2f_FullFunction_PAIRWISE_SPARSEMATRIX_SIMPLIFIED_TEST_Multiprocess_With_Network_Make_Figures_With_Frequencies_AllVariants.py 25
