#!/bin/bash

#Change the following variables to the proper path/to/directories as indicated below:
#  FOLDER: full path to the directory containing all EpiTrack files
#  DATASET_LATEST_BUILD_FOLDER: full path to the directory where the example data file, titled Example_Data.csv, should be stored
#  DATASET_LATEST_BUILD_FILENAME: name of the example data file (Example_Data.csv)
#  ORIGINAL_FILE: full path to the directory where the peptide specific files are stored (following processing by the Generating_Peptide_Files.sh script)
#  WORKING_DIRCT: full path to the directory where the EpiTrack suite is being run

FOLDER="/home/dhamelin/projects/ctb-hussinju/dhamelin/MHCVal_FinalScripts_GITHUB/Mutational_dynamics_NOSBATCH"
DATASET_LATEST_BUILD_FOLDER="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/GISAID/2023_10_24"
DATASET_LATEST_BUILD_FILENAME="VirusSeq_SubSampling10K.data" #"msaCodon_*_final.data"
ORIGINAL_FILES="/home/dhamelin/projects/ctb-hussinju/shared/covid-19/GISAID/2023_10_24/PeptideExtracted_Viruseq_2/" #PeptideExtracted
WORKING_DIRCT="/home/dhamelin/projects/ctb-hussinju/dhamelin/MHCVal_FinalScripts_GITHUB/MUT_DYN_T_VIRUSSEQ_Test_two" #MUT_DYN_T_NOSBATCH


####################################### DO NOT CHANGE CODE BELOW THIS LINE ################################

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    echo ""
    echo "__________________________________________________"
    echo ""
    echo "                   INSTRUCTIONS"
    echo "__________________________________________________"
    echo ""
    echo "-s | --script Alternative_peptide_tracker.sh"
    echo ""
    echo "        inputs:"
    echo "            -z | --Zoomed: ."
    echo "            -o | --Output_File: Name of directory where results files will be stored"
    echo ""
    echo ""
    echo "-s | --script Peptide_Lineage_Tracker.sh"
    echo ""
    echo "        inputs:"
    echo "            -o | --Output_File: Name of directory where results files will be stored"
    echo ""
    echo ""
    echo "-s | --script Peptide_Map_Generator.sh"
    echo ""
    echo "        inputs:"
    echo "            -g | --Geography: netMHCOutput"
    echo "            -p | --Peptide_Specific: no to analyse all alternative epitopes for every single epitope given; If yes, provide sequence of alternative"
    echo "                                     peptide followed by its unmutated counterpart for the command -p; example: "NALRITFGGP,NAPRITFGGP" "
    echo "                                     ****** The Peptide_Specific function will only work if there is more than one alternative peptide ********"
    echo "            -o | --Output_File: Name of directory where results files will be stored"
    echo ""
    echo ""
    echo "-s | --script ExtractPeptide_annotated_Custom_output.sh"   #ExtractPeptide_annotated.sh
    echo ""
    echo "        inputs:"
    echo "            -l | --PeptideList: list of peptides (no commas, in quotation marks) separated by  space. Ex: "KLPDDFTGC TLNDLNETL NAPRITFGGP VPYNMRVI..." "
    echo "            -o | --Output_File: Name of directory where results files will be stored"
    echo ""
    echo "-s | --script Pandemic_Specific_Conservation.sh"
    echo ""
    echo "        inputs:"
    echo "            -l | --PeptideList: list of peptides (no commas, in quotation marks) separated by  space. Ex: "KLPDDFTGC TLNDLNETL NAPRITFGGP VPYNMRVI..." "
    echo "            -o | --Output_File: Name of directory where results files will be stored"
    echo ""
    echo "-s | --script Junction_driven_nonCanonical_epitopes.sh"
    echo ""
    echo ""
    
     
    

    #EXTENSION="$2"
    shift # past argument
    #shift # past value
    ;;
    -s|--script)
    Script="$2"
    shift # past argument
    shift # past value
    ;;

    #########################################Alternative_peptide_tracker
    -z|--Zoomed)
    ZOOMED="$2"
    shift # past argument
    shift # past value
    ;;

    #########################################Peptide_Map_Generator
    -g|--Geography)
    GEOG="$2"
    shift # past argument
    shift # past value
    ;;

    -p|--Peptide_Specific)
    SPECIF="$2"
    shift # past argument
    shift # past value
    ;;

    #########################################Generate_Peptide_Files
    -l|--PeptideList)
    PEPLIST="$2"
    shift # past argument
    shift # past value
    ;;

    -o|--Output_File)
    OUTP="$2"
    shift # past argument
    shift # past value
    ;;

    


    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

#echo "Script  = ${Script}"
#echo "First word     = ${WordOne}"
#echo "Second Word    = ${WordTwo}"
#echo "DEFAULT         = ${DEFAULT}"
#echo "Number files in SEARCH PATH with EXTENSION:" $(ls -1 "${SEARCHPATH}"/*."${EXTENSION}" | wc -l)
#if [[ -n $1 ]]; then
#    echo "Last line of file specified as non-opt/last argument:"
#    tail -1 "$1"
#fi

if [ ! -z "${Script}"  ] ; then




    if [ $Script == "Alternative_peptide_tracker.sh" ]; then
        echo "Running Alternative_peptide_tracker.sh"
        if [ ! -z "${ZOOMED}"  ] && [ ! -z "${OUTP}" ]; then #[ ! -z "${mutFile1}"  ] && [ ! -z "${HLAs}" ] && [ ! -z "${outputFilePATH}" ] && [ ! -z "${netMHCpanPATH}" ]
            echo "We're here"
            echo "${ZOOMED}"

            $FOLDER/scripts/Alternative_peptide_tracker/$Script $ZOOMED $ORIGINAL_FILES $WORKING_DIRCT $OUTP $FOLDER     #sbatch --export=ARG1=$ZOOMED,ARG2=$ORIGINAL_FILES,ARG3=$WORKING_DIRCT,ARG4=$OUTP,ARG5=$FOLDER
        else
            echo "one or more arguments are missing"

        fi



    elif [ $Script == "Peptide_Lineage_Tracker.sh" ]; then
        echo "Running Peptide_Lineage_Tracker.sh"
        if [ ! -z "${OUTP}" ]; then #[ ! -z "${mutFile1}"  ] && [ ! -z "${HLAs}" ] && [ ! -z "${outputFilePATH}" ] && [ ! -z "${netMHCpanPATH}" ]
            echo "We're here"
            echo "${OUTP}"

            $FOLDER/scripts/Peptide_Lineage_Tracker/$Script $ORIGINAL_FILES $WORKING_DIRCT $OUTP $FOLDER   #sbatch --export=ARG1=$ORIGINAL_FILES,ARG2=$WORKING_DIRCT,ARG3=$OUTP,ARG4=$FOLDER

        else
            echo "one or more arguments are missing"

        fi
    


    elif [ $Script == "Peptide_Map_Generator.sh" ]; then
        echo "Peptide_Map_Generator.sh"
        if [ ! -z "${GEOG}" ] && [ ! -z "${OUTP}" ]; then
            echo "We're here"
            
            $FOLDER/scripts/Peptide_Map_Generator/$Script $GEOG $ORIGINAL_FILES $WORKING_DIRCT $OUTP $SPECIF $FOLDER

        else
            echo "one or more arguments are missing"

        fi

    


    elif [ $Script == "ExtractPeptide_annotated_Custom_output.sh" ]; then #ExtractPeptide_annotated.sh
        echo "Running ExtractPeptide_annotated_Custom_output.sh" #ExtractPeptide_annotated.sh
        if [ ! -z "${PEPLIST}" ] && [ ! -z "${OUTP}" ]; then # && [ ! -z "${OUTP}" ]
            echo "We're here"

            for i in $PEPLIST
            do
                echo "$i"
                # or do whatever with individual element of the array
                $FOLDER/scripts/Generate_Peptide_Files/$Script $DATASET_LATEST_BUILD_FOLDER $FOLDER $i $DATASET_LATEST_BUILD_FILENAME $OUTP
            done
            
             

        else
            echo "one or more arguments are missing"

        fi
    

    elif [ $Script == "Pandemic_Specific_Conservation.sh" ]; then
        echo "Running Pandemic_Specific_Conservation.sh"
        if [ ! -z "${OUTP}" ]; then #[ ! -z "${mutFile1}"  ] && [ ! -z "${HLAs}" ] && [ ! -z "${outputFilePATH}" ] && [ ! -z "${netMHCpanPATH}" ]
            echo "We're here"
            echo "${OUTP}"

            $FOLDER/scripts/Pandemic_Specific_Conservation/$Script $ORIGINAL_FILES $WORKING_DIRCT $OUTP $FOLDER   #sbatch --export=ARG1=$ORIGINAL_FILES,ARG2=$WORKING_DIRCT,ARG3=$OUTP,ARG4=$FOLDER

        else
            echo "one or more arguments are missing"

        fi
    

    elif [ $Script == "Junction_driven_nonCanonical_epitopes.sh" ]; then
        echo "Running Junction_driven_nonCanonical_epitopes.sh"
        if [ ! -z "${OUTP}" ]; then #[ ! -z "${mutFile1}"  ] && [ ! -z "${HLAs}" ] && [ ! -z "${outputFilePATH}" ] && [ ! -z "${netMHCpanPATH}" ]
            echo "We're here"
            echo "${OUTP}"

            $FOLDER/scripts/Junction_driven_nonCanonical_epitopes/$Script $ORIGINAL_FILES $WORKING_DIRCT $OUTP $FOLDER   #sbatch --export=ARG1=$ORIGINAL_FILES,ARG2=$WORKING_DIRCT,ARG3=$OUTP,ARG4=$FOLDER

        else
            echo "one or more arguments are missing"

        fi




    else echo "option does not correspond to any script."

    fi

else
    echo ""
    echo ""
    echo "No scripts were provided"

fi










