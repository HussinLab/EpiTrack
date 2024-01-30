#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --account=def-carone
#SBATCH --job-name=19_Compute_Allresults_NUC.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-user=david.hamelin.1@umontreal.ca
#SBATCH --mail-type=FAIL

#declare -a arr=(KLPDDFTGC,"TLNDLNETL","NAPRITFGGP","VPYNMRVI","RANNTKGSL","GPMVLRGLIT","STTTNIVTR","TGSNVFQTR","HTTDPSFLGR","KTIQPRVEK","TTDPSFLGRYM","PTDNYITTY","YLFDESGEFKL","LPKEITVAT","TTDPSFLGRY","RTIKVFTTV")

#for i in "${arr[@]}"

for i in KLPDDFTGC TLNDLNETL NAPRITFGGP VPYNMRVI RANNTKGSL GPMVLRGLIT STTTNIVTR TGSNVFQTR HTTDPSFLGR KTIQPRVEK TTDPSFLGRYM PTDNYITTY YLFDESGEFKL LPKEITVAT TTDPSFLGRY RTIKVFTTV
do
   echo "$i"
   # or do whatever with individual element of the array
   /lustre06/project/6065672/shared/covid-19/GISAID/code/cov19-MSA-tools/scripts/ExtractPeptide.sh $i
done

