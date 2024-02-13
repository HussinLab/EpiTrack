#!/bin/bash

wp=/lustre06/project/6065672/shared/covid-19/GISAID

lastgisaid=$(find $wp | grep PeptideExtracted$ | rev | cut -f2- -d/ | rev | tail -n 1) #Can change to location of 

scripts=$wp/code/cov19-MSA-tools/scripts #Change This to file with script
ref=$scripts/NC_045512.2.fasta #put this is file with script

peptide_ref=$1


for i in 1 2 3; do
	$scripts/translate $(cat $ref | tail -n 1 | cut -c$i-) |
	grep -bo $peptide_ref | cut -f1 -d: | awk -v i=$i '{print $1*3+i}';
done > $lastgisaid/PeptideExtracted/.temp$peptide_ref.lookup

nbfound=$(wc -l $lastgisaid/PeptideExtracted/.temp$peptide_ref.lookup | cut -f1 -d' ')

echo $peptide_ref was found $nbfound time in the reference: > /dev/stderr

if [ $nbfound -ne 1 ]
then
	echo $peptide_ref found $nbfound times in the reference > /dev/stderr
	echo STOP - please ask RAF > /dev/stderr
	exit;
fi


cuttodo=$( (echo $peptide_ref; cat  $lastgisaid/PeptideExtracted/.temp$peptide_ref.lookup) |
	awk 'NR==1{l=length($1)*3}NR==2{printf "%i-%i",$1,$1+l-1}' )

peptide_wp=$lastgisaid/PeptideExtracted/$peptide_ref



if test -d "$peptide_wp"; then
    echo "Directory exists, remove it in order to run this script"  > /dev/stderr
	echo $peptide_wp
	echo STOP - or ask RAF > /dev/stderr
	exit
else
	echo Creating result directory > /dev/stderr
	echo $peptide_wp
fi

mkdir $peptide_wp


#Create the file
echo Extracting the positions $cuttodo from the last GISAID version in : > /dev/stderr
echo $peptide_wp > /dev/stderr
echo  > /dev/stderr
cat $lastgisaid/msaCodon_*_final.data | cut -c$cuttodo,29892- | awk 'NR%1000000==0{printf "%i lines processed\n",NR > "/dev/stderr"}{print}' > $peptide_wp/full.data


#Translate all nucleotides
echo translating all nucleotides sequences in : > /dev/stderr
echo $peptide_wp/$peptide_ref.msu > /dev/stderr
echo  > /dev/stderr
for nuc_seq in $(cat $peptide_wp/full.data | cut -f1 | grep -v [RYSWKMBDHVN] | awk 'NR%3000000==0{printf "%i lines processed\n",NR > "/dev/stderr"}{t[$1]++}END{for(i in t){print t[i]","i}}' ); do
	echo $nuc_seq $($scripts/translate $(echo $nuc_seq| cut -f2 -d,)) | tr ',' ' ';
done > $peptide_wp/$peptide_ref.msu


#Create the file
echo Keep all amino acid sequences present at least 1000 times : > /dev/stderr
echo $peptide_wp/$peptide_ref.msu.best > /dev/stderr
echo  > /dev/stderr
cat $peptide_wp/$peptide_ref.msu | awk '{t[$3]+=$1;seq[$3]=seq[$3]"|^"$2}END{for(i in t){print i,t[i],seq[i]}}' | awk '$2>1000' |sort -k2,2nr > $peptide_wp/$peptide_ref.msu.best; 


pid_to_wait=""
for alt in $(cat $peptide_wp/$peptide_ref.msu.best | awk '{print $1","$3}' ); do
	egrep "$(echo $alt | cut -f2 -d, | sed 's/|//')" $peptide_wp/full.data > $peptide_wp/$(echo $alt | cut -f1 -d, | sed 's/\*/\\*/').data &
	pid_to_wait=$pid_to_wait" "$!
done

for pid in $pid_to_wait; do
    wait $pid
done

for alt in $(ls $peptide_wp | grep data$ | sed 's/.data$//') ; do
    cat $peptide_wp/$alt.data | cut -f3 | sed 's/\(202[0-3]-\)\([0-9]-\)/\10\2/' | awk -v i=$alt  '{print i,substr($1,1,7)}' | grep -v 00 ;
done |
$scripts/2colTable |
awk 'NR==1{for(i=1;i<=NF;i++){if($i=="full"){ifull=i;$i="other"}}print}NR>1{s=0;for(i=2;i<=NF;i++){if(i!=ifull){s+=$i}};$ifull=$ifull-s;print}' |
sort -nk1,1 | tr ' ' '\t' > $peptide_wp/permonth.tab

echo ================= > /dev/stderr
echo ==== DONE ! ===== > /dev/stderr
echo ================= > /dev/stderr
echo > /dev/stderr
echo check $peptide_wp/permonth.tab > /dev/stderr
