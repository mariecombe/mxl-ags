#!/bin/bash

###### USER DEFINED : ######

Listoffiles="namoptions_P4*"
Outputdir="DROUGHT"
nbdays=21

############################

echo
cd ../
#find Files_scripts/create_sensitivity/ -name ${Listoffiles} | xargs -I '{}' mv {} . # copy all namoptions files to the current directory
if [ -e SensitivityRuns.txt ]  #if the SensitivityRuns.txt file exists
then
rm SensitivityRuns.txt
fi
cp namoptions namoptions.saved

for i in ${Listoffiles} # for each namoptions_P400X file
do
for j in $(jot ${nbdays} 1) # loop from j=1 to j=nbdays
do
fileext="$(echo ${i}|cut -d'_' -f2)" #e.g. P4001, P4002, P4003, etc.
ext="$(printf "%3.3d" $j)"     #e.g. 001, 002, 003, etc. (counts of wg) 
if [[ $fileext != "" ]]
then

if [ $ext == "001" ]  # if this is the first day of the run
then
#3- create a new namoptions file with updated output directory name
cp ${i} ${i}wg${ext}
sed -i '.saved' 4s/"OBS04AUG2007"/"SENS_${fileext}wg${ext}"/1 ${i}wg${ext}
rm -f ${i}wg${ext}.saved
fi

# 1- execute a run
echo "Executing Sensitivity run ${fileext}wg${ext}"
echo "===================================" >> SensitivityRuns.txt
echo " " >> SensitivityRuns.txt
cp ${i}wg${ext} namoptions
./expanded_MXL.exe >> SensitivityRuns.txt
echo " " >> SensitivityRuns.txt
echo " " >> SensitivityRuns.txt
echo " " >> SensitivityRuns.txt

#2- retrieve start and final soil moisture of the day
#"sed -n '/outdir/p' namoptions" prints the line of the namoptions file which contains the string "outdir"
#"sed 's/  */ /g'" will replace all multiple spaces with single spaces in that line
#"cut -d' '" will split the line using a single space delimiter, and "-f3" will select the 3rd value
#"tail -1 output_sca" will print the last line of file named output_sca
#rundir="$(sed -n '/outdir/p' namoptions | sed 's/  */ /g' | cut -d' ' -f3)"
startwg="$(sed -n '/wg /p' namoptions | sed 's/  */ /g' | cut -d' ' -f3)"
rundir="SENS_${fileext}wg${ext}"
cd $rundir
dummy="$(tail -1 output_sca | sed 's/  */ /g' | cut -d' ' -f45)" 
endwg="$(printf "%6.4f" $dummy)"
echo $startwg
echo $endwg
cd ../

if [ $ext != "101" ]
then
#3- create a namoptions file with updated soil moisture for the next iteration
next="$(printf "%3.3d" $(($j+1)))"  #e.g. 002, 003, 004 etc.
cp ${i}wg${ext} ${i}wg${next}
#sed -i '.saved' 123s/"$startwg"/"$endwg"/1 ${i}wg${next} 
sed -i '.saved' 124s/"$startwg"/"$endwg"/1 ${i}wg${next}  
sed -i '.saved' 4s/"SENS_${fileext}wg${ext}"/"SENS_${fileext}wg${next}"/1 ${i}wg${next}
rm -f ${i}wg${next}.saved
fi

fi
done
#rm -f $i
done

cp namoptions.saved namoptions
rm namoptions.saved

# Moving all output to a single directory
Listofdirs="SENS_*"
if [ -d "$Outputdir" ]  #if the directory exists
then
echo "Removing old directory $Outputdir and all its content"
rm -r ${Outputdir}
fi
echo "Creating new empty directory $Outputdir"
mkdir ${Outputdir}
echo "Moving output to directory $Outputdir"
mv ${Listofdirs} ${Outputdir}/
rm -f namoptions_P4*wg*
