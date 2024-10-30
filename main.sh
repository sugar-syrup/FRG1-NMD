#!/bin/bash

# This script is used to run the RNA-Seq Analysis for the FRG1-NMD project.
# For any queries please contact: shithij.t@niser.ac.in
# Usage:
# ./main.sh -f config.txt <Flag>
# -S : Salmon
# -T : Stringtie
# -I : IsoformSwitchAnalyzeR

if [ $# -eq 0 ]; then
  echo "Oops! You forgot to add options."
  echo "Hint: run $0 -h for help!"
  echo ""
  exit 0
fi

ECHO="false"

while getopts "STIf:" OPTION; do
  case $OPTION in
    T)
      Stringtie_enable="true"
      ;;
    S)
      Salmon_enable="true"
      ;;
    I)
      ISAR_enable="true"
      ;;
    f)
      config_file=$OPTARG
      if [ -e "${config_file}" ]; then
        # Read in design file and variables therein
        source $config_file
        # Check if variables are defined
        if [[ -z $reference || -z $myname || -z $mydir || -z $files_file ]]; then
          echo "One or more variables are undefined"
        else
          ECHO="true"
        fi
      else
        echo "Design file ${config_file} does not exist"
      fi
      ;;
    h)
      echo "Usage: ./main.sh -f config.txt <Flag>"
      echo "Flags:"
      echo "  -S : Salmon"
      echo "  -T : Stringtie"
      echo "  -I : IsoformSwitchAnalyzeR"
      exit 0
      ;;
  esac
done

# Proceed if all variables are defined
if [ $ECHO = "true" ]
then

# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${mydir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/Logs/"
  mkdir ${mydir}/Logs/
fi
log_file="${mydir}/Logs/log_$date"
exec &> >(tee -a "$log_file")

echo "*******************************************************" 
echo " Welcome to the RNA-Seq Analysis Script " 
echo "*******************************************************" 
echo " This script will guide you through the RNA-Seq analysis for the FRG1-NMD project. " 
echo " Please ensure that you have provided all necessary configurations. " 
echo " For assistance, contact: shithij.t@niser.ac.in " 
echo "*******************************************************"
echo "Reference used: $reference"
echo "Name of the study: $myname"
echo "Location of project folder: $mydir"
echo "Experiment file used for this analysis: $files_file"
echo ""

declare -a samples

# Load experiment file, parse for samples and save them into the array
let p=0
while read -r f1 f2; do
  samples[p]="${f1}"
  ((++p))
done < $files_file

# Declare condition arrays
declare -a cond

# Load experiment file, parse for conditions and save unique conditions into array. Note: Condition1 is always "control"
let i=1
cond[0]="control"
while read -r f1 f2; do
  if [[ " ${cond[*]} " == *"$f2"* ]];then
  continue
else
  cond[i]="${f2}"
  ((++i))
fi
done < $files_file

# Declare individual condition arrays
arr_length="$((${#cond[@]}-1))"
for i in $( eval echo {0..${arr_length}} );do
  declare -a cond${i}
done

# Load experiment file again, parse for conditions and save filenames into condition-specific arrays.
while read -r f1 f2; do
  for i in $( eval echo {0..${arr_length}} );do
  if [[ "$f2" == "${cond[i]}" ]];then
    eval cond${i}[cond${i}count]="${f1}"
    ((++cond${i}count))
  fi
done
done < $files_file

if  [ "$Salmon_enable" = "true" ]; then
echo "------------------------------------------"
echo "           Starting Salmon Analysis       "
echo "------------------------------------------"
echo ""
echo -n "Salmon version: "
salmon -v
echo ""

echo "Entering NMD environment"
source /home/shithij/miniconda3/bin/activate NMD

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Start Salmon analysis ${i}"

  if [ ! -d "${mydir}/Salmon/" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Make directory ${mydir}/Salmon/"
    mkdir ${mydir}/Salmon/
  fi

  if [ -d "$mydir/Salmon/${i}" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Sample ${i} already analyzed"
  else
    salmon quant -i /home/shithij/projects/SMG5-SMG7/reference/Transcriptome/$reference -l ISF \
    -1 $mydir/FASTQ/${i}1_001.fastq.gz \
    -2 $mydir/FASTQ/${i}2_001.fastq.gz \
    -p 15 \
    --useVBOpt \
    --gcBias \
    --seqBias \
    --validateMappings \
    -o $mydir/Salmon/${i}
fi
done

echo "Exiting NMD environment"
conda deactivate

echo "------------------------------------------"
echo "           Salmon Analysis Completed      "
echo "------------------------------------------"
echo "Results can be found in the Salmon directory."
echo "See Logs for more details."

fi

if  [ "$ISAR_enable" = "true" ]; then


echo "------------------------------------------"
echo "             Starting ISAR Analysis      "
echo "------------------------------------------"

echo "zentering ISAR environment"
source /home/shithij/miniconda3/bin/activate ISAR

# Generate ISAR folder if necessary
if [ ! -d "${mydir}/ISAR/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/ISAR/"
  mkdir ${mydir}/ISAR/
fi

# Organize conditions from array in a single string
cond_string=$(printf "%s," "${cond[@]}" | cut -d "," -f 1-${#cond[@]})

# Run self-made ISAR script using mydir and the cond_string as positional arguments
$mydir/ISAR_Script.R ${mydir} ${cond_string}

echo "Exiting ISAR environment"
conda deactivate

echo "------------------------------------------"
echo "            ISAR Analysis Completed      "
echo "------------------------------------------"


