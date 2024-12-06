#!/bin/bash

# We have used https://github.com/boehmv/SMG5-SMG7/tree/main as a reference to write the script
# It is a complete pipeline for RNA-seq analysis, from which we have used the Salmon and ISAR analysis part
# We have modified the script to suit our needs

if [ $# -eq 0 ]
then
  echo "Missing options!"
  echo "(run $0 -h for help)"
  echo ""
  exit 0
fi

ECHO="false"

while getopts "ahMCDILiQf:" OPTION; do
  case $OPTION in

    a)
    Salmon_enable="true"
    ISAR_enable="true"
    ;;

    C)
    Salmon_enable="true"
    ;;

    I)
    ISAR_enable="true"
    ;;

    Q)
    QC_enable="true"
    ;;

    f)
    design_file=$OPTARG
    if [ -e "${design_file}" ]; then

	# Read in design file and variables therein
	source $design_file

	# Check if variables are defined
	if [[ -z $seq_design || -z $myname || -z $mydir || -z $experiment_file ]]; then
		echo "One or more variable are undefined"
	else
		ECHO="true"
	fi
    else
	echo "Design file ${design_file} does not exist"
    fi    
    ;;

    h)
    echo "Usage:"
    echo "runAnalysis.sh -f /path/to/design.txt [-OPTION] "
    echo ""
    echo "   -f path/to/design.txt	MANDATORY: give this path to the design file; please see example design file for guidance"
    echo "   -a     to execute the full analysis (all other options enabled)" 
    echo "   -h     help (this output)"
    echo "   -C     to execute Salmon"
    echo "   -I     to execute ISAR"
    echo "   -Q     to execute QC"
    exit 0
    ;;

  esac
done

####################################
#
# Start of script
#
####################################

# Only proceed if acceptable options were given (NOTE: this is one big 'if' command until the end of the script!)

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


# Show design parameters in log
echo "######################"
echo "## General settings ##"
echo "######################"
echo ""
echo "Sequencing design (paired/single): $seq_design"
echo "Name of the study: $myname"
echo "Location of files on server: $srvdir"
echo "Location of project folder: $mydir"
echo "Experiment file used for this analysis: $experiment_file"
echo ""

# Declare samples array
declare -a samples

# Load experiment file, parse for samples and save them into the array
let p=0
while read -r f1 f2; do
  samples[p]="${f1}"
  ((++p))
done < $experiment_file

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
done < $experiment_file

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
done < $experiment_file

# State the conditions and samples for this analysis
echo "############################"
echo "## Conditions and samples ##"
echo "############################"
echo ""

arr_length="$((${#cond[@]}-1))"
for i in $( eval echo {0..${arr_length}} );do
  echo -e "cond${i} \t ${cond[i]} \t $(eval echo \${cond$i[*]})"
done

#Generate folders if necessary

# Create complete sample text file
> ${mydir}/Samples.txt
echo -e "sample \t condition" >> ${mydir}/Samples.txt
cat $experiment_file >> ${mydir}/Samples.txt

# Return state of analysis options
echo ""
echo "######################"
echo "## Analysis options ##"
echo "######################"
echo ""
echo "STAR_enable (-M) = ${STAR_enable}"
echo "Salmon_enable (-C) = ${Salmon_enable}"
echo "DESeq2_enable (-D) = ${DESeq2_enable}"
echo "ISAR_enable (-I) = ${ISAR_enable}"
echo "leafcutter_enable (-L) = ${leafcutter_enable}"
echo "IRFinder_enable (-i) = ${IRFinder_enable}"
echo "QC_enable (-Q) = ${QC_enable}"
echo ""

######
######################
#
# Perform salmon analysis
#
######################
######

if  [ "$Salmon_enable" = "true" ]; then

source /home/shithij/miniconda3/bin/activate NMD

echo ""
echo "###################"
echo "## Perform salmon analysis"
echo "###################"
echo ""
echo -n "Salmon version: "
$mydir/salmon-1.7.0_linux_x86_64/bin/salmon -v
echo ""


# Iterate over fastq files and quantify with Salmon

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

    # Paired end reads
   if  [ "$seq_design" = "paired" ]; then
    $mydir/salmon-1.7.0_linux_x86_64/bin/salmon quant -i /home/shithij/projects/FRG1-NMD/salmon_index -l ISF \
    -1 $mydir/FASTQ/${i}1_001.fastq.gz \
    -2 $mydir/FASTQ/${i}2_001.fastq.gz \
    -p 15 \
    --useVBOpt \
    --gcBias \
    --seqBias \
    --validateMappings \
    -o $mydir/Salmon/${i}

  # Single end reads
  elif  [ "$seq_design" = "single" ]; then
    $mydir/salmon-1.7.0_linux_x86_64/bin/salmon quant -i /home/shithij/projects/FRG1-NMD/salmon_index -l A \
    -r $mydir/FASTQ/${i}.fastq.gz \
    -p 15 \
    --useVBOpt \
    --gcBias \
    --seqBias \
    --validateMappings \
    -o $mydir/Salmon/${i}
  fi
fi
done

conda deactivate

fi

######
####################################
#
# Perform Isoform Switch AnalyseR (ISAR) analysis
#
####################################
######

if  [ "$ISAR_enable" = "true" ]; then


echo ""
echo "###################"
echo "## Perform ISAR analysis"
echo "###################"
echo ""
source /home/shithij/miniconda3/bin/activate ISAR
echo "Entered ISAR environment"


# Generate ISAR folder if necessary
if [ ! -d "${mydir}/ISAR/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/ISAR/"
  mkdir ${mydir}/ISAR/
fi

# Organize conditions from array in a single string
cond_string=$(printf "%s," "${cond[@]}" | cut -d "," -f 1-${#cond[@]})

# Run self-made ISAR script using mydir and the cond_string as positional arguments
$mydir/ISAR1.R ${mydir} ${cond_string}

echo "############################################################################################"
echo "ISAR1.R completed running"
echo "############################################################################################"

# Generate ISAR-specific condition table for running the comparison script
> ${mydir}/ISAR/ConditionTable.txt
echo -e "cond\tpath" >> ${mydir}/ISAR/ConditionTable.txt
echo -e "${cond[@]:1}\t${mydir}/ISAR/SwitchList_filt_Analyzed.csv" >> ${mydir}/ISAR/ConditionTable.txt

# Run ISAR Comparison script
$mydir/ISAR2.R ${mydir} ${mydir}/ISAR/ConditionTable.txt
fi

######
####################################
#
# Perform multiQC analysis
#
####################################
######

if  [ "$QC_enable" = "true" ]; then

echo ""
echo "###################"
echo "## Perform MultiQC analysis"
echo "###################"
echo ""

# Generate QC folder if necessary
if [ ! -d "${mydir}/QC/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/QC/"
  mkdir ${mydir}/QC/
fi

# Execute multiQC
multiqc ${mydir}/FASTQ ${mydir}/Salmon ${srvdir}/${myname} -v -f -o ${mydir}/QC/ --sample-names ${mydir}/Samples.txt

fi

fi