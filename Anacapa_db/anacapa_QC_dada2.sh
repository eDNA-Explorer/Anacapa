#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/anacapa_QC_dada2.sh -i <input_dir> -o <out_dir> -d <database_directory> -f <fasta file of forward primers> -r <fasta file of reverse primers> -u <HPC_user_name>, -g (add flag (-g), no text required, if fastq files are not compressed), -c change cut adapt error 3' and 5' trimming .3 default need value, -p change cut adapt error primer sorting / trimming .3 default need value, -q minimum quality score, -x additional base pairs trimmed from 5' end forward read, -y additional base pairs trimmed from 5' end reverse read, -b number of times an ASV must be found to retain after dada2, -e file path to minimum overlap for user determined F and R primers, -k path to custom HPC job submission header
HELP=""
IN=""
OUT=""
DB=""
FP=""
RP=""
GUNZIPED=""
CTADE=""
PCTADE=""
QUALS=""
MILEN=""
FETRIM=""
RETRIM=""
MINTIMES_ASV=""
PHIX_REF="/usr/local/bin/resources/phix174_ill.ref.fa.gz"
CONTAMINANTS="/usr/local/bin/resources/adapters.fa"

while getopts "h?:i:o:d:f:r:g?:c:p:q:m:x:y:b:t:" opt; do
    case $opt in
        h) HELP="TRUE"
        ;;
        i) IN="$OPTARG" # path to raw .fastq.gz files
        ;;
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        f) FP="$OPTARG"  # need forward reads for cutadapt
        ;;
        r) RP="$OPTARG"  # need reverse reads for cutadapt
        ;;
        g) GUNZIPED="TRUE" #reads not compressed
        ;;
        c) CTADE="$OPTARG"  # cutadapt error for 3' adapter and 5' primer adapter trimming
        ;;
        p) PCTADE="$OPTARG"  # cutadapt error 3' primer sorting and trimming
        ;;
        q) QUALS="$OPTARG"  # Minimum Quality score
        ;;
        m) MILEN="$OPTARG"  # Minimum Length after trimming
        ;;
        x) FETRIM="$OPTARG"  # Additional 5' trimming forward read
        ;;
        y) RETRIM="$OPTARG"  # Additional 5' trimming reverse read
        ;;
        b) MINTIMES_ASV="$OPTARG" # number of occurances required to keep ASV
        ;;
        t) METADATA="$OPTARG"
        ;;
    esac
done

# if $DB has a / at the end remove it!
case "$DB" in
*/)
    DB=${DB%/}
    ;;
*)
    echo ""
    ;;
esac

if [ "${HELP}" = "TRUE" ]
then
  printf "<<< Anacapa: Sequence QC and ASV Parsing >>>\n\nThe purpose of these script is to process raw fastq or fastq.gz files from an Illumina HiSeq or MiSeq.  It removes 3' and 5' sequencing artifacts and 5' metabarcode primers (cutadapt), removes low quality base pairs and short reads (fastX-toolkit), sorts reads by 3' metabarcode primers prior to trimming (cutadapt), and uses dada2 to denoise, dereplicate, merge and remove chimeric reads\n\n	For successful implementation \n		1. Make sure you have all of the dependencies and correct paths in the anacapa_config.sh file\n		2. Add the Metabarcode locus specific CRUX reference libraries to the Anacapa_db folder\n		3. All parameters can be modified using the arguments below.  Alternatively, all parameters can be altered in the anacapa_vars.sh folder\n\nArguments:\n- Required for either mode:\n	-i	path to .fastq.gz files, if files are already compressed use flag -g (see below)\n	-o	path to output directory\n	-d	path to Anacapa_db\n	-t	Illumina Platform: HiSeq (2 x 150) or MiSeq ( >= 2 x 250)\n    \n- Optional:\n 	-u	If running on an HPC (e.g. UCLA's Hoffman2 cluster), this is your username: e.g. eecurd\n 	-f	path to file with forward primers in fasta format \n    		e.g.	>16s\n    			GTGYCAGCMGCCGCGGTAA\n			>18S\n			GTACACACCGCCCGTC\n	-r	path to file with forward primers in fasta format \n    		e.g. 	>16s\n    			GGACTACNVGGGTWTCTAAT\n    			>18S\n			TGATCCTTCTGCAGGTTCACCTAC\n	-g	If .fastq read are not compressed: -g (no argument need)\n	-c	To modify the allowed cutadapt error for 3' adapter and 5' primer adapter trimming: 0.0 to 1.0 (default 0.3)\n	-p	To modify the allowed cutadapt error 3' primer sorting and trimming: 0.0 to 1.0 (default 0.3)\n	-q	To modify the minimum quality score allowed: 0 - 40 (default 35)\n	-m	To modify the minimum length after quality trimming: 0 - 300 (default 100)\n	-x	To modify the additional 5' trimming of forward reads: 0 - 300 (default HiSeq 10, default MiSeq 20)\n	-y	To modify the additional 5' trimming of reverse reads: 0 - 300 (default HiSeq 25, default MiSeq 50)\n	-b	To modify the number of occurrences required to keep an ASV: 0 - any integer (default 0)\n	-e	File path to a list of minimum length(s) reqired for paired F and R reads to overlap \n		(length of the locus - primer length + 20 bp). The user should take into account variability in amplicon \n		region (e.g.The amplicon size for 18S 1389f-1510r is ~260 +/- 50 bp) and make appropriate allowances.\n		e.g.	LENGTH_16S="235"\n			LENGTH_18S="200"\n	-k	Path to file with alternate HPC job submission parameters:  \n		default file = ~/Anacapa_db/scripts/Hoffman2_HPC_header.sh\n		modifiable template file = ~/Anacapa_db/scripts/anacapa_qsub_templates.sh\n\n\n-Other:\n	-h	Shows program usage then quits\n\n\n\n"
  exit
else
  echo ""
fi

####################################script & software
# This pipeline was developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Gaurav Kandlikar (gkandlikar@ucla.edu), and Baochen Shi (biosbc@gmail.com), and with contributions from Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Last Updated 11-18-2017
#
# The purpose of these script is to process raw fastq.gz files from an Illumina sequencing and generate summarized taxonomic assignment tables for multiple metabarcoding targets.
#
# This script is currently designed to run on UCLA's Hoffman2 cluster.  Please adjust the code to work with your computing resources. (e.g. module / path names to programs, submitting jobs for processing if you have a cluster, etc)
#
# This script runs in two phases, the first is a QC and dada2 seqeunce dereplication, denoising, mergeing (if reads are paired) and chimera detection.  The second phase runs bowtie2 and our blowtie2 modified blca run_scripts.
#
######################################

###Check that User has the correct set of primer files
if [[ ! -e ${FP} && ! -e ${RP} ]];
then
  echo "Using Default Primers"
elif [[ -e ${FP} && -e ${RP} ]];
then
  echo "Using User Defined Primers"
elif [[ -e ${FP} && ! -e ${RP} ]];
then
  echo "Using User Defined Primers"
  echo "Missing File For Reverse Primer"
  echo ""
  exit
elif [[ ! -e ${FP} && -e ${RP} ]];
then
  echo "Using User Defined Primers"
  echo "Missing File For Forward Primer!"
  echo ""
  exit
fi
######

### Check that user has a valid metadata file
if [ -f "$METADATA" ] && [[ "$METADATA" == *.csv ]]; then
    echo "Metadata file exists and is a CSV file."
else
    echo "Metadata file does not exist or is not a CSV file."
    exit
fi

# Check that PHIX_REF and CONTAMINANTS files exists
if [ -f "$PHIX_REF" ] && [ -f "$CONTAMINANTS" ]; then
  echo "PHIX_REF and CONTAMINANTS file exist."
else
  echo "PHIX_REF and/or CONTAMINANTS file do not exist."
  exit
fi

#Check that user has all of the default flags set
if [[ -e ${IN} && ! -z ${OUT} && -e ${DB} ]];
then
  echo "Required Arguments Given"
  echo ""
else
  echo "Required Arguments Missing:"
  echo "check that you included arguments or correct paths for -i -o -d -a"
  echo ""
  exit
fi

# location of the config and var files
source $DB/scripts/anacapa_vars.sh  # edit to change variables and parameters
source $DB/scripts/anacapa_config.sh # edit for proper configuration

###

################################
# Preprocessing .fastq files
################################
suffix1=R1_001.fastq.gz
suffix2=R2_001.fastq.gz
###################################
mkdir -p ${OUT}
mkdir -p ${OUT}/Run_info
mkdir -p ${OUT}/Run_info/run_logs
mkdir -p ${OUT}/QC
mkdir -p ${OUT}/QC/fastq

echo " "
date
echo " "
echo "Preprocessing: 1) Generate an md5sum file"  # user can check for file corruption
md5sum ${IN}/*fastq.gz > ${OUT}/Run_info/raw_fastq.md5sum
date
echo "Preprocessing: 2) Change file suffixes"
is_gzipped() {
  local file_path="$1"
  local magic_number=$(hexdump -n 2 -v -e '/1 "%02x"' "$file_path")
  [[ "$magic_number" == "1f8b" ]]
}

# Read the file names into an array
readarray -t filename_pairs < <(awk -F',' '
    BEGIN {OFS=","} 
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i == "Fastq Forward Reads Filename") f_col=i;
            if ($i == "Fastq Reverse Reads Filename") r_col=i;
        }
    }
    NR > 1 {
        # Replace all "-" with "_" in the forward and reverse reads filenames
        gsub("-", "_", $f_col);
        gsub("-", "_", $r_col);

        # Change file extension from "fq.gz" to "fastq.gz"
        sub(/fq\.gz$/, "fastq.gz", $f_col);
        sub(/fq\.gz$/, "fastq.gz", $r_col);

        # Append fastq.gz if it doesnt end with it
        if ($f_col !~ /fastq\.gz$/) $f_col = $f_col ".fastq.gz";
        if ($r_col !~ /fastq\.gz$/) $r_col = $r_col ".fastq.gz";

        print $f_col, $r_col;
    }
' "$METADATA")

declare -a updated_filename_pairs
for pair in "${filename_pairs[@]}"
do
  echo "Processing pair: $pair"
  IFS=',' read -r forward_file reverse_file <<< "$pair"
  # remove file extensions
  forward_file_wo_ext=$forward_file
  reverse_file_wo_ext=$reverse_file
  while [[ "$forward_file_wo_ext" == *.* ]]; do
    forward_file_wo_ext=$(basename "$forward_file_wo_ext" .${forward_file_wo_ext##*.})
  done
  while [[ "$reverse_file_wo_ext" == *.* ]]; do
    reverse_file_wo_ext=$(basename "$reverse_file_wo_ext" .${reverse_file_wo_ext##*.})
  done
  # replace all underscores
  forward_file_mod="${forward_file_wo_ext//_/-}_1.fastq.gz"
  reverse_file_mod="${reverse_file_wo_ext//_/-}_2.fastq.gz"

  # check if files are gzipped
  forward_file_type=$(file -b "${IN}/${forward_file}")
  reverse_file_type=$(file -b "${IN}/${reverse_file}")

  if [[ "${forward_file_type}" == *"gzip"* ]]; then
    # file is gzipped, copy directly
    cp ${IN}/${forward_file} ${OUT}/QC/fastq/${forward_file_mod}
  else
    # file is not gzipped
    # rm ".gz" suffix and gzip the file.
    echo "File ${IN}/${forward_file} is not gzipped. Removing '.gz' suffix and compressing file..."
    mv ${IN}/${forward_file} ${IN}/${forward_file%.gz}
    gzip ${IN}/${forward_file%.gz}
    # copy
    cp ${IN}/${forward_file} ${OUT}/QC/fastq/${forward_file_mod}
    echo "Done."
  fi
  if [[ "${reverse_file_type}" == *"gzip"* ]]; then
    # file is gzipped, copy directly
    cp ${IN}/${reverse_file} ${OUT}/QC/fastq/${reverse_file_mod}
  else
    # file is not gzipped
    # rm ".gz" suffix and gzip the file.
    echo "File ${IN}/${reverse_file} is not gzipped. Removing '.gz' suffix and compressing file..."
    mv ${IN}/${reverse_file} ${IN}/${reverse_file%.gz}
    gzip ${IN}/${reverse_file%.gz}
    # copy
    cp ${IN}/${reverse_file} ${OUT}/QC/fastq/${reverse_file_mod}
    echo "Done."
  fi

  # Update filename_pairs with new locations
  updated_filename_pairs+=("${OUT}/QC/fastq/${forward_file_mod},${OUT}/QC/fastq/${reverse_file_mod}")
done

# Now, replace the original filename_pairs with the updated_filename_pairs
filename_pairs=("${updated_filename_pairs[@]}")
updated_filename_pairs=() # empty array
date

###

################################
# QC the preprocessed .fastq files
#############################

echo "QC: 1) Run cutadapt to remove 5' sequencing adapters and 3' primers + sequencing adapters, sort for length, and quality."

# Generate cut adapt primer files -> merge reverse complemented primers with adapters for cutting 3'end sequencing past the end of the metabarcode region, and add cutadapt specific characters to primers and primer/adapter combos so that the appropriate ends of reads are trimmed
mkdir -p ${OUT}/Run_info/cutadapt_primers_and_adapters

echo " "
echo "Generating Primer and Primer + Adapter files for cutadapt steps."
cp ${DB}/adapters_and_PrimAdapt_rc/*_truseq_*_adapter.txt ${OUT}/Run_info/cutadapt_primers_and_adapters  # make a copy of the appropriate adapter file in your ourput directory
cp ${DB}/adapters_and_PrimAdapt_rc/*_nextera_*_adapter.txt ${OUT}/Run_info/cutadapt_primers_and_adapters  # make a copy of the appropriate adapter file in your ourput directory

cat ${OUT}/Run_info/cutadapt_primers_and_adapters/g_nextera_Forward_adapter.txt >> ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt
sed -i "s/F_adapt/F_adapt_nextera/" ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt
echo "" >> ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt
cat ${OUT}/Run_info/cutadapt_primers_and_adapters/g_truseq_Forward_adapter.txt >> ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt

cat ${OUT}/Run_info/cutadapt_primers_and_adapters/G_nextera_Reverse_adapter.txt >> ${OUT}/Run_info/cutadapt_primers_and_adapters/G_Reverse_adapter.txt
sed -i "s/R_adapt/R_adapt_nextera/" ${OUT}/Run_info/cutadapt_primers_and_adapters/G_Reverse_adapter.txt
echo "" >> ${OUT}/Run_info/cutadapt_primers_and_adapters/G_Reverse_adapter.txt
cat ${OUT}/Run_info/cutadapt_primers_and_adapters/G_truseq_Reverse_adapter.txt >> ${OUT}/Run_info/cutadapt_primers_and_adapters/G_Reverse_adapter.txt

python3 ${DB}/scripts/anacapa_format_primers_cutadapt.py ${FP:=$FP_PATH} ${RP:=$RP_PATH} ${OUT}/Run_info/cutadapt_primers_and_adapters  # given your adapter and primer sets, make cutadapt readable fasta files for trimming adapter / primer reads

# now use the formated cutadapt primer file to trim fastq reads
mkdir -p ${OUT}/QC/fastp_cleaned
mkdir -p ${OUT}/QC/fastp_logs
mkdir -p ${OUT}/QC/bbduk_cleaned
mkdir -p ${OUT}/QC/bbduk_logs
mkdir -p ${OUT}/QC/cutadapt_fastq
mkdir -p ${OUT}/QC/cutadapt_fastq/untrimmed
mkdir -p ${OUT}/QC/cutadapt_fastq/primer_sort
mkdir -p ${OUT}/Run_info/cutadapt_out
###
for pair in "${filename_pairs[@]}"
do
  echo "Processing pair: $pair"
  IFS=',' read -r forward_file reverse_file <<< "$pair"
  j1=${forward_file%_*}
  j1=${j1#${OUT}/QC/fastq/}
  j2=${reverse_file%_*}
  j2=${j2#${OUT}/QC/fastq/}
  echo ${j1} "..."
  echo ${j2} "..."
  
  # fastp processing
  echo "Running fastp for quality and adapter trimming..."
  fastp -i ${forward_file} -I ${reverse_file} \
        -o ${OUT}/QC/fastp_cleaned/${j1}_clean_1.fastq.gz -O ${OUT}/QC/fastp_cleaned/${j2}_clean_2.fastq.gz \
        -h ${OUT}/QC/fastp_logs/fastp_report.html -j ${OUT}/QC/fastp_logs/fastp_report.json \
        --detect_adapter_for_pe \
        --thread 8 \
        --cut_right \
        --cut_right_window_size 4 --cut_right_mean_quality ${QUALS:=$MIN_QUAL} \
        --cut_tail \
        --cut_tail_window_size 4 --cut_tail_mean_quality ${QUALS:=$MIN_QUAL} \
        --low_complexity_filter \
        --complexity_threshold ${QUALS:=$MIN_QUAL} \
        --correction \
        --length_required ${MILEN:=$MIN_LEN} \
        --html \
        --json
  # remove intermediate files
  rm ${forward_file} ${reverse_file}
  
  # bbmap processing for further quality control
  echo "Running bbmap for additional error correction and quality filtering..."
  bbduk.sh in1=${OUT}/QC/fastp_cleaned/${j1}_clean_1.fastq.gz in2=${OUT}/QC/fastp_cleaned/${j2}_clean_2.fastq.gz \
          out1=${OUT}/QC/bbduk_cleaned/${j1}_clean_1.fastq.gz out2=${OUT}/QC/bbduk_cleaned/${j2}_clean_2.fastq.gz \
          ref=${PHIX_REF},${CONTAMINANTS} \
          ktrim=r k=23 mink=11 hdist=1 tpe tbo \
          qtrim=r trimq=30 \
          minlen=50 \
          maxns=1 \
          -Xmx4g stats=${OUT}/QC/bbduk_logs/bbmap_stats.txt
  # remove intermediate files
  rm ${OUT}/QC/fastp_cleaned/${j1}_clean_1.fastq.gz ${OUT}/QC/fastp_cleaned/${j2}_clean_2.fastq.gz
  # chop off the 5' adapter and 3' adapter and primer combo (reverse complemented)
  # this step removes all primers and adapters with the exception of the 5' forward and reverse primers.  These are needed in a later step to sort reads by primer set.  Leaving 3' primers and 5' or 3' adapters can affect read merging and taxonomic assignment.
  # this cutadapt command allows a certain amount of error/missmatch (-e) between the query (seqeuncing read) and the primer and adapter.  It searches for and trims off all of the 5' forward adapter (-g) and the 3' reverse complement reverse primer / reverse complement reverse adapter (-a) or the 5' reverse adapter (-G) and the 3' reverse complement forward primer / reverse complement forward adapter (-A).  It processes read pairs, and results in two files one for each read pair.
  ${CUTADAPT} -e ${CTADE:=$ERROR_QC1} -g ${F_ADAPT} -a ${Rrc_PRIM_ADAPT} -G ${R_ADAPT} -A ${Frc_PRIM_ADAPT} --minimum-length 1 -o ${OUT}/QC/cutadapt_fastq/${j1}_qcPaired_1.fastq -p ${OUT}/QC/cutadapt_fastq/${j2}_qcPaired_2.fastq ${OUT}/QC/bbduk_cleaned/${j1}_clean_1.fastq.gz ${OUT}/QC/bbduk_cleaned/${j2}_clean_2.fastq.gz >> ${OUT}/Run_info/cutadapt_out/cutadapt-report.txt
  rm ${OUT}/QC/bbduk_cleaned/${j1}_clean_1.fastq.gz # remove intermediate files
  rm ${OUT}/QC/bbduk_cleaned/${j2}_clean_2.fastq.gz # remove intermediate files
  # sort by metabarcode but run additional trimming.  It makes a differnce in merging reads in dada2.  Trimming varies based on sequencing platform.
  echo "forward..."
   # use cut adapt to search 5' end of forward reads for forward primers.  These are then sorted by primer name.  We do an additional trimming step analagous to the trimming step in the dada2 tutorial.  Because these a forward reads an tend to be higher quality we only trim  20 bp from the end by default for the MiSeq (longer Reads). Users can modify all parameters in the vars file.
  ${CUTADAPT} -e ${PCTADE:=$ERROR_PS} -g ${F_PRIM} -u -${FETRIM:=$MS_F_TRIM} --minimum-length 1 -o ${OUT}/QC/cutadapt_fastq/primer_sort/{name}_${j1}_Paired_1.fastq  ${OUT}/QC/cutadapt_fastq/${j1}_qcPaired_1.fastq >> ${OUT}/Run_info/cutadapt_out/cutadapt-report.txt
  echo "check"
  echo "reverse..."
  # use cut adapt to search 5' end of reverse reads for reverse primers.  These are then sorted by primer name.  We do an additional trimming step analagous to the trimming step in the dada2 tutorial.  Because these a reverse reads an tend to be lower quality we only trim  50 bp from the end by default for the MiSeq (longer Reads). Users can modify all parameters in the vars file.
  ${CUTADAPT} -e ${PCTADE:=$ERROR_PS} -g ${R_PRIM}  -u -${RETRIM:=$MS_R_TRIM} -o ${OUT}/QC/cutadapt_fastq/primer_sort/{name}_${j2}_Paired_2.fastq   ${OUT}/QC/cutadapt_fastq/${j2}_qcPaired_2.fastq >> ${OUT}/Run_info/cutadapt_out/cutadapt-report.txt
  echo "check"
  rm ${OUT}/QC/cutadapt_fastq/${j1}_qcPaired_1.fastq # remove intermediate files
  rm ${OUT}/QC/cutadapt_fastq/${j2}_qcPaired_2.fastq # remove intermediate files

  # Update filename_pairs with new locations
  metabarcode="$( ls -l | grep -o '>.*' ${FP:=$FP_PATH} | cut -c 2- | tr '\n' ' ' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' )"
  updated_filename_pairs+=("${OUT}/QC/cutadapt_fastq/primer_sort/${metabarcode}_${j1}_Paired_1.fastq,${OUT}/QC/cutadapt_fastq/primer_sort/${metabarcode}_${j2}_Paired_2.fastq")
done
# Now, replace the original filename_pairs with the updated_filename_pairs
filename_pairs=("${updated_filename_pairs[@]}")
updated_filename_pairs=()
date
###

###############################
# Make sure unassembled reads are still paired
###############################
mkdir -p ${OUT}/Run_info
mkdir -p ${OUT}/Run_info/run_scripts

metabarcode="$( ls -l | grep -o '>.*' ${FP:=$FP_PATH} | cut -c 2- | tr '\n' ' ' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' )"


echo ${metabarcode}

echo " "
echo "Checking that Paired reads are still paired:"
for j in ${metabarcode}
do
  echo " "
  echo ${j} "..."
  # make directories for the soon to be sorted reads.  Reads will be sorted by primer and then by paired or unpaired read status.
  shopt -s nullglob
  files=(${OUT}/QC/cutadapt_fastq/primer_sort/${j}*)
  # now check the size of the array
  if (( ${#files[@]} == 0 )); then
      echo "No reads matched the ${j} 3' primer sequences"
  else
    mkdir -p ${OUT}/${j}
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type/paired/
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type/unpaired_F/
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type/unpaired_R/
    
    echo ${OUT}/QC/cutadapt_fastq/primer_sort/${j}_*_Paired_1.fastq
    for pair in "${filename_pairs[@]}"
    do
      echo "Processing pair: $pair"
      IFS=',' read -r forward_file reverse_file <<< "$pair"
      echo "$forward_file $reverse_file"
      # Remove everything after first space in ID row
      sed -i 's/^\(@[^ ]*\) .*/\1/' $forward_file
      sed -i 's/^\(@[^ ]*\) .*/\1/' $reverse_file
      # For each sample and each metabarcode, this python script checks to see if the forward and reverse files have read pairs, or singleton F or R reads.  Reads are then sorted into the directories generated above.
      python3 ${DB}/scripts/check_paired.py $forward_file $reverse_file ${OUT}/${j}/${j}_sort_by_read_type/paired ${OUT}/${j}/${j}_sort_by_read_type/unpaired_F/ ${OUT}/${j}/${j}_sort_by_read_type/unpaired_R/
      echo ${k} "...check!"
     done
     date
   fi
done

# rm -r ${OUT}/QC/fastq # remove intermediate directories
# rm -r ${OUT}/QC/cutadapt_fastq # remove intermediate directories
# rm -r ${OUT}/QC

###############################
# Make sure unassembled reads are still paired and submit dada2 jobs
###############################

echo " "
echo "Process metabarcode reads for with dada2"
for j in `ls ${OUT}/`
do
 if [[ "${j}" != "QC" && "${j}" != "Run_info" ]]; # ignore non-metabarcode folders...
 then
    #make folders for the metabarcode specific output of dada2 and bowtie2
  echo ""
  echo "${j}"
	mkdir -p ${OUT}/${j}/${j}dada2_out

  echo "Running Dada2"
  printf "#!/bin/bash\n ${RUNNER} ${DB}/scripts/run_dada2.sh -o ${OUT} -d ${DB} -m ${j} -t paired -b ${MINTIMES_ASV:=$MIN_ASV_ABUNDANCE} \n" > ${OUT}/Run_info/run_scripts/${j}_dada2_paired_job.sh
  printf "#!/bin/bash\n ${RUNNER} ${DB}/scripts/run_dada2.sh -o ${OUT} -d ${DB} -m ${j} -t forward -b ${MINTIMES_ASV:=$MIN_ASV_ABUNDANCE} \n" > ${OUT}/Run_info/run_scripts/${j}_dada2_F_job.sh
  printf "#!/bin/bash\n ${RUNNER} ${DB}/scripts/run_dada2.sh -o ${OUT} -d ${DB} -m ${j} -t reverse -b ${MINTIMES_ASV:=$MIN_ASV_ABUNDANCE} \n" > ${OUT}/Run_info/run_scripts/${j}_dada2_R_job.sh
  ${RUNNER} ${OUT}/Run_info/run_scripts/${j}_dada2_paired_job.sh
  date
  ${RUNNER} ${OUT}/Run_info/run_scripts/${j}_dada2_F_job.sh
  date
  ${RUNNER} ${OUT}/Run_info/run_scripts/${j}_dada2_R_job.sh
  date
 fi
done

echo "All done"