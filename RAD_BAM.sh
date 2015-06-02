#! /bin/bash

USAGE="Usage: `basename $0` [-hv] [-d deduplicate] -i input_dir -o output_dir -q min_mapq (0-42) -y output_type (stacks or paired) -p /path/to/picard -s /path/to/samtools -r /path/to/scripts"

####################################################
## Standard getopt argument parser and help line ##
####################################################

# Parse command line options.

while getopts hvi:o:q:y:p:s:r:d OPT; do
    case "$OPT" in
        h)
            echo -ne "\n\n\tRAD-sequencing BAM pipeline\n\n\tThis wrapper pipes RAD sequencing reads that have been mapped to a reference genome into BAM files that are directly usable for downstream analysis. This may be catalog building in Stacks, SNP calling in GATK, or SNP-free analysis in ANGSD. The process goes as follows: 1) paired reads are filtered for mapping quality, 2) orphaned reads are discarded and only proper pairs retained, 3) PCR and optical duplicates are removed and 4) selected reads are output.\n\n\tAs an input, we need SAM files that have been typically produced through alignement of paired-end Illumina fastq files onto a reference genome (this may be a RADome produced from a draft catalog assembly), using Bowtie2, with the following options: [--fr -I 250 -X 1000 --no-discordant --no-mixed --phred33 --no-unal]. This ensures properly paired, mapped and cleaned SAM files.\n\n\tIt is assumed that you have SAMtools and PicardTools installed, as well as the SAM_KeepOnlyPairs.R script. You need to specify install paths using the -p, -s and -r options.\n\n\t" $USAGE "\n\n"
            exit 0
            ;;
        v)
            echo "`basename $0` version 0.1"
            exit 0
            ;;
        d)
            DEDUP=1
            ;;
        i)
            INPUT_DIR=$OPTARG
            ;;
        o)
            OUTPUT_DIR=$OPTARG
            ;;
        q)
            MAPQ=$OPTARG
            ;;
        y)
            TYPE=$OPTARG
            ;;
        p)
            PICARD=$OPTARG
            ;;
        s)
            SAMTOOLS=$OPTARG
            ;;
        r)
            SCRIPTS=$OPTARG
            ;;
        \?)
            # getopts issues an error message
            echo $USAGE >&2
            exit 1
            ;;
    esac
done

# Remove the switches we parsed above.
shift `expr $OPTIND - 1`

# We want at least one non-option argument. 
# Remove this block if you don't need it.

#if [ $# -eq 0 ]; then
#    echo $USAGE >&2
#    exit 1
#fi

# Access additional arguments as usual through 
# variables $@, $*, $1, $2, etc. or using this loop:

for PARAM; do
    echo $PARAM
done

#####################################################################################
#Clean up the path names to remove the final slash
INPUT_DIR=`echo $INPUT_DIR | sed -E 's/[/]$//g'`
OUTPUT_DIR=`echo $OUTPUT_DIR | sed -E 's/[/]$//g'`
PICARD=`echo $PICARD | sed -E 's/[/]$//g'`
SAMTOOLS=`echo $SAMTOOLS | sed -E 's/[/]$//g'`
SCRIPTS=`echo $SCRIPTS | sed -E 's/[/]$//g'`
#####################################################################################
#Print a summary of the options
#####################################################################################


#Listing the SAM files in the input directory
SAM_LIST=`ls -1 $INPUT_DIR | grep '.sam' | sed -E 's/([[:alnum:]]+)(.sam)/\1/g'`

#Setting up the work environment
mkdir -p $OUTPUT_DIR/BAM
mkdir -p $OUTPUT_DIR/log

#####################################################################################
#ADD A FILE CHECKUP STEP TO AVOID LOOSING TIME
#ADD VERBOSITY
#PRINT RESULTS OF CHECK. NB OF FILES FOUND
#PIPE LOG FILES AND PRINT STATISTICS OF READ LOSS ALONG THE PIPELINE
#####################################################################################

#Beginning a loop through these files. SAM files will be processed one by one into finished BAM.
for SAMPLE in $SAM_LIST ; do
echo -e "\n\n\tProcessing sample" $SAMPLE"..."


#####################################################################################
#Filtering with minimum MAPQ
	#Takes in "SAMPLE.sam", takes out "SAMPLE.mapq.sam"
echo -ne "\t\tFiltering out MAPQ under" $MAPQ"..."
$SAMTOOLS/samtools view -h -q $MAPQ -S $INPUT_DIR/$SAMPLE.sam -o $OUTPUT_DIR/BAM/$SAMPLE.mapq.sam
echo -e "done."


#####################################################################################
#Filtering out orphaned reads
	#Takes in "SAMPLE.mapq.sam", takes out "SAMPLE.pairs.sam"
echo -ne "\t\tFiltering out orphaned reads..."

$SCRIPTS/SAM_KeepOnlyPairs.R --S=$OUTPUT_DIR/BAM/$SAMPLE.mapq.sam --out=$OUTPUT_DIR/BAM/ --trim --embed 2>&1 | tee $OUTPUT_DIR/log/$SAMPLE.orphans.log

mv $OUTPUT_DIR/BAM/$SAMPLE.mapq.sam.pairs $OUTPUT_DIR/BAM/$SAMPLE.pairs.sam
	rm $OUTPUT_DIR/BAM/$SAMPLE.mapq.sam
echo -e "done."


#####################################################################################
#Converting SAM to BAM
	#Takes in "SAMPLE.pairs.sam", takes out "SAMPLE.pairs.bam"
echo -ne "\t\tConverting to BAM..."
$SAMTOOLS/samtools view -b -S $OUTPUT_DIR/BAM/$SAMPLE.pairs.sam -o $OUTPUT_DIR/BAM/$SAMPLE.pairs.bam
	rm $OUTPUT_DIR/BAM/$SAMPLE.pairs.sam
echo -e "done."


#####################################################################################
#Adding read groups
#Need to find a way to add this to PATH for portability
	#Takes in "SAMPLES.pairs.bam", takes out "SAMPLE.group.bam"
echo -ne "\t\tAdding read groups..."
java -jar $PICARD/AddOrReplaceReadGroups.jar  I= $OUTPUT_DIR/BAM/$SAMPLE.pairs.bam O= $OUTPUT_DIR/BAM/$SAMPLE.group.bam LB= RAD-SAMPLE PL= ILLUMINA PU= RADSEQ SM= $SAMPLE QUIET=TRUE VERBOSITY=ERROR
	rm $OUTPUT_DIR/BAM/$SAMPLE.pairs.bam
echo -e "done."


#####################################################################################
#Sorting reads
	#Takes in "SAMPLE.group.bam", takes out "SAMPLE.sort.bam"
echo -ne "\t\tSorting reads by coordinates..."
$SAMTOOLS/samtools sort $OUTPUT_DIR/BAM/$SAMPLE.group.bam $OUTPUT_DIR/BAM/$SAMPLE.sort
	rm $OUTPUT_DIR/BAM/$SAMPLE.group.bam
echo -e "done."


#####################################################################################
#Removing duplicates
if [ $DEDUP -eq 1 ] ; then
#Need to find a way to add this to PATH for portability
	#Takes in "SAMPLE.sort.bam", takes out "SAMPLE.dedup.bam"
echo -ne "\t\tRemoving duplicates..."
java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTPUT_DIR/BAM/$SAMPLE.sort.bam OUTPUT=$OUTPUT_DIR/BAM/$SAMPLE.dedup.bam METRICS_FILE=$OUTPUT_DIR/log/$SAMPLE.dedup.metrics REMOVE_DUPLICATES=true READ_NAME_REGEX='[0-9]_([0-9]+)_([0-9]+)_([0-9]+)_paired' QUIET=TRUE VERBOSITY=ERROR
	rm $OUTPUT_DIR/BAM/$SAMPLE.sort.bam
fi
echo -e "done."


#####################################################################################
#Export only the selected reads
	#Takes in "SAMPLE.dedup.bam", takes out "SAMPLE.TYPE.bam"
#Two cases: if -d is on, we import $SAMPLE.dedup.bam ; id -d is off, we import $SAMPLE.sort.bam

echo -ne "\t\tExporting for" $TYPE"..."
if [ "$TYPE" == stacks ] ; then EXPORT_FILTER=64 ; fi
if [ "$TYPE" == paired ] ; then EXPORT_FILTER=1 ; fi
if [ $DEDUP -eq 1 ]
then $SAMTOOLS/samtools view -f $EXPORT_FILTER -b $OUTPUT_DIR/BAM/$SAMPLE.dedup.bam -o $OUTPUT_DIR/BAM/$SAMPLE.$TYPE.bam
else $SAMTOOLS/samtools view -f $EXPORT_FILTER -b $OUTPUT_DIR/BAM/$SAMPLE.sort.bam -o $OUTPUT_DIR/BAM/$SAMPLE.$TYPE.bam
fi
echo -e "done."


#####################################################################################
#Index the final BAM file
	#Makes a "SAMPLE.TYPE.bam.bai" file.
echo -ne "\t\tIndexing..."
$SAMTOOLS/samtools index $OUTPUT_DIR/BAM/$SAMPLE.$TYPE.bam
echo -e "done."


#####################################################################################
#End of the main processing loop
done

echo -ne "\t\tFinished."

