#! /usr/bin/Rscript
args <- commandArgs(TRUE)

#Exception if no argument is provided
if(length(args) < 1) {
  args <- "--help"
}

###########################################################################################

## Help section
if("--help" %in% args) {
  cat("
	RADOME PARSER
	Builds a synthetic genome composed of scaffolds beginning with a single-end RAD-tag
	and ending with the corresponding paired-end contig, padded out with Ns up to the
	average fragment length.
	Input files are:
		- single-end: output from Stacks export_sql.pl script, tsv format, with -F pe=1
		- paired-end: output from Stacks exec_velvet.pl script
	(Requires package seqinr)
 
	Arguments:
	--single=	path to the single-end export file
	--paired=	path to the paired-end contigs file
	--out=		path to output file
	--rev_comp	if specified, the paired-end contigs will be reverse-complemented
	--minDP=	minimum coverage to retain the paired-end contig [default=1]
	--max_length=	maximum allowed paired-end contig length [default: no limit]
	--length=	scaffold final length [default=5 N insert]
	--force_match	if specified, name of paired-ends will not be checked (this may lead to false assembly)
	--help          print this text
 
	Example:
	./RADgenome_parse.R --single=~/path/to/export.tsv --paired=~/path/to/collated.fa --out=~/path/to/genome.fa--minDP=1 --length=750\n\n")
 
  q(save="no")
}

###########################################################################################

#Standard argument parser
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## --single error
if(is.null(argsL$single)) {
  stop("No single-end input file has been provided")}
## --paired error
if(is.null(argsL$paired)) {
  stop("No paired-end input file has been provided")}
## --out error
if(is.null(argsL$out)) {
  stop("No output file has been provided")}
##--rev_comp default
if(is.null(argsL$rev_comp)) {
  RC_flag <- 0
} else { RC_flag <- 1 }
## --minDP default
if(is.null(argsL$minDP)) {
  minDP <- 1
} else { as.numeric(argsL$minDP) -> minDP }

## --length default
if(is.null(argsL$length)) {
  fragment_length <- 1500
  length_switch <- 0
} else { as.numeric(argsL$length) -> fragment_length
	 length_switch <- 1 }

## --max_length default
if(is.null(argsL$max_length)) {
	max_length <- 0
} else { as.numeric(argsL$max_length) -> max_length }

##--force_match default
if(is.null(argsL$force_match)) {
  force_flag <- 0
} else { force_flag <- 1 }

argsL$single -> path_single.raw
argsL$paired -> path_paired.raw
argsL$out -> path_out

###########################################################################################

#Check whether the data files exist
if(file.exists(eval(path_single.raw))==FALSE){
	stop("Single-end input file doesn't exist")}

if(file.exists(eval(path_paired.raw))==FALSE){
	stop("Paired-end input file doesn't exist")}

###########################################################################################

#Display recap, conditionnal on whether arguments have been provided or take their default values.

cat(paste("\n\tSingle-end data file:\t\t", eval(path_single.raw),"\n\tPaired-end data file:\t\t", eval(path_paired.raw),"\n\tOutput file:\t\t\t", eval(path_out), "\n", sep="", collapse=""))

if(is.null(argsL$rev_comp)){
	cat("\tPaired-end contigs:\t\talready in-line\n")
} else {
	cat("\tPaired-end contigs:\t\treverse-complement\n")}

if(is.null(argsL$minDP)) {
	cat("\tMinimum paired-end coverage:\t1x (default)\n")
} else {
	cat(paste("\tMin. paired-end coverage:\t", eval(minDP), "x\n", sep="", collapse=""))}

if(is.null(argsL$max_length)) {
  cat("\tMaximum contig length:\t\tno limit\n")
} else {
	cat(paste("\tMax. paired-end contig length:\t", eval(max_length), "bp\n", sep="", collapse=""))}

if(is.null(argsL$length)) {
	cat("\tOutput scaffold length:\t\tvariable\n")
} else {
	cat(paste("\tOutput scaffold length:\t\t", eval(fragment_length), "bp\n", sep="", collapse=""))}

if(is.null(argsL$force_match)){
cat("\tCheck name matches:\t\tyes\n")
} else {
	cat("\tCheck name matches:\t\twill not check names\n")}


###########################################################################################
#Function "rev.comp" applied on a DNA sequence as a string. Requires the "comp" function from seqinr package.
#Used only if --rev_comp flag is used
library(seqinr)
rev.comp <- function(x){paste(rev(comp(unlist(strsplit(x,"")), forceToLower=F)), collapse="")}
###########################################################################################

#Setting up the work environment
dirname(path_out)->workdir
system(paste("mkdir -p ", workdir, "/tmp", sep="", collapse=""))	#creates a "tmp" folder in the output directory

#Formatting the collated.fa file into a standard no-newline fasta format
system(paste("cat ", eval(path_paired.raw), " | tr -d '\\n' | sed -E 's/>/\\n>/g' | sed -E 's/([0-9]+)([ACTG]{20})/\\1\\t\\2/g' | tr -d '>' | tr '|' '\\t' | tr '_' '\\t' | cut -f 1,3,5,7,8 >", eval(workdir), "/tmp/collated.df", sep="", collapse=""))
#Formatting the export_sql.pl file. Don't forget the -F pe=1 flag to export only reads with a pe contig, to speed things up (shouldn't crash otherwise, however)
system(paste("cat ", eval(path_single.raw), " | cut -f 1,5 > ", eval(workdir), "/tmp/single.df", sep="", collapse=""))

paste(eval(workdir), "/tmp/single.df", sep="", collapse="")->path_single
paste(eval(workdir), "/tmp/collated.df", sep="", collapse="")->path_paired

###########################################################################################
#PARAMETERS FOR HARDCODING / EXECUTING FROM R CONSOLE
#Get the path of work directory and both input files (if harcoded and run from R)
#path_single <- "~/paired/king/tmp/single.df"
#path_paired <- "~/paired/king/tmp/collated.df"
#Set parameters (if harcoded and run from R)
#Min sequencing coverage to retain a paired-end contig
#minDP <- 5
#Output fragment length
#fragment_length <- 1024
#Extracts the folder above the working folder (assumed to be a temp folder)
#output_dir <- paste(strsplit(dirname(path_single),"/")[[1]][-(length(strsplit(dirname(path_single),"/")[[1]]))], collapse="/")
#path_output <- paste(output_dir, "/rad_genome.fa", sep="", collpase="")
###########################################################################################

cat("\n\tFormatting input files...")
#Import the single-end reads in R from export_sql.pl.
#Prior scan to speed up the import
tab2rows <- read.table(eval(path_single), header = FALSE, sep='\t', nrows = 2)
classes <- c('numeric','character')
nrow_single<-(as.numeric(strsplit(system(paste("wc -l", eval(path_single), collapse=""), intern=T),' ')[[1]][1]))-1

#Import and format
read.table(eval(path_single), header=FALSE, sep='\t', colClasses = classes, nrows=nrow_single)->single
colnames(single)<-c('id','sequence')
toupper(single$sequence)->single$sequence	#corrects the occasional case mistakes

#Import the paired-end contigs in R
#Prior scan to speed up the import
tab2rows <- read.table(eval(path_paired), header = FALSE, sep='\t', nrows = 2)
classes <- sapply(tab2rows, class)
classes[5]<-'character'
nrow_paired<-(as.numeric(strsplit(system(paste("wc -l", eval(path_paired), collapse=""), intern=T),' ')[[1]][1]))-1

#Import and format
read.table(eval(path_paired), header=FALSE, sep='\t',colClasses = classes, nrows=nrow_paired)->paired
colnames(paired)<-c('id','node','contig_length','coverage','sequence')
nchar(paired$sequence)->paired$contig_length	#corrects the discrepancy between Stacks output contig length and actual contig length

###########################################################################################

if(RC_flag==1){
#Reverse-complement the sequences
cat("done.\n\n\tReverse-complementing the paired-end contigs...")
unname(sapply(paired[,5], rev.comp))->RC_paired
RC_paired->paired$sequence
}

toupper(paired$sequence)->paired$sequence	#corrects the occasional case mistakes

cat("done.\n")
###########################################################################################

#Create a vector with the index of duplicated id's in the 'paired' table.
which(duplicated(paired$id)==TRUE)->dup_in_paired	#contains the IDs for the duplicated rows in the unsorted 'paired' table
unique(paired[dup_in_paired,1])->dup_id			#contains the matching locus IDs
sort(paired$id[!(paired$id %in% dup_id)])->pass_id	#extracts IDs for the loci which have only 1 paired-end contig
cat("done.\n\n\tRemoving loci with paired-end coverage under threshold...")

###########################################################################################

#Add the optional filter by paired-contig minimum coverage
paired$id[which(paired$coverage < minDP)]->low_cov_id
pass_id[!(pass_id %in% low_cov_id)]->pass_id

cat("done.\n\n")

###########################################################################################

#Add the optional filter by paired-end-contig maximum length
if(is.null(argsL$max_length)) {
  max(nchar(paired$sequence))->max_length }

paired$id[which(paired$contig_length > max_length)]->long_contigs_id
pass_id[!(pass_id %in% long_contigs_id)]->pass_id

###########################################################################################

cat("\tRetained ", length(pass_id), " loci after filtering.\n\n")

###########################################################################################

#To avoid a crash during the parsing, check whether fragment_length will always be long enough

max_scaffold_length <- max(nchar(single$sequence))+max(nchar(paired$sequence))

if(max_scaffold_length >= fragment_length){
	stop("Maximum scaffold length longer than fragment length. Set a higher output fragment length or filter out long contigs.")}


###########################################################################################
cat("\tBuilding the scaffolds...\n\n")
#Keep only the paired-end contigs that passed all test:
paired[paired$id %in% pass_id,] -> paired_keep

#Merge the two dataframes:
merge(single, paired_keep, by='id')->merged


#Core loop to actually pull out the sequences matching these IDs and parse them together with Ns to make pseudo-contigs.
system(paste("rm -f", eval(path_out), collapse=""))			#remove any eventual previous output file


nrow(merged)->n.loc
pb <- txtProgressBar(min = 0, max = n.loc, style = 3)		#set the progress bar

for(locus in 1:n.loc){
	read1 <- merged$sequence.x[locus]
	read2 <- merged$sequence.y[locus]
	contig_length <- merged$contig_length[locus]
	coverage <- merged$coverage[locus]

	seqid <- paste('>Scaffold_', as.character(merged$id[locus]), '\n', sep='')
	
	if(length_switch==1){
	N <- (fragment_length - nchar(read1) - nchar(read2))
	} else { N <- 5 }

	Ns <- paste(sample('N', N, replace=T), collapse='')
	fragment <- paste(seqid,read1, Ns, read2, sep='', collapse='')
	write.table(fragment, eval(path_out), append=T, quote=F, row.names=F, col.names=F)

	setTxtProgressBar(pb,locus)
}

#cat(paste("\n\n\tOutput in ", eval(path_out), "\n\n", sep="", collapse=""))

###########################################################################################
system(paste("rm -r ", eval(workdir), "/tmp", sep="", collapse=""))

cat('\n\n\tΤΑ ΔΗ ΝΥΝ ΠΑΝΤΑ ΤΕΛΕΙΤΑΙ\n\n')
q(save="no")

