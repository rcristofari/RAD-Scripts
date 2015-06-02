#! /cluster/software/R/bin/Rscript
args <- commandArgs(TRUE)

#Exception if no argument is provided
if(length(args) < 1) {
  args <- "--help"
}

###########################################################################################
#	#HELP SECTION
###########################################################################################

if("--help" %in% args) {
  cat("
	SAM_KeepOnlyPairs

	Filters a SAM file in order to retain only properly paired reads. Paired reads may become orphaned after MAPQ filtering and removal of one mate with poor mapping. Removing these improves the duplicate filtering efficiency in Picard MarkDuplicates. KeepOnlyPairs should be applied after MAPQ filtering and sorting (eg in SAMtools sort).
	
	Input requires either a list of SAM files to be processed, one per line (absolute path), or a single SAM file.
 
	Arguments:
	--list=		path to the SAM input files list
	--S=		path to a single SAM input file (once only)
	--out=		path to write the output files
	--trim		trim the trailing part of sequence names (_1 and _2)
	--discards	write discarded reads to a separate SAM file
	--embed		sets the verbosity for embedding in a pipeline
	--no-header	specify if SAM file does not have header (default NULL)
	--help          print this text
 
	Example:
	./SAM_KeepOnlyPairs.R --SAM=~/path/to/SAMlist.txt --out=~/path/to/outdir --no-header --discards\n\n")
 
q(save="no")}

###########################################################################################
#	#STANDARD ARGUMENT PARSER
###########################################################################################

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## --list error
if((is.null(argsL$list)&is.null(argsL$S))) {
  stop("No input file has been provided")}
## --out error
if(is.null(argsL$out)) {
  stop("No output folder has been provided")}
##--no-header default
if(is.null(argsL$no_header)) {
  header_flag <- 0
} else { header_flag <- 1 }

##--orphans default
if(is.null(argsL$discards)) {
  orphans_flag <- 0
} else { orphans_flag <- 1 }
##--embed default
if(is.null(argsL$embed)) {
  embed_flag <- 0
} else { embed_flag <- 1 }
##--trim default
if(is.null(argsL$trim)) {
  trim_flag <- 0
} else { trim_flag <- 1 }

if(is.null(argsL$list)){
argsL$S -> path_list}

if(is.null(argsL$S)){
argsL$list -> path_list}

argsL$out -> path_out

###########################################################################################
#	#DISPLAY RECAP
###########################################################################################

if (embed_flag==0){
if(is.null(argsL$list)){
	cat("\n\t\tSAM input file:\t\t", eval(path_list),"\n", sep="", collapse="")
} else {
	cat("\t\tSAM file list:\t\t", eval(path_list), "\n", sep="", collapse="")}

	cat("\t\tOutput directory:\t", eval(path_out), "\n", sep="", collapse="")

if(is.null(argsL$no_header)){
	cat("\t\tSAM files contain standard header\n")
} else {
	cat("\t\tSAM files do not contain header\n")}
if(is.null(argsL$discards)) {
	cat("\t\tDiscarded reads will be erased\n\n")
} else {
	cat("\t\tDiscarded reads will be written to a separate file\n\n")}
}

###########################################################################################
#	#SETTING UP THE PATHS
###########################################################################################

#	#CASE 1: SAM LIST

if(is.null(argsL$S)){

#Check whether the data file exist
	if(file.exists(eval(path_list))==FALSE){
		stop("SAM input file list not found")}

#Check whether all SAM files exist
	samlist <- as.character((read.table(eval(path_list), header=FALSE))[,1])
	nfiles <- length(samlist)
	
	for(check in 1:nfiles){
		if(file.exists(eval(samlist[check]))==FALSE){
			stop("At least 1 file in SAM list could not be found")}}

cat(paste("Input contains ", nfiles, " SAM files.\n\n", sep="", collapse=""))
}

#	#CASE 2: SAM FILE

if(is.null(argsL$list)){
	if(file.exists(eval(path_list))==FALSE){
		stop("SAM input file list not found")}
samlist <- path_list
nfiles <- 1
}

###########################################################################################
#	#SETTING UP THE WORKING ENVIRONMENT
###########################################################################################

#Setting up the work environment
	system(paste("mkdir -p ", path_out, "/tmp", sep="", collapse="")) #creates a "tmp" folder in the output directory
	suppressPackageStartupMessages(suppressWarnings(library(sqldf)))
	suppressPackageStartupMessages(suppressWarnings(library(tcltk)))

#Then the loop proper
	for(sample in 1:nfiles){
if (embed_flag==0){
cat("Processing sample ", sample, " of ", nfiles, " (", basename(eval(samlist[sample])), ")...\n", sep="", collapse="")
}

#Formatting the SAM file for import:
#The final flags in the SAM file aren't used at all in the RAD pipeline, so we can just get rid if them here and save some bytes (and trouble):
if(header_flag==0){
	system(paste("cat ", eval(samlist[sample]), " | LC_ALL=C grep -v '^@' | cut -f 1-11 > ", eval(path_out), "/tmp/", basename(eval(samlist[sample])), ".df", sep="", collapse="")) } else {
	system(paste("cat ", eval(samlist[sample]), " | cut -f 1-11 > ", eval(path_out), "/tmp/", basename(eval(samlist[sample])), ".df", sep="", collapse=""))}

###########################################################################################
#	#IMPORTING DATA INTO R
###########################################################################################

#Reading the file into R
if (embed_flag==0){cat("\tReading sample...")}
#Prior scan to speed up the import. Setting the expected classes.
	df_path <- paste(eval(path_out), "/tmp/", basename(eval(samlist[sample])), ".df", sep="", collapse="")
	classes <- c('character', 'numeric', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character')
	nrow_sam<-(as.numeric(strsplit(system(paste("wc -l", eval(df_path), collapse=""), intern=T),' ')[[1]][1]))

#Importing the file itself.
#	read.csv(eval(df_path), comment.char='', header=FALSE, sep='\t', colClasses = classes, nrows=nrow_sam)->sam
	data <- file(eval(df_path))
	sam <- sqldf("select * from data", file.format = list(header=FALSE, sep='\t', colClasses = classes))
	colnames(sam)<-c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL')
if (embed_flag==0){cat("done.\n")}

#Fixing the sequence names in the SAM file by removing the trailing part
if (trim_flag==1){
	fix.name <- function(x){strtrim(x, (nchar(x)-2))}
	fix.name(sam$QNAME)->sam$QNAME
}

###########################################################################################
#	#IDENTIFYING PROPER PAIRS
###########################################################################################

#Identifying the second-in-pair reads (indices)
if (embed_flag==0){cat("\tIdentifying proper pairs...")}
	which(duplicated(sam$QNAME)==TRUE)->pairs_2
	sam$QNAME[pairs_2]->pairs_id

#Discarded reads:
	n_discards <- length(sam$QNAME)-2*length(pairs_id)
if (embed_flag==0){cat("done.\n\tDiscarded ", n_discards, " orphan reads.\n", sep="", collapse="")}

###########################################################################################
#	#PULLING OUT PROPER PAIRS
###########################################################################################

#In case "--orphans" option is set to 1:
	if(orphans_flag==1){
if (embed_flag==0){cat("\tWriting discarded reads to file...")}	
	
	#Selecting only reads counted once
	table(sam$QNAME)->read_count
		which(read_count==1)->orphans
		names(orphans)->orphans_id
		sam[(sam$QNAME %in% orphans_id),]->orphaned_reads
		write.table(orphaned_reads, eval(paste(path_out, basename(eval(samlist[sample])), ".discards", sep="", collapse="")), sep='\t', quote=F, row.names=F, col.names=F)
		if (embed_flag==0){cat("...done.\n")}}

	sam[(sam$QNAME %in% pairs_id),]->complete_pairs
if (embed_flag==0){cat("\tWriting to output...")}

###########################################################################################
#	#PRINTING OUT THE CLEANED SAM FILE
###########################################################################################

	write.table(complete_pairs, eval(paste(path_out, "/tmp/", basename(eval(samlist[sample])), ".data", sep="", collapse="")), sep='\t', quote=F, row.names=F, col.names=F)

#Extracting the SAM header into a separate file:
if(header_flag==0){
	system(paste("cat ", eval(samlist[sample]), " | grep '^@' > ", eval(path_out), "/tmp/", basename(eval(samlist[sample])), ".header", sep="", collapse=""))

#Concatenating the data file and the header, cleaning up extra tabs:
	system(paste("cat ", path_out, "/tmp/", basename(eval(samlist[sample])), ".header ", path_out, "/tmp/", basename(eval(samlist[sample])), ".data | sed -E 's/[\t]+/\t/g' | sed -E 's/[\t]$//g'> ",  path_out, basename(eval(samlist[sample])), ".pairs", sep="", collapse=""))

} else {
	system(paste("mv ", path_out, "/tmp/", basename(eval(samlist[sample])), ".data ",  path_out, "/", basename(eval(samlist[sample])), ".pairs", sep="", collapse=""))}

if (embed_flag==0){cat("done.\n\n")}

#Closing the main loop:

}

###########################################################################################
#	#CLEANING UP AND FINISHING
###########################################################################################

	system(paste("rm -r ", eval(path_out), "/tmp", sep="", collapse=""))
	if (embed_flag==0){cat("Finished. Output is in", path_out, "\n\n")}


q(save="no")
