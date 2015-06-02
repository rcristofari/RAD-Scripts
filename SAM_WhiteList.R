#! /usr/bin/Rscript
options(warn=-1) #suppress warnings during execution

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
	SAM_WhiteList

	Filters a SAM file in order to retain only loci included in a whitelist. This list can be extracted, for example, from facetting of a Stacks database. Output BAM are typically suitable for analysis in ANGSD.
	
	Input requires either a list of SAM files to be processed, one per line (absolute path), or a single SAM file. The whitelist of loci must be given in the following format: one row per locus, scaffold ID and mapping position, tab-separated. A tsv file from Stacks export_sql.pl script can also be input directly.
 
	Arguments:
	--list=		path to the SAM input files list
	--S=		path to a single SAM input file (once only)
	--out=		path to write the output files
	--wl=		whitelist of loci to be retained (one column for scaffold and one for start position)
	--export=	export file from Stacks export_sql.pl script
	--help      print this text
 
	Example:
	./SAM_WhiteList.R --S=~/path/to/SAMfile.sam --out=~/path/to/outdir --wl=~/path/to/whitelist.wl\n\n")
 
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
##--wl error  
if((is.null(argsL$wl)&is.null(argsL$export))) {
	stop("No whitelist has been provided")}

###########################################################################################
#	#SETTING UP THE PATHS
###########################################################################################

require(proto, quietly=TRUE)
require(RSQLite, quietly=TRUE)
require(sqldf, quietly=TRUE)
require(tcltk, quietly=TRUE)

#Setting the path_list
if(is.null(argsL$list)){
	argsL$S -> path_list}
if(is.null(argsL$S)){
	argsL$list -> path_list}

#Setting the path_out
argsL$out -> path_out

#Setting the path_whitelist
if(is.null(argsL$export)){
	argsL$wl -> path_whitelist}
if(is.null(argsL$wl)){
	argsL$export -> path_whitelist}

###########################################################################################
#	#DISPLAY RECAP
###########################################################################################	

if(is.null(argsL$list)){
	cat("\n\t\tSAM input file:\t\t", eval(path_list),"\n", sep="", collapse="")
} else {
	cat("\n\t\tSAM file list:\t\t", eval(path_list), "\n", sep="", collapse="")}
	
if(is.null(argsL$export)){
	cat("\t\tLocus whitelist:\t", eval(path_whitelist),"\n", sep="", collapse="")
} else {
	cat("\t\tStacks export file:\t", eval(path_whitelist), "\n", sep="", collapse="")}

	cat("\t\tOutput directory:\t", eval(path_out), "\n\n", sep="", collapse="")
	
###########################################################################################		
#	#CHECKING AND FORMATTING THE INPUT
###########################################################################################

#Check the sample input
if(is.null(argsL$list)){
	if(file.exists(eval(path_list))==FALSE){
		stop("SAM input file not found.")}}
if(is.null(argsL$S)){
	if(file.exists(eval(path_list))==FALSE){
		stop("SAM input file list not found.")}}

#Check the whitelist input
if(is.null(argsL$export)){
	if(file.exists(eval(path_whitelist))==FALSE){
		stop("Locus whitelist not found.")}}
if(is.null(argsL$wl)){
	if(file.exists(eval(path_whitelist))==FALSE){
		stop("Stacks export file not found.")}}
	
#Check the output directory
if(file.exists(eval(path_out))==FALSE){
		stop("Output directory not found")}

#Reading in the sample list
if(is.null(argsL$S)){
	read.csv(eval(path_list), header=F)->sam.list} else {
	path_list->sam.list}

#Reading the whitelist
if(is.null(argsL$export)){
	read.csv(eval(path_whitelist), header=F, sep='\t')->whitelist
} else {
	paste("count=`head ", path_whitelist, " -n 1 | cut -f 13- | tr '\t' '\n' | wc -l` ; c=$(($count+1)) ; cat ", path_whitelist, " | grep -v '#' | cut -f 3 | head -n -$c", sep='', collapse='')->command
	system(eval(command), intern=TRUE)->whitelist.1
	paste("count=`head ", path_whitelist, " -n 1 | cut -f 13- | tr '\t' '\n' | wc -l` ; c=$(($count+1)) ; cat ", path_whitelist, " | grep -v '#' | cut -f 4 | head -n -$c", sep='', collapse='')->command
	system(eval(command), intern=TRUE)->whitelist.2
	data.frame(as.character(whitelist.1), as.numeric(whitelist.2))->whitelist
	names(whitelist)<-c('V1', 'V2')
}
	as.numeric(whitelist$V2)+1->whitelist$V2
	wl.index<-paste(whitelist$V1, '_', whitelist$V2, sep='')	#Composite list including scaffold name and mapping position

#Measuring the lists
wl.length<-nrow(whitelist)
sam.length<-nrow(sam.list)

#Check input files one by one
for(c in 1:sam.length){
	as.character((sam.list[c,1]))->sam.file
	if(file.exists(eval(sam.file))==FALSE){
		stop("At least 1 file in SAM list could not be found.")}}

	if(sam.length==1){
cat("\t\tFound 1 SAM file.\n", sep='')} else {cat("\t\tFound ", sam.length, " SAM files.\n", sep='')}
cat("\t\tWhitelist contains ", wl.length, " loci.\n\n", sep='')

###########################################################################################	
#	#BEGINNING OF THE SAMPLE LOOP
###########################################################################################

for(s in 1:sam.length){

#For each sample, three handles : SAMPLE.full (full path), SAMPLE.path (path only), SAMPLE.base (basename)
	as.character(sam.list[s,1])->SAMPLE.full
	dirname(as.character(sam.list[s,1]))->SAMPLE.path
	strsplit((basename(as.character(sam.list[s,1]))), '\\.')[[1]][1]->SAMPLE.base

cat("Processing sample ", s, " of ", sam.length, " (", SAMPLE.base, ")...", sep='', collapse='')

#Importing the file itself.
cat('\n\tImporting data file...')
	data <- file(SAMPLE.full)
	classes <- c('character', 'numeric', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character')
	sam <- sqldf("select * from data", file.format = list(header=FALSE, sep='\t', comment.char='', colClasses = classes))
	colnames(sam)<-c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'OPT1', 'OPT2')
cat('done.\n')

#Indexing the tables
cat('\tExtracting whitelisted entries...')	
	sam.index<-paste(sam$RNAME, '_', sam$POS, sep='')	#Composite name including the scaffold name and mapping position

#Extracting the corresponding entries
	keep.index <- which(sam.index %in% wl.index)
	sam[keep.index,]->out
cat('done.\n')

#Writing out to SAM file
cat('\tWriting to ', SAMPLE.base,'.wl.sam...', sep='')
	paste(path_out, "/", SAMPLE.base,".wl.sam", sep='', collapse='')->path
	write.table(out, eval(path), col.names=F, row.names=F, quote=F, sep='\t')	
cat('done.\n\n')

}

###########################################################################################
#	#END OF THE SAMPLE LOOP
###########################################################################################

quit('no')