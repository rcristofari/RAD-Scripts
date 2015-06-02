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
	MIGRATE-n INPUT FILE BUILDER
	
	Formats a series of block-files ouput by the RAD_Haplotypes.R script with 'migrate' option into a single MIGRATE-n input file.
	The script assumes that all loci have an equal number of populations and haplotypes (true if the blocks come from
	RAD_Haplotypes.R). A list can easily be obtained with \" ls -1 | sort -R | head -n 50 \" for example.
 
	Arguments:
	--list=		a list of block files to be combined, one per line.
	--out=		output file path and name.
	--title=	a brief description of the data.
	--help         	print this text.
 
	Example:
	
	./parse.fasta_engine-1.0.R --list=~/random.list --out=./migrate.infile --title='Sequence data - 50 random loci'
	
	
")	
  q(save="no")

}

###########################################################################################
#Standard argument parser

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## --list error
if(is.null(argsL$list)) {
  stop("You forgot to provide an input file")}

## --out default
if(is.null(argsL$out)) {
  paste(dirname(argsL$list), '/migrate.infile', sep='', collapse='')->path_out
} else { argsL$out ->  path_out }

## --title default
if(is.null(argsL$title)) {
  file.string <- ''
} else { argsL$title -> file.string }

argsL$list -> path_filelist
file.string<- paste(file.string, ' (generated with MIGRATE.infile.R)', sep='', collapse='')


###########################################################################################
#Check whether the list file exist
if(file.exists(eval(path_filelist))==FALSE){
	stop("File list not found.")}
	
############################################################################################################################

read.table(eval(path_filelist), header=F)->filelist
as.character(filelist$V1)->filelist
length(filelist)->nfiles

###########################################################################################
#Check whether individual input files exist

for(i in 1:nfiles){

if(file.exists(eval(filelist[i]))==FALSE){
	stop(paste("Input file ", filelist[i], " could not be found.", sep='', collapse=''))}
}

############################################################################################################################
#Extract basic stats from the first file in the list

read.table(eval(filelist[1]), sep='\n',header=F)->trailer
trailer[seq(1,nrow(trailer), 2),]->header
trailer[seq(2,nrow(trailer), 2),]->sequences
gsub('>','',header)->header
strsplit(header,'_')->header.list
data.frame(matrix(unlist(header.list), ncol=2, byrow=T))->header.cols
data.frame(header.cols$X1, header.cols$X2, sequences)->locus
names(locus)<-c('pop','hap','sequence')

length(unique(locus$pop))->npop
nchar(as.character(locus$sequence[1]))->seq.length
paste(npop,nfiles,file.string)->header.line1
paste(as.character(rep(seq.length, nfiles)), collapse=" ")->header.line2
newline<-'\n'

############################################################################################################################
#Initiate the output file
rbind(header.line1, header.line2)->file.header
write.table(file.header, eval(path_out), row.names=F, col.names=F, quote=F)

############################################################################################################################
#Read in all migrate block files
#(in the form of a loop through the provided list) and format them into lists

data.list <- list()

for(i in 1:nfiles){
read.table(eval(filelist[i]), sep='\n',header=F)->fasta
fasta[seq(1,nrow(fasta), 2),]->header
fasta[seq(2,nrow(fasta), 2),]->sequences
gsub('>','',header)->header
strsplit(header,'_')->header.list
data.frame(matrix(unlist(header.list), ncol=2, byrow=T))->header.cols
data.frame(header.cols$X1, header.cols$X2, sequences)->locus
names(locus)<-c('pop','hap','sequence')
split(locus, locus$pop)->locus.list
data.list[[i]]<-locus.list
}


############################################################################################################################
#Interleave and append to the output file

cat('Writing to', path_out)

for(p in 1:npop){
	pop.df<-data.frame()
	#set a header line for the population block and print it
	nrow(data.list[[1]][[p]][1])->nind
	as.character(data.list[[1]][[p]][[1]][1])->pop.name
	paste(newline, nind, ' ', pop.name, sep='', collapse='')->pop.header
	write.table(pop.header, eval(path_out), row.names=F, col.names=F, quote=F, append=T)
	
	for(l in 1:nfiles){
		as.data.frame(data.list[[l]][p])->df
		names(df)<-c('pop','hap','sequence')
		rbind(pop.df,df)->pop.df
		}
	write.table(pop.df, eval(path_out), row.names=F, col.names=F, quote=F, append=T)
}

cat('\nSuccess..!\n\n')
q(save="no")
		
