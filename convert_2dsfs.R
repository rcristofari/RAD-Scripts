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
	_______________________________________________________________________________________________

	Convert 2d-SFS do natural exponential (for piping).
	
	Converts a log-scaled 2d-sfs generated in ANGSD to natural exponential. Use --norm to scale
	it up to something else than 1.
 
	Arguments:
	--sfs=		a log-scaled joint-sfs produced by realSFS.
	--norm=		normalise the sfs to a sum [1]
	--out=		output file path and name.
	--help         	print this text.
 
	Example:
	
	./convert-2dsfs.R --sfs=~/POP1_POP2.2dsfs --norm=1000000 --out=./POP1_POP2.exp.2dsfs
	_______________________________________________________________________________________________	

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
if(is.null(argsL$sfs)) {
  stop("You forgot to provide an input file")}

## --norm default
if(is.null(argsL$norm)) {
	norm.out <- 1
} else { norm.out <- as.numeric(argsL$norm) }
	
## --out default
if(is.null(argsL$out)) {
  paste(dirname(argsL$sfs), '/', basename(argsL$sfs), '.exp', sep='', collapse='')->path_out
} else { argsL$out ->  path_out }

argsL$sfs -> path_sfs

###########################################################################################
#Check whether the list file exist
if(file.exists(eval(path_sfs))==FALSE){
	stop("2d-SFS file not found.")}
	
############################################################################################################################

read.table(eval(path_sfs), sep=' ', header=F)->sfs
Filter(function(x)!all(is.na(x)), sfs)->sfs
exp(sfs)*norm.out->sfs

write.table(sfs, eval(path_out), quote=F, row.names=F, col.names=F)
q(save="no")
		
