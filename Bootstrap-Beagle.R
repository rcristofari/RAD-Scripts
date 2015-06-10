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

	Bootstrap Beagle genotype likelihoods
	_______________________________________________________________________________________________

	Resample a genotype likelihood file to produce bootstrap replicates in ngsAdmix.
 
	Arguments:
	--likes=		beagle-format genotype likelihoods from ANGSD
	--bootstrap=	number of bootstrap replicates [10]
	--out=		output file path and prefix.
	--offset=		ID number of the first replicate [1]
	--help         	print this text.
 
	Example:
	
	./Bootstrap-Beagle.R --likes=./beagle.lhoods --bootstrap=100 --out=./replicates
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
 
## --likes error
if(is.null(argsL$likes)) {
  stop("You forgot to provide an input file")
} else { argsL$likes -> path_likes }

## --bootstrap default
if(is.null(argsL$bootstrap)) {
	bootstrap <- 10
} else { bootstrap <- as.numeric(argsL$bootstrap) }

## --offset default
if(is.null(argsL$offset)) {
	offset <- 0
} else { offset <- (as.numeric(argsL$offset)-1) }
	
## --out default
if(is.null(argsL$out)) {
  paste(dirname(argsL$likes), '/', basename(argsL$bootstrap), sep='', collapse='')->path_out
} else { argsL$out ->  path_out }

###########################################################################################
#Setup environment

	suppressPackageStartupMessages(suppressWarnings(library(sqldf)))
	suppressPackageStartupMessages(suppressWarnings(library(tcltk)))

###########################################################################################
#Read the beagle file

	cat("\n\n\tReading likelihoods file...")
	data <- file(eval(path_likes))
	likes <- sqldf("select * from data", file.format = list(header=TRUE, sep='\t'))
	cat("done.")
	
###########################################################################################
#Prepare bootstrapped indices	
	nrow(likes)->n.likes
	seq(1,n.likes,1)->likes.id
	
###########################################################################################
#Perform bootstrapping and print out
	
	cat("\n\tBootstrapping...\n")
pb <- txtProgressBar(min = 1, max = bootstrap, style = 3)

	for(B in 1:bootstrap){
		sample(likes.id, n.likes, replace=T)->indices
		likes[indices,]->likes.boots
		write.table(likes.boots, paste(path_out, "_", (B+offset), ".lhoods", sep="", collapse=""), sep='\t', quote=F, col.names=T, row.names=F)
		setTxtProgressBar(pb,B)
	}


cat('\n\n\tΤΑ ΔΗ ΝΥΝ ΠΑΝΤΑ ΤΕΛΕΙΤΑΙ\n\n')
q(save="no")

