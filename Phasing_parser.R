#! /usr/bin/Rscript
options(warn=-1) #suppress warnings during execution
args <- commandArgs(TRUE)
#Exception if no argument is provided
if(length(args) < 1) {
  args <- "--help"
}

	suppressPackageStartupMessages(suppressWarnings(library(sqldf)))
	suppressPackageStartupMessages(suppressWarnings(library(tcltk)))

###########################################################################################
#	#HELP SECTION
###########################################################################################

if("--help" %in% args) {
  cat("
	_____________________________________________________________________________________________________

	HapCompass phasing parser
	_____________________________________________________________________________________________________

	Takes the output of HapCompass physical phasing (short haplotype building) and integrates it as phased
	haplotypes into the full VCF file, for downstram analysis.
	Will only work on single-sample, single-chromosome files at a time.
	_____________________________________________________________________________________________________

	Arguments:
	--vcf=				VCF file containing a single sample.
	--phase=			HapCompass phasing output for that sample.
	--out=				Path and filename to write the output.
	
	--help      		print this text
 
	Example:
	./Phasing_parser.R --vcf=scaffold_unphased.vcf --phase=scaffold_MWER_sxolution.txt --out=scaffold
	
 	_____________________________________________________________________________________________________\n\n")
q(save="no")}

###########################################################################################
#	#STANDARD ARGUMENT PARSER
###########################################################################################

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

#required arguments:
#path_vcf <- 'King-1_Scaffold15.vcf'
#path_phase <- 'King-1_Scaffold15_MWER_solution.txt'
#path_out
#chr <- "Scaffold15"

#--path_vcf error:
if(is.null(argsL$vcf)){
	stop("You didn't provide any input VCF")
} else {
	argsL$vcf->path_vcf
}

#--path_phase error:
if(is.null(argsL$vcf)){
	stop("You didn't provide any input VCF")
} else {
	argsL$phase->path_phase
}

#--path_out default:
if(is.null(argsL$out)){
	path_out <- paste(dirname(path_vcf), "/", strsplit(path_vcf, "\\.")[[1]][1], sep='')
} else {
	argsL$out->path_out
}

#--debug default:
if(is.null(argsL$debug)){
	debug <- 0
} else {
	debug <- 1
}

cat("\n\tVCF file: ", path_vcf, sep='')
cat("\n\tPhase file: ", path_phase, sep='')
cat("\n\tOutput file: ", path_out, sep='')


#########################################################################################################################
#VCF FILE AND PAHSE FILE IMPORT AND FORMAT
#########################################################################################################################

##Read the VCF file:
cat("\n\n\tReading VCF file...")
	paste("cat ", path_vcf, " | grep -v '##' > ", path_vcf, ".tmp", sep='', collapse='')->cmd
	system(cmd)
	paste(path_vcf, '.tmp', sep='', collapse='') -> path_tmp
	read.table(eval(path_tmp), sep='\t', header=T, comment.char='')->vcf # en sqldf
	paste("rm ", path_vcf, ".tmp", sep='', collapse='')->cmd
	system(cmd)
	as.character(vcf[,9])->vcf[,9]
	'GT'->vcf[,9]
	as.character(vcf[,10])->vcf[,10]
	unlist(strsplit(as.character(vcf[,10]), ':'))->split
	split[seq(1,length(split),2)]->vcf[,10]
	vcf[vcf[,10]=="1/1",10]<-"1|1"
	vcf[vcf[,10]=="0/0",10]<-"0|0"
	n.sites <- nrow(vcf)
cat("done.\n")

##Read the phasing file and parse it into blocks:
#Split the phase file into phase blocks.
cat("\tParsing the phase file...")
	paste("mkdir .temp -p ; rm -f .temp/* ; cat ", path_phase, " | tr '\n' '&' | sed -E 's/&&/\\n/g' | tr '&' '\n' | awk -v RS=\"BLOCK\" 'NR > 1 { print RS $0 > \"temp\" (NR-1); close(\"temp\" (NR-1)) }' ; mv temp* .temp/", sep='')->cmd
	system(cmd)
#List the phase-block files
	list.files('.temp')->block.files
	length(block.files)->n.blocks
cat("done.\n")

#########################################################################################################################
#LOOP THROUGH THE BLOCKS AND STREAM OUT THE VCF FILE
#########################################################################################################################

#Initialise the VCF file by copying in the header information:
	paste("cat ", path_vcf, " | grep '#' > ", path_out, ".phased.vcf", sep="")->cmd
	system(cmd)

#Get the start position of the first block:
	read.table(".temp/temp1", comment.char="B", header=F)->first.block
	which(vcf$POS==first.block$V2[1])-1 ->first.block.start

#Write the beginning of the (unphased) data:
	if(first.block.start>0){
	write.table(vcf[1:first.block.start,], paste(path_out, ".phased.vcf", sep=""), quote=F, row.names=F, col.names=F, sep='\t', append=T)
	}

cat("\tWriting phased VCF...\n\n")

#Set the progress bar
	pb <- txtProgressBar(min = 0, max = n.blocks, style = 3)		#set the progress bar

#Start looping through the block files
for(b in 1:n.blocks){

#Read and format the phasing block
	read.table(paste(".temp/temp", b, sep=''), comment.char="B", header=F)->block
	data.frame(POS=block$V2, Hap1=block$V4, Hap2=block$V5)->block
	block$POS[1]->block.start
	block$POS[length(block$POS)]->block.end

#Also read the next block to get its start position:
	if(b<n.blocks){
		read.table(paste(".temp/temp", b+1, sep=''), comment.char="B", header=F)->next.block
		which(vcf$POS==next.block$V2[1])-1->next.block.start
	} else { next.block.start <- n.sites }

#Turn it so the first row is always 0/1, and the rest phased accordingly.
#Add a condition for the last block.
	if(block$Hap1[1]==1){ 
		paste(block$Hap2, block$Hap1, sep='/')->PHASE
	} else if(block$Hap1[1]==0){
		paste(block$Hap1, block$Hap2, sep='/')->PHASE}
		gsub("/","|", PHASE[2:length(PHASE)])->PHASE[2:length(PHASE)]


#Extract the corresponding (longer) block from the VCF file:
	which(vcf$POS==block.start)->vcf.block.start
	which(vcf$POS==block.end)->vcf.block.end
	vcf[vcf.block.start:vcf.block.end,]->vcf.block
	gsub('/','|',vcf.block[,10])->vcf.block[,10]
	vcf.block[vcf.block$POS %in% block$POS,10]<-as.character(PHASE)

vcf.block->vcf[vcf.block.start:vcf.block.end,]

#To account for cases of overlap, switch the next call as unphased.
gsub('\\|', '/', vcf[(vcf.block.end+1), 10])-> vcf[(vcf.block.end+1), 10]


#In order to write by blocks. More efficient memory-wise, but problem with overlapping / included blocks.
#if(debug==1){
#		beacon<-paste("Phased start: ",vcf.block.start,", phased stop: ",vcf.block.end, sep='')
#		write.table(beacon, paste(path_out, ".phased.vcf", sep=""), quote=F, row.names=F, col.names=F, sep='\t', append=T)
#}
#
##Write the phased data:
#	write.table(vcf.block, paste(path_out, ".phased.vcf", sep=""), quote=F, row.names=F, col.names=F, sep='\t', append=T)
#if(debug==1){
#		beacon<-paste("Unphased start: ",vcf.block.end+1,", unphased stop: ",next.block.start, sep='')
#		write.table(beacon, paste(path_out, ".phased.vcf", sep=""), quote=F, row.names=F, col.names=F, sep='\t', append=T)
#}
#
##Write the buffering unphased data:
#	if(next.block.start >= (vcf.block.end+1)){
#	write.table(vcf[(vcf.block.end+1):next.block.start,], paste(path_out, ".phased.vcf", sep=""), quote=F, row.names=F, col.names=F, sep='\t', append=T)
#	if(debug==1){
#		beacon<-paste("End of block")
#		write.table(beacon, paste(path_out, ".phased.vcf", sep=""), quote=F, row.names=F, col.names=F, sep='\t', append=T)
#	}
#	}
	
setTxtProgressBar(pb,b)

}

write.table(vcf, paste(path_out, ".phased.vcf", sep=""), quote=F, row.names=F, col.names=F, sep='\t', append=T)
system("rm -r .temp/")

cat('\n\tΤΑ ΔΗ ΝΥΝ ΠΑΝΤΑ ΤΕΛΕΙΤΑΙ\n\n')

quit('no')
