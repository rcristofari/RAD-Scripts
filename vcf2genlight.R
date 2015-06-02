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
	vcf2genlight

	Converts a genotype call file produced by ANGSD into a VCF format, used for analysis in other packages, or for SNP call comparison using VCFtools perl scripts. Note that it will set a quality score of 30 by default.
	Genotypes should be in 012 format, as output by ANGSD using the -doGeno 3 option (Major, minor and calls).
 
	Arguments:
	--vcf=		path to the VCF file to convert
	--out=		output file name and path
	--popmap=	population map (optional)
	--help      	print this text
 
	Example:
	./vcf2genlight.R --vcf=~/path/to/data.vcf --out=~/path/to/outfile.snp\n\n")
 
q(save="no")}

###########################################################################################
#	#STANDARD ARGUMENT PARSER
###########################################################################################

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
 ##--vcf error
if(is.null(argsL$vcf)) {
	stop("No VCF file has been provided")}
## --out error
if(is.null(argsL$out)) {
  stop("No output file has been provided")}
##--popmap option
if(is.null(argsL$popmap)==F){
	argsL$popmap -> path_popmap
	pop_flag<-1
} else{ pop_flag <- 0 }

###########################################################################################
#	#SETTING UP THE PATHS
###########################################################################################
argsL$vcf -> path_vcf
argsL$out -> path_out
	
###########################################################################################		
#	#CHECKING INPUT
###########################################################################################

#Check the input
if(file.exists(eval(path_vcf))==FALSE){stop("VCF file list not found.")}
#Check the output directory
if(file.exists(dirname(eval(path_out)))==FALSE){stop("Output directory not found")}

###########################################################################################		
#	#CONVERTING FORMAT
###########################################################################################

cat("\n\tReading VCF file...")
paste("cat ", path_vcf, " | grep -v '##' | cut -f 10- > ", path_vcf, ".tmp", sep='', collapse='')->cmd
system(cmd)
paste(path_vcf, '.tmp', sep='', collapse='') -> path_tmp
read.table(eval(path_tmp), sep='\t', header=T, comment.char='')->vcf
paste("rm ", path_vcf, ".tmp", sep='', collapse='')->cmd
system(cmd)
cat("done.\n")

if(pop_flag==1){
cat("\tReading the population map...")
read.table(path_popmap, sep='\t', header=F)->popmap
cat("done.\n")
}

cat("\tParsing data...")
#A complicated way to keep only the genotypes for each individual
as.list(vcf)->vcf.list
vcf.split <- function(x){strsplit(as.character(x), ':')[[1]][1]}
vcf.apply <- function(x){unlist(lapply(x, vcf.split))}
lapply(vcf.list, vcf.apply)->gen.list
cat("done.\n")

#Then convert them to 012 format
sub.<-function(x){gsub('./.', '-', x)}
sub0<-function(x){gsub('0/0','0',x)}
sub1<-function(x){gsub('0/1','1',x)}
sub2<-function(x){gsub('1/1','2',x)}

cat("\tFormatting SNP file...\n")
cat('\t|##        |\t20%')
lapply(gen.list, sub0)->gen.list
cat('\r\t|####      |\t40%')
lapply(gen.list, sub1)->gen.list
cat('\r\t|######    |\t60%')
lapply(gen.list, sub2)->gen.list
cat('\r\t|########  |\t80%')
lapply(gen.list, sub.)->gen.list
cat('\r\t|##########|\t100%\t\n\tLoading data...')
lapply(gen.list, as.numeric)->gen.list
#And transform that into a genotype matrix
do.call(cbind.data.frame, gen.list)->gen.df
t(as.matrix(gen.df))->geno
rep(">", length(rownames(geno)))->prefix
rep("\n", length(rownames(geno)))->return
paste(prefix, rownames(geno), return, sep='')->rownames(geno)
cat('done.\n\tWriting genlight file...')

line1 <- ">>>> begin comments - do not remove this line <<<<"
line2<-paste("Generated with vcf2genlight.R from ", path_vcf, ".", sep='', collapse='')
line3<-">>>> end comments - do not remove this line <<<<"
line6<-">> ploidy\n2"

if(pop_flag==0){
rbind(line1,line2,line3, line6)->header

} else if (pop_flag==1){
line4<-">> population"
line5<-paste(as.character(popmap[,2]), sep=' ', collapse=' ')
rbind(line1,line2,line3,line4,line5, line6)->header
}
write.table(header, path_out, quote=F, sep='', row.names=F, col.names=F, na='-', append=F)
write.table(geno, path_out, quote=F, sep='', col.names=F, na='-', append=T)
cat("done.\n")
cat('\n\tΤΑ ΔΗ ΝΥΝ ΠΑΝΤΑ ΤΕΛΕΙΤΑΙ\n\n')
