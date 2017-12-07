#! /usr/bin/Rscript
args <- commandArgs(TRUE)
options (warn=-1)

#Exception if no argument is provided
if(length(args) < 1) {
  args <- "--help"
}

############################################################################################################################
## HELP SECTION
############################################################################################################################

if("--help" %in% args) {
  cat("
  	_____________________________________________________________________________________________________

	RAD Haplotype tool
	_____________________________________________________________________________________________________
	
	Parses a Stacks population fasta file, in order to 1) export locus files for downstream analysis, and
	2) perform locus-based summary statistic caclulation (number of haplotypes, Tajima's D, nucleotide
	diversity, pairwise mismatch distribution, segregating site density distribution).
	
	Output can be fasta, nexus, or alternatively fasta-like blocks suitable for building a MIGRATE-n input
	file using the MIGRATE_infile.R script. Output can also be set to 'silent'. The tool is designed to be
	used in a population-coalescent framework: sequences are drawn randomly from the population haplotype
	pool. If coverage is not rock-solid, it is possible to downsample to haploid individuals (--ploidy=1).
 	
 	Population map is a tab-separated text file: first field contains unique individual names, second field
 	3-letter population codes, tab-separated. Use ignore_pop to ignore that information and output a single
 	file per locus. 'Keep' list is a simple column of individual names to retain if analysis is restricted
 	to a subset of samples.
 	
 	You can choose the number of samples to output using nind. In case of multiple populations, specify the
 	number of samples per population (coma-separated). Use '-1' to keep all available haplotypes. Specified
 	numbers can be exact, or a minimum number (--min flag). If only one number is specified, it is applied
 	to all populations. If ignore_pop is used, only one value is needed. 
 	_____________________________________________________________________________________________________
	
	Arguments:
	
	--fasta=	the input whole-database fasta file from Stacks
	--map=		a *complete* population map (all samples must be present)
	--keep=		an optional subset of individuals to restrict the analysis
	--nind=		number of invididuals to keep per locus and per population (coma-separated) [-1]
	--min		is nind a minimum number or individuals (as opposed to exact)?
	--ploidy=	number of chromosomes to randomly sample per individual [1]	
	--stats		will calculate and output summary statistics on a per-locus basis
	--out=		path to output file [same as fasta]
	--save=		prefix to save the data output as R lists to be loaded in an R session [null]
	--type=		'fasta', 'nexus', 'migrate' or 'silent' [nexus, or silent in stats mode]
	--ignore_pop=	logical. If 1, population information will be discarded (fasta or nexus only) [0]
	
	--help         	print this text
 
	Example:
	
	./RAD_Haplotypes.R --fasta=batch_1.fa --map=~/paired/emperor/emp.popmap --out=./ --type='migrate' 
	--nind=8,4,-1,8 --min --ploidy=1 --stats
 	_____________________________________________________________________________________________________
	
	
")	
  q(save="no")

}

############################################################################################################################
## STANDARD ARGUMENT PARSER
############################################################################################################################

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## --fasta error
if(is.null(argsL$fasta)) {
  stop("You forgot to provide an input file")}
  
## --map error
if(is.null(argsL$map)) {
  stop("You forgot to provide a population map")}
  
## --nind error
if(is.null(argsL$nind)) {
	stop("You must provide a number of individuals to output")
} else { as.numeric(unlist(strsplit(argsL$nind,','))) -> nind }  

## --min default
if(is.null(argsL$min)) {
	min_nind <- 0
} else { 1 -> min_nind }  

## --out default
if(is.null(argsL$out)) {
  cbind(dirname(argsL$fasta), '/RADhap', sep='')->path_out
} else { argsL$out ->  path_out }

## --save default
if(is.null(argsL$save)) {
  do_save <- 0
} else { 
	do_save <- 1
	argsL$save ->  path_save }

## --type default
if(is.null(argsL$type) & is.null(argsL$stats)) {
  'nexus'->out_type
} else if(is.null(argsL$type) & is.null(argsL$stats)==F) {
	'silent'->out_type
} else { argsL$type ->  out_type }

## --ploidy default
if(is.null(argsL$ploidy)) {
  ploidy <- 1
} else { as.numeric(argsL$ploidy) -> ploidy }

## --keep default
if(is.null(argsL$keep)){
	keep<-0
} else { keep <- 1 
	argsL$keep->path_keep
}

## --ignore_pop default
if(is.null(argsL$ignore_pop)) {
  ignore <- 0
} else { as.numeric(argsL$ignore_pop) -> ignore }

## --ignore_pop error
if(ignore==1 & out_type=='migrate') {
	stop("You can't ignore population information for migrate blocks")}

## --stats default
if(is.null(argsL$stats)){
	stats<-0
} else { stats <- 1 }

argsL$fasta -> path_in
argsL$map -> path_map

############################################################################################################################
## DISPLAY TASK SUMMARY
############################################################################################################################

cat('_____________________________________________________________________________________________________\n')
cat('\nInput fasta: ', basename(path_in), sep='')
cat('\nInput populations: ', basename(path_map), sep='')
if(keep==1){
cat('\nKeep samples in: ', basename(path_keep), sep='')}
if(ploidy==1){
cat('\nDownsampling to haploid individuals')}
if(ignore==1){
cat('\nIgnoring population information')}
if(stats==1){
cat('\nComputing per-locus summary statistics')}
if(out_type != 'silent'){
cat('\nWriting to ', out_type, ' format', sep='')
} else if (out_type=='silent'){
cat('\nSuppressing locus alignment file output')}
cat('\nWriting out files to ', path_out, sep='')
if(do_save==1){
cat('\nSaving to Rdata file')}
cat('\n_____________________________________________________________________________________________________\n')

############################################################################################################################
## SET UP THE ENVIRONMENT
############################################################################################################################

##########
#Check whether the data files exist
	if(file.exists(eval(path_in))==FALSE){
		stop("Fasta input file doesn't exist")}
	if(file.exists(eval(path_map))==FALSE){
		stop("Population map file doesn't exist")}
	if(file.exists(eval(path_out))==FALSE){
		stop("Output directory doesn't exist")}

##########
#Load the required packages
	suppressMessages(require(adegenet, quietly=T))
	suppressMessages(require(ape, quietly=T))
	suppressMessages(require(pegas, quietly=T))
	suppressMessages(require(seqinr, quietly=T))
	suppressPackageStartupMessages(suppressWarnings(library(sqldf)))
	suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
	suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
	suppressPackageStartupMessages(suppressWarnings(library(gridExtra)))

##########
#Correct the ape function

write.nexus.correct <- function(x, file, format = "dna", datablock = TRUE, interleaved = TRUE, charsperline = NULL, gap = NULL, missing = NULL)  {
    format <- match.arg(toupper(format), c("DNA", "PROTEIN"))
    indent <- "  "
    maxtax <- 5
    defcharsperline <- 80
    defgap <- "-"
    defmissing <- "?"
    if (is.matrix(x)) {
        if (inherits(x, "DNAbin")) 
            x <- as.list(x)
        else {
            xbak <- x
            x <- vector("list", nrow(xbak))
            for (i in seq_along(x)) x[[i]] <- xbak[i, ]
            names(x) <- rownames(xbak)
            rm(xbak)
        }
    }
    ntax <- length(x)
    nchars <- length(x[[1]])
    zz <- file(file, "w")
    if (is.null(names(x))) 
        names(x) <- as.character(1:ntax)
    fcat <- function(..., file = zz) cat(..., file = file, sep = "", 
        append = TRUE)
    find.max.length <- function(x) max(nchar(x))
    print.matrix <- function(x, dindent = "    ") {
        Names <- names(x)
        printlength <- find.max.length(Names) + 2
        if (!interleaved) {
            for (i in seq_along(x)) {
                sequence <- paste(x[[i]], collapse = "")
                taxon <- Names[i]
                thestring <- sprintf("%-*s%s%s", printlength, 
                  taxon, dindent, sequence)
                fcat(indent, indent, thestring, "\n")
            }
        }
        else {
            ntimes <- ceiling(nchars/charsperline)
            start <- 1
            end <- charsperline
            for (j in seq_len(ntimes)) {
                for (i in seq_along(x)) {
                  sequence <- paste(x[[i]][start:end], collapse = "")
                  taxon <- Names[i]
                  thestring <- sprintf("%-*s%s%s", printlength, 
                    taxon, dindent, sequence)
                  fcat(indent, indent, thestring, "\n")
                }
                if (j < ntimes) 
                  fcat("\n")
                start <- start + charsperline
                end <- end + charsperline
                if (end > nchars) 
                  end <- nchars
            }
        }
    }
    fcat("#NEXUS\n[Data written by write.nexus.data.R, ", date(), 
        "]\n")
    NCHAR <- paste("NCHAR=", nchars, sep = "")
    NTAX <- paste0("NTAX=", ntax)
    DATATYPE <- paste0("DATATYPE=", format)
    if (is.null(charsperline)) {
        if (nchars <= defcharsperline) {
            charsperline <- nchars
            interleaved <- FALSE
        }
        else charsperline <- defcharsperline
    }
    if (is.null(missing)) 
        missing <- defmissing
    MISSING <- paste0("MISSING=", missing)
    if (is.null(gap)) 
        gap <- defgap
    GAP <- paste0("GAP=", gap)
    INTERLEAVE <- if (interleaved) 
        "INTERLEAVE=YES"
    else "INTERLEAVE=NO"
    if (datablock) {
        fcat("BEGIN DATA;\n")
        fcat(indent, "DIMENSIONS ", NTAX, " ", NCHAR, ";\n")
        fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", 
            GAP, " ", INTERLEAVE, ";\n")
        fcat(indent, "MATRIX\n")
        print.matrix(x)
        fcat(indent, ";\nEND;\n\n")
    }
    else {
        fcat("BEGIN TAXA;\n")
        fcat(indent, "DIMENSIONS", " ", NTAX, ";\n")
        fcat(indent, "TAXLABELS\n")
        fcat(indent, indent)
        j <- 0
        for (i in seq_len(ntax)) {
            fcat(names(x[i]), " ")
            j <- j + 1
            if (j == maxtax) {
                fcat("\n", indent, indent)
                j <- 0
            }
        }
        fcat("\n", indent, ";\n")
        fcat("END;\n\nBEGIN CHARACTERS;\n")
        fcat(indent, "DIMENSIONS", " ", NCHAR, ";\n")
        fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", DATATYPE, 
            " ", INTERLEAVE, ";\n")
        fcat(indent, "MATRIX\n")
        print.matrix(x)
        fcat(indent, ";\nEND;\n\n")
    }
    close(zz)
}


############################################################################################################################
## TASK 0: READ THE DATA FILES (FASTA, MAP, KEEP LIST)
############################################################################################################################

##########
#Read Fasta file
	cat("\nLoading fasta file...")
	data <- file(eval(path_in))
	sqldf("select * from data", file.format = list(header=FALSE, sep='\n'))->fasta
	closeAllConnections()

##########
#Parse the data
	cat("done.\nFormatting fasta file...")
	fasta[seq(1,nrow(fasta), 2),]->header
	gsub('>CLocus_|\\[|\\]\\,', '', header)->header
	fasta[seq(2,nrow(fasta), 2),]->sequences
	strsplit(as.character(header), '_| ')->splits
	lapply(splits, "[", c(1:7))->splits
	data.frame(matrix(unlist(splits), nrow=length(header), byrow=T))->metadata
##########
#Combine it into a dataframe
	data.frame(metadata$X1,metadata$X3,metadata$X7, sequences)->data
	names(data)<-c("Locus","Sample","Allele","Sequence")
	as.numeric(levels(data$Locus))[data$Locus]->data$Locus
	as.numeric(levels(data$Sample))[data$Sample]->data$Sample
	as.numeric(levels(data$Allele))[data$Allele]->data$Allele

##########
#Fixing the real sample and population names from the population map
	cat("done.\nParsing populations...")
	read.csv(eval(path_map), sep='\t', header=F) -> map
		#In case the ignore_pop flag is on, population info is discarded 
		#and we pick the sum of all asked sequences at random.
			if(ignore==1){map[,2]<-'ALL'
			nind <- nind[1]}
	as.character(map[,1])->Name
	as.character(map[,2])->Pop
	seq(1,length(Name),1)->Sample
	data.frame(Name, Pop, Sample) -> map

##########
#Additional filtering in case we wish to restrict the analysis to some individuals only
	if(keep==1){
		read.table(path_keep)->keep.list
		as.data.frame(keep.list)->keep.list
		nrow(keep.list)->nkeep
		cat("\rParsing populations: keeping ", nkeep, " samples...", sep='')
		names(keep.list)<-'Name'
		getwd()
		merge(map, keep.list)->map
	}

##########
#Merging the data and the population map
	merge(data, map, by="Sample")->data
	length(unique(map[,2])) -> map.npop
	
##########
#Checking if the correct number of nind has been provided and correct nind if needed
if(ignore==0){
	if(length(nind)<map.npop){
	rep(nind[1], map.npop)->nind}	
}

##########
#Splitting the data into la list, per locus
	cat("done.\nSorting out individual loci...")
	split(data, data$Locus)->fasta.list
	length(names(fasta.list))->nloci
	cat("done.\n")
#Clear up some memory before proceeding
rm(fasta)
rm(data)
rm(sequences)
rm(metadata)
rm(header)
rm(splits)
#gc()

############################################################################################################################
## END OF TASK 0
############################################################################################################################


############################################################################################################################
## TASK 1: BEGINNING OF THE CORE LOOP THROUGH LOCI
############################################################################################################################

cat('\n')

##########
#Initialize the list to store alignments
	align.list<-list()
	loc.names<-vector()
	
##########
#Initialize counting variables
	notenough <- 0
	popdrop <- 0
	polypl <- 0
	snps <- numeric(96)
	align.index <- 0

##########
#Begin the locus loop
	for(l in 1:nloci){
	cat('\rParsing locus ', l ,' of ', nloci , ', discarded ', (notenough+popdrop+polypl) , '.', sep='', collapse='')
	locus <- fasta.list[[l]]
	align.index <- align.index+1

##########
#Split the locus into a list, by population
		split(locus, as.character(locus$Pop))->pop.list
		length(names(pop.list))->npop
		numeric(npop)->nhap.pop
		for(p in 1:npop){length(unique(pop.list[[p]]$Name))->nhap.pop[p]}

##########
#Control whether each population has enough data
	if(npop<map.npop){			
	#this means that this locus is not present in every mapped population.
		popdrop <- popdrop+1
		align.index-1 -> align.index
	} else if(sum(nind>nhap.pop)>0) {
	#this means at least one pop does not have enough sequences.
		notenough <- notenough+1
		align.index-1 -> align.index
	} else {

##########
#Remove pseudopolyploid loci: Stacks numbers alleles [0,1,..,n] but we want only [0,1].
	if(length(which(locus$Allele>1))<1){

	##########
	#Initialize the locus-list
	align.list[[align.index]]<-list()

	##########
	#Now we can parse a name for that locus
	paste("Locus_",locus$Locus[1], "_", locus$Chromosome[1] ,":", locus$Position[1], sep='')->loc.name
	c(loc.names, loc.name)->loc.names
	
	########################################################################################################################
	## EMBEDDED POPULATION LOOP - CORE TASK
	########################################################################################################################

	##########
	#Prepare a vector to receive concatenated population data
	align.list[[align.index]][[npop+1]] <- as.alignment()
	align.list[[align.index]][[npop+1]]$nb<-0

	##########
	#Start the loop
	for(p in 1:npop){
	locus <- pop.list[[p]]


	##########
	#Re-establish the missing alleles for homozygotes
		table(locus$Sample)->tab
		tab[tab==1]->tab #these are the homozygotes (they have just one recorded sequence)
		locus[which(locus$Sample %in% as.numeric(names(tab))),]->haploid
		rbind(locus, haploid)->diploid
		diploid[order(diploid$Sample),]->fixed

	##########
	#Sample one allele at random if required
	if(ploidy==1){
		nrow(fixed)/2->n.haplo
		seq(1,2*n.haplo,2)->index
		sample(c(1,2), 1, replace=F)->random
		for(i in 1:nind[p]){sample(c(0,1), 1)+index[i]->index[i]}
		fixed[index,]->haploid_subset
		##########
		#Sample the specified number of sequences from the haploid pool
			if(nind[p]>=0){
				if(min_nind==0){
			sample(seq(1,n.haplo,1), nind[p], replace=F)->nind_subset
			haploid_subset[sort(nind_subset),]->subset
				} else if (min_nind==1){
				haploid_subset -> subset			
				}
			} else if (nind[p] < 0) { haploid_subset -> subset }

	##########
	#Alternatively, keep the whole diploid sample	
	} else if(ploidy==2){
		seq(1,nrow(fixed),1)->index
		##########
		#Sample the specified number of sequences from the diploid pool
		if(nind[p]>=0){
			if(min_nind==0){
		sample(index, nind, replace=F)->random
		fixed[random,] -> subset
			} else if (min_nind==1){
			fixed -> subset
			}
		} else if (nind[p] < 0) { fixed -> subset }
	}
		
	##########
	#Format that locus into a printable DNA block		
		parse.name <- function(x){paste('>Locus_',as.numeric(x[2]),'_',x[8], sep='', collapse='')}
		apply(subset, 1, parse.name)->fasta.names
		data.frame(names=fasta.names, sequences=subset$Sequence, com=subset$Pop)->output

	##########
	#Make it a DNA object
		length(output$names)->nseq
		as.alignment(nseq, as.character(output$names), as.character(output$sequences), as.character(output$com))->align

	##########
	#Fit it into the general result list
		align.list[[align.index]][[p]]<-align
		align.list[[align.index]][[npop+1]]->all.align
		as.alignment(as.numeric(all.align$nb)+align$nb, c(all.align$nam, align$nam), c(all.align$seq, align$seq), c(all.align$com, align$com))->align.list[[align.index]][[npop+1]]

	} # end of the population loop.
	names(align.list[[align.index]])<-c(names(pop.list),'ALL')
	########################################################################################################################
	## END OF EMBEDDED POPULATION LOOP
	########################################################################################################################

##########
#End of the missing-pop and pseudopolyploid exceptions
	} else { 
	polypl <- polypl+1
	align.index-1 -> align.index
	}
} # end of the 'missing locus' conditions
##########

} # end of the locus loop.

names(align.list)<-loc.names
c(as.character(unique(sort(map$Pop))),'ALL')->pop.names

cat('\rParsed ', nloci, ' loci, discarded ', polypl, ' as polyploid, ', popdrop, ' as missing in some populations, and ', notenough, ' as not having enough sequences. Kept ', length(align.list), ' loci.\n', sep='')

############################################################################################################################
## END OF THE CORE LOOP THROUGH LOCI - align.list now contains all loci split by pop
############################################################################################################################

#The connexion goes through the align.list object

############################################################################################################################
## TASK 2: EXTRACTING SUMMARY STATISTICS
############################################################################################################################
if(stats==1){
cat('\nCalculating summary statistics.\n')
#We have align.list, a list of 'nloci' loci, each including 'npop' populations + all combined, in seqinr alignment format.

#Summary statistics will be formatted in a different way: a list with one dataframe per pop+all, with always the same
#number of nows (=nloci), and NA if test can't be performed / data is missing.
#Locus names are row names, and we have: V1 number of variable loci, V2 the number of haplotypes, V3-V5 Tajima test, V6 nucleotide diversity

#Mismatch distribution will be calculated in two steps: (1) a list of locus-pop mismatches, and (2) quantiles from that
#distribution, split by number of variable sites in the locus.

##########
#These need to be reset
length(pop.names)->npop
length(align.list)->nloci

##########
#Initialize the sumstats list
sumstats <- list()
for(p in 1:npop){sumstats[[p]]<-vector()}
names(sumstats)<-pop.names

##########
#Initialize the mismatch distribution list
mismatch <- list()
for(p in 1:npop){mismatch[[p]]<-vector()}
names(mismatch)<-pop.names

##########
#Loop through the locus list to retrieve alignments
#set the progress bar
pb <- txtProgressBar(min = 1, max = nloci, style = 3)

for(l in 1:nloci){
	for(p in 1:npop){
		align.list[[l]][[p]]->locus.align
		as.DNAbin(locus.align)->locus.bin
	
		#Find the number of sequences
		locus.align$nb -> n_seq
	
		#Find the number of mismatches
		if(sum(dist.alignment(locus.align))==0){n_mismatches<-0
		} else { alignment2genind(locus.align, polyThres=0)->dna
		length(dna@all.names)->n_mismatches}
	
		#Find the number of haplotypes
		length(unique(locus.align$seq))->n_haplo
		
		#Perform Tajima's test
		if(n_mismatches > 0 & locus.align$nb >=4){
		as.numeric(tajima.test(locus.bin))->taj_test
		} else { c(NA,NA,NA)->taj_test}
		
		#Compute nucleootide diversity
		nuc.div(locus.bin)->div
		
		#Make a result vector
		c(n_seq, n_mismatches, n_haplo, taj_test, div)->l.p.stat
	
		rbind(sumstats[[p]], l.p.stat)->sumstats[[p]]
	
		#Calculate the mismatch distribution as a vector of 95+1 cells
		table(round(dist.dna(locus.bin)*95))->pairwise
		pairwise*(1/sum(pairwise))->pairwise
		numeric(97)->m_dist
		pairwise->m_dist[1:length(pairwise)]
		m_dist[97]<-n_mismatches
		rbind(mismatch[[p]], m_dist)->mismatch[[p]]
	
	} #end of the pop loop
	setTxtProgressBar(pb,l)
} #end of the locus loop

for(p in 1:npop){as.data.frame(sumstats[[p]], row.names=loc.names)->sumstats[[p]] ; names(sumstats[[p]])<-c('N_SEQ','N_SEG','N_HAPLO','TAJ_D','TAJ_Pn','TAJ_Pb','PI')}
for(p in 1:npop){as.data.frame(mismatch[[p]], row.names=loc.names)->mismatch[[p]]}

##########
#Writing out to sumstat files
for(p in 1:npop){
	as.data.frame(sumstats[[p]])->df.out
	outname<-paste(path_out, '_', pop.names[p], '.sumstats.tsv', sep='')
	write.table(df.out, outname, quote=F, sep='\t', row.names=T, col.names=T)
	}


##########
#Split the mismatch list according to the number of variant sites
split_mismatch <- function(x){split(x[,1:96], x$V97)}
lapply(mismatch, split_mismatch)->mismatch.split


############################################################
#Figure out mismatch stats per pop
############################################################
cat('\nCalculating pairwise mismatch distributions.\n')
pb <- txtProgressBar(min = 1, max = npop, style = 3)
##########
#A couple functions to sum up the distributions
quantile_c <- function(x){quantile(x, prob=c(0.025, 0.05, 0.5, 0.95, 0.975))}
parse_mismatch <- function(x) {sapply(x, mean)->MEAN ; sapply(x, quantile_c)->quant ; rbind(N_MIS=seq(0,95,1), MEAN, quant)}
fix_mismatches <- function(x){as.data.frame(t(as.matrix(x)))}

distributions<-list()
for(p in 1:npop){
	mismatch.split[[p]]->pop_mm.list
	lapply(pop_mm.list, parse_mismatch)->pop.mm.stats
	lapply(pop.mm.stats, fix_mismatches)->pop.mm.stats
		for(i in 1:length(pop.mm.stats)){names(pop.mm.stats[[i]])<-c('N_MIS','MEAN','L95','L90','MEDIAN','U90','U95')}
	distributions[[p]]<-pop.mm.stats
	setTxtProgressBar(pb,p)
}

names(distributions)<-pop.names

cat('\nWriting out to PDF.')
##########
#Loop through everything to create a list of plots
mismatch.plots<-list()

for(p in 1:npop){
	mismatch.plots[[p]]<-list()
	length(mismatch.split[[p]])->nclass
	for(c in 2:nclass){
		distributions[[p]][[c]]->data
		names(distributions[[p]])[c]->S
		names(distributions)[p]->P
		nrow(mismatch.split[[p]][[c]])->N
		
##########
#ggPlot definition
plot.title <- paste('Population ', P, ', mismatch distribution for category ', S, '\n(based on ', N,' loci)\n', sep='')
ggplot(data, aes(x=N_MIS)) +
	geom_line(aes(y=MEDIAN), col='brown') +
	geom_ribbon(aes(ymin=L95, ymax=U95), fill='tan', alpha=0.25) +
	geom_line(aes(y=L95), col='brown', size=0.1, linetype='dashed') +
	geom_line(aes(y=U95), col='brown', size=0.1, linetype='dashed') +
	coord_cartesian(xlim=c(0,8), ylim=c(0,1)) +
	scale_y_continuous(breaks=c(seq(from = 0, to = 1, by = 0.1))) +
	scale_x_continuous(breaks=c(seq(from = 0, to = 8, by = 1))) +
	labs(list(title=plot.title, x='Number of mismatches', y='Frequency'), vjust=1) +
	theme(plot.title=element_text(size=16,lineheight=.8,vjust=1.2), axis.text.x=element_text(size=12,angle=90), axis.text.y=element_text(size=12), axis.title.y=element_text(size=12,lineheight=.8,vjust=0.8)) +
	theme_bw()->g
	
	g->mismatch.plots[[p]][[c-1]]
}
}

#Plot as PDF files
for(p in 1:npop){
length(mismatch.plots[[p]])->n.plots
mismatch.plots[[p]]->plot.list
suppressMessages(pdf(paste(path_out, '_', pop.names[p], '.pdf', sep=''), width=21, height=27))
	i<-1
	plot<-list() 
	for (n in 1:n.plots){
  	plot[[i]]<-plot.list[[n]]
  		if(i %% 6 == 0) {
    		do.call(grid.arrange, plot)
    		plot<-list()
    		i<-0}
  	i<-i+1}  	
if (length(plot) != 0) {do.call(grid.arrange, plot)}
suppressMessages(dev.off())
}


############################################################
#Figure out segregating sites
############################################################
cat('\nCalculating segregating site distribution.\n')
pb <- txtProgressBar(min = 1, max = npop, style = 3)
seg.plots<-list()
for(p in 1:npop){
	sumstats[[p]]$N_SEG->segsites
	table(segsites)->obs
	numeric(96)->observed
	table(segsites)->observed[1:length(obs)]
	mean(segsites)->lambda.obs
	observed*(1/sum(observed))->observed

	matrix(nrow=10000,ncol=96, data=0)->res
	for(i in 1:10000){as.numeric(table(rpois(length(segsites),lambda.obs)))->tmp ; tmp*(1/sum(tmp))->tmp ; tmp -> res[i,1:length(tmp)]}
	matrix(nrow=3, ncol=96, data=0)->quant
	for(i in 1:96){mean(res[,i])->quant[2,i]}
	for(i in 1:96){quantile(res[,i], 0.975)->quant[1,i]}
	for(i in 1:96){quantile(res[,i], 0.025)->quant[3,i]}
	as.data.frame(t(quant))->res
	cbind(seq(0,95,1),res,observed)->res
	names(res)<-c('seq','lower95','mean','upper95','observed')

title <- paste(pop.names[p], " - SNP density distribution - observed (red) vs. Poisson (black)", sep='')
ggplot(res, aes(x=seq))+geom_line(aes(y=mean), size=0.4) +
	geom_ribbon(aes(ymin=lower95,ymax=upper95), alpha=0.2) +
	geom_line(aes(y=observed), col='red', size=0.4) +
	coord_cartesian(xlim=c(0,10)) +
	scale_y_continuous(breaks=c(seq(from = 0, to = 1, by = 0.1))) +
	scale_x_continuous(breaks=c(seq(from = 0, to = 10, by = 1))) +
	xlab("Number of SNPs per locus")+ylab("Number of loci") +
	ggtitle(title) +
	theme(plot.title=element_text(size=12,lineheight=.8,vjust=1),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=12,lineheight=.8,vjust=-0.4),axis.title.y=element_text(size=12,lineheight=.8,vjust=0.8)) -> g
	g -> seg.plots[[p]]
	
	setTxtProgressBar(pb,p)

} #end of poisson loop

cat('\nWriting out to PDF.')
#Plot as PDF files
length(seg.plots)->n.plots
seg.plots->plot.list
pdf(paste(path_out, '_segsites.dist.pdf', sep=''), width=21, height=27)
	i<-1
	plot<-list() 
	for (n in 1:n.plots){
  	plot[[i]]<-plot.list[[n]]
  		if(i %% 6 == 0) {
    		do.call(grid.arrange, plot)
    		plot<-list()
    		i<-0}
  	i<-i+1}  	
if (length(plot) != 0) {do.call(grid.arrange, plot)}
#suppressMessages(dev.off())

} #end of Stats task
############################################################################################################################
## END OF TASK 2
############################################################################################################################

#We still go back to the align.list data structure

############################################################################################################################
## TASK 3: WRITING TO SELECTED OUTPUT TYPE
############################################################################################################################
#We will again need to loop through loci. If output is set to Migrate, we will take the $ALL pop

if(out_type != 'silent'){

cat('\nWriting out loci.\n')
pb <- txtProgressBar(min = 1, max = nloci, style = 3)


##########
#These need to be reset
length(pop.names)->npop
length(align.list)->nloci

matrix(unlist(strsplit(names(align.list), '_')), ncol=3, byrow=T)[,2]->loc.names

##########
#Case 1, fasta or nexus (locus-pop)

if(out_type != 'migrate'){

	for(l in 1:nloci){
		for(p in 1:(npop-1)){
	
			align.list[[l]][[p]]->align
	
			##########
			#Get the number of mismatches in that locus
			if(sum(dist.alignment(align))==0){
				n_mismatches<-0
			} else {
				alignment2genind(align, polyThres=0)->dna
				length(dna@all.names)->n_mismatches
			}
			
		##########
		#Reformat the alignment into a simple dataframe
		data.frame(names=align$nam, sequences=align$seq)->output

		##########
		#Write a Nexus file
		if(out_type=='nexus'){
		
			as.character(output$sequences)->nex
			names(nex)<-output$names
			as.list(nex)->nex
			sapply(nex, function(x){as.vector(strsplit(x,''))})->nex
					if(npop>1){	write.nexus.correct(nex, paste(dirname(path_out), '/locus_', loc.names[l],'_',pop.names[p],'.',n_mismatches,'snps.subset.nex', sep='', collapse=''), format = "dna", interleaved=F)
					} else { write.nexus.correct(nex, paste(dirname(path_out), '/locus_', loc.names[l],'.',n_mismatches,'snps.subset.nex', sep='', collapse=''), format = "dna", interleaved=F)}

		##########
		#Write a Fasta file				
		} else if(out_type=='fasta') {	#Writing out as Fasta file
		
					if(npop>1){ write.table(output, paste(dirname(path_out), '/locus_', loc.names[l],'_',pop.names[p],'.',n_mismatches,'snps.subset.fa', sep='', collapse=''), sep='\n', quote=F, row.names=F, col.names=F)
					} else {	write.table(output, paste(dirname(path_out), '/locus_', loc.names[l],'.',n_mismatches,'snps.subset.fa', sep='', collapse=''), sep='\n', quote=F, row.names=F, col.names=F)}
		}
	} #end of population loop
		setTxtProgressBar(pb,l)
	} #end of the locus loop


##########
#Case 2, migrate (all pops concatenated)

} else if (out_type == 'migrate'){

for(l in 1:nloci){
		align.list[[l]][[npop]]->align

		##########
		#Get the number of mismatches in that locus
		if(sum(dist.alignment(align))==0){
			n_mismatches<-0
		} else {
			alignment2genind(align, polyThres=0)->dna
			length(dna@all.names)->n_mismatches
		}

		align$nb->nseq
		migrate.names<-vector()
		sprintf("_HAP.%02d", seq(1,nseq,1))->hapname
		for(s in 1:nseq){
		paste(">", align$com[s], hapname[s], sep='', collapse='') -> migrate.name
		c(migrate.names, migrate.name)->migrate.names
		}

		cbind(migrate.names, align$seq)->migrate.out
		write.table(migrate.out, paste(dirname(path_out), '/locus_', loc.names[l],'.', n_mismatches,'snps.subset.migrate.fa', sep='', collapse=''), sep='\n', quote=F, row.names=F, col.names=F)

	setTxtProgressBar(pb,l)
} #end of the locus loop
}


} #end of the Silent condition
############################################################################################################################
## END OF TASK 3
############################################################################################################################

if(do_save==1){
cat('\nSaving to Rdata object...')
	save(list=c('align.list', 'sumstats', 'distributions','mismatch.split'), file=paste(path_save, '.RADhap.Rdata', sep=''))
cat('done.\n')
}


cat('\n\tΤΑ ΔΗ ΝΥΝ ΠΑΝΤΑ ΤΕΛΕΙΤΑΙ\n\n')
q(save="no")
