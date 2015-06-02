#! /usr/bin/Rscript
args <- commandArgs(TRUE)
options (warn=-1)

#Exception if no argument is provided
if(length(args) < 1) {
  args <- "--help"
}

###########################################################################################
## Help section

if("--help" %in% args) {
  cat("
  	_____________________________________________________________________________________________________

	ngsFST Harvester
	_____________________________________________________________________________________________________
	
	A tool to summarize and filter ngsFST output files. It takes as an input either a single fst file
	(single-file mode), or a population matrix file (matrix mode). A population matrix is a 3-column
	tab-separated file: the first column contains the path to the pairwise Fst file, and the 2nd and 3rd
	columns the names of the populations involved.
	
	ngsFST files contain all sequenced positions, including non variable. In order to keep only variable
	sites, a threshold can be set with --pval on the probability of each site being variable. It can help
	to run --choose_p in order to test different probability limits on an individual file first (outputs
	the mean Fst for that filtering level, and the number of sites used).
	
	Averaging is done with non-overlapping sliding windows (width set by --window), which is the correct
	way, as opposed to averaging per-site Fst values.
	
	Retained sites can be output together with their position using the --print_sites option - be careful,
	the sites file must be exactly the one output by angsd (*.saf.pos.gz) and NOT the site selection file,
	which may be sorted in a different order. These sites can then be checked against a different SNP-call.
	
	It is aslo possible to provide a list a sites to be retained in the analysis, both in single-file and
	in matrix mode. Together with this file must be provided a full sites file (--sites). In that case again,
	the sites file must be the one output by angsd when generating the .saf files used in ngsFST. In matrix
	mode, all pairwise comparisons must share exactly the same list of loci.
	
	All runs (except choose_p) will dump a log file with the mean and standard deviation of the Fst. It is
	also possible to output the Fst estimates across populations in graphical mode with --print_means, and
	distribution of Fst across windows for each pair of population with --print_dist.
	_____________________________________________________________________________________________________
	
 
	Arguments:
	--fst=		one Fst outfile from ngsTools\' ngsFST program (single-file mode)
	--matrix=	population matrix (matrix mode - see details above)
	--out=		output file path and prefix [same as input]
	--pval=		threshold p-value for considering a site as variable [0.95]
	--window=	window size, in bp [10000]
	--list=		retain only a whitelisted loci for the analysis (requires the --sites option)
	
	Output options:
	--choose_p	print fst for different p-val thresholds (in single-file mode only)
	--print_dist	print the distribution of Fst across windows
	--print_means	print the mean and stdev for each pair of populations
	--print_sites	prints the sites above pval together with their position (in single-file mode only)
	--sites=	sites used to generate the Fst file (required for --print_sites)
	
	--help         	print this text.
	_____________________________________________________________________________________________________

	Example:
	
	./ngsFST_Harvester.R --fst=./pop1_pop2.fst --out=./pop1_pop2_filtered --pval=0.99 --window=20000
	--print_dist --print_sites --sites=./pop1.saf.pos
	_____________________________________________________________________________________________________
	

")	
  q(save="no")

}

###########################################################################################
#Standard argument parser

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
#The program takes either a single Fst file, or a matrix, but not both.

## Input error
if(is.null(argsL$fst) & is.null(argsL$matrix)){
	stop("You forgot to provide an input file")
} else if(is.null(argsL$fst)==F & is.null(argsL$matrix)==F){
	stop("You can not provide a single Fst file and a population matrix at the same time")
} else if(is.null(argsL$fst) & is.null(argsL$matrix)==F){	#will run in matrix mode
	run_mode <- 1
	argsL$matrix -> path_pops
	argsL$matrix -> path_out.tmp
} else if(is.null(argsL$fst)==F & is.null(argsL$matrix)){	#will run in single-file mode
	run_mode <- 2
	argsL$fst -> path_fst
	argsL$fst -> path_out.tmp
}

## --out default
if(is.null(argsL$out)) {
  paste(dirname(path_out.tmp), '/', basename(path_out.tmp), sep='', collapse='') -> path_out
} else { argsL$out ->  path_out }

## --choose_p incompatibility
if(is.null(argsL$choose_p)==F & is.null(argsL$matrix)==F){
	stop("--choose_p option only runs in single-file mode")}

## --choose_p default
if(is.null(argsL$choose_p)) {
  choose_p <- 0
} else { choose_p <- 1 }

## --pval default
if(is.null(argsL$pval)) {
  pval <- 0.95
} else { as.numeric(argsL$pval) -> pval }

if(is.null(argsL$choose_p)==F & is.null(argsL$pval)==F){
	cat("--pval argument is ignored in --choose_p mode.\n")}

## --window default
if(is.null(argsL$window)) {
  window.size <- 10000
} else { as.numeric(argsL$window) -> window.size }

## --print_means default
if(is.null(argsL$print_means)) {
  print_means <- 0
} else { print_means <- 1 }

## --print_means incompatibility
if(is.null(argsL$print_means)==F & run_mode==2) {
	cat("--print_means argument is ignored in single-file mode.\n")}

## --print_dist default
if(is.null(argsL$print_dist)) {
  print_dist <- 0
} else { print_dist <- 1 }

##Disentangle the varisou sites options:

## --sites default
if(is.null(argsL$sites)){
	sites_mode <- 0
} else { sites_mode <- 1 }

## --list default
if(is.null(argsL$list)){
	list_mode <- 0
} else { list_mode <- 1 }

## --print_sites default
if(is.null(argsL$print_sites)){
	print_sites <- 0
} else { print_sites <- 1 }

##Check the compatibilisty of options
## --sites --print_sites = print_sites mode 1 (emit the positions for pval sites)
## --sites --list = list_mode 1 (do the calculations on the selected sites)
## --sites --list --print_sites = print_sites mode 2 (print data for sites in list)
## --list --print sites : error
## --print_sites : error
## --list : error
## --sites : ignored

if(sites_mode==1 & print_sites==1 & list_mode==0){
	print_sites <- 1
} else if (sites_mode==1 & print_sites==0 & list_mode==1){
	list_mode <- 1
} else if (sites_mode==1 & print_sites==1 & list_mode==1){
	print_sites <- 2
} else if (sites_mode==0 & print_sites==1 & list_mode==1){
	stop("A site file is required for the --print_sites option")
} else if (sites_mode==0 & print_sites==0 & list_mode==1){
	stop("A site file is required for the --list option")
} else if (sites_mode==0 & print_sites==1 & list_mode==0){
	stop("A site file is required for the --print_sites option")
} else if (sites_mode==1 & print_sites==0 & list_mode==0){
	sites_mode <- 0
}

if(sites_mode==1){argsL$sites -> path_sites}
if(list_mode==1){argsL$list ->path_list}

if(list_mode==1 & is.null(argsL$pval)==F){
	cat("--pval argument is ignored in --list mode.\n")}


############################################################################################################################
#Print a recap of the chose options
	cat('\n')

	if(run_mode==1){cat("\tRunning in population matrix mode.\n")
	} else if(run_mode==2){cat("\tRunning in single-file mode.\n")}
	
	cat('\tInput file:\t\t', basename(path_out.tmp),'\n', sep='')
	
	cat('\tOutput file:\t\t', basename(path_out),'\n', sep='')
	
	if(is.null(argsL$window)) {
  			cat('\tWindow size:\t\t10000bp\n')
	} else {cat('\tWindow size:\t\t',window.size,'bp\n', sep='')}

	if(is.null(argsL$pval) & list_mode==0) {
  			cat('\tSite prob. limit:\t0.95\n')
	} else if (is.null(argsL$pval)==F & list_mode==0){
			cat('\tSite prob. limit:\t',pval,'\n', sep='')}
			
	if(list_mode==1){cat('\tSite filtering after:\t', basename(path_list),'\n', sep='')}

	if(sites_mode==1){cat('\tSites file:\t\t', basename(path_sites), '\n', sep='')}

	if(print_sites==1){cat('\tWill print the positions and Fst for sites passing pval threshold.\n')}
	
	if(print_sites==2){cat('\tWill print the positions and Fst for sites included in the list file.\n')}

	if(is.null(argsL$choose_p)==F) {cat('\t--choose_p mode: will only output Fst~site-probability distribution.\n')}

	if(is.null(argsL$print_dist)==F){cat('\tWill print the Fst distribution across windows.\n')}
	
	if(is.null(argsL$print_means)==F){cat('\tWill print the mean pairwise Fst across populations.\n')}


############################################################################################################################

suppressPackageStartupMessages(suppressWarnings(library(sqldf)))
suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(gridExtra)))
#suppressPackageStartupMessages(suppressWarnings(library(extrafonts)))


############################################################################################################################
##	IF USING A SITES FILE (sites_mode ==1)
############################################################################################################################
if(sites_mode==1){

#Reading in the sites file
cat('\n\tReading sites file...')
	data <- file(eval(path_sites))
	sqldf("select * from data", file.format = list(header=FALSE, sep='\t', comment.char=''))->sites
	closeAllConnections()
	cat('done.\n')

		####################################################################################################################
		##IF IN LIST MODE (list_mode==1)
		#We need to import the sites file, and the list file, and to build a flag vector.
		#This vector will be substituted to pval in the filtering. pval threshold becomes 1
		
		if(list_mode==1){
		
		#Reading in the list file
		cat('\tReading list file...')
			data <- file(eval(path_list))
			sqldf("select * from data", file.format = list(header=FALSE, sep='\t', comment.char=''))->list
			closeAllConnections()
			cat('done.\n')
		
		cbind(sites, V3=seq(1,nrow(sites),1))->sites
		merge(sites, list, by=c("V1","V2"), sort=F)->list.merged
		list.merged$V3->index
		
		list.flag<-vector(length=nrow(sites))	#a vector to bind with the fst files
		list.flag[]<-0							#sites not in list receive '0'
		
		###################
		list.flag[index]<-1						#sites in list receive '1'
		###################
		
		sum(list.flag)->n.snps
		nrow(list)->l.list
		
		cat("\tThere are ", l.list, " sites in the list, ", n.snps, " of which were found in the Fst file.\n", sep="")
		
		} #end for the flag list building
		####################################################################################################################
	
} # End of the sites-file reading task.
############################################################################################################################


############################################################################################################################
##	IF RUNNING IN MATRIX MODE
############################################################################################################################
	if(run_mode == 1){
	
############################################################################################################################
#Read population matrix

#Check whether the matrix file exist
if(file.exists(eval(path_pops))==FALSE){
	stop("Population matrix file not found.")}
	
read.table(eval(path_pops), sep='\t', header=F)->pops
levels(as.factor(c(as.character(pops$V2),as.character(pops$V3))))->pops.names
length(unique(as.numeric(as.factor(c(as.character(pops$V2),as.character(pops$V3))))))->n.pops
matrix(nrow=n.pops, ncol=n.pops)->res
rownames(res)<-pops.names
colnames(res)<-pops.names

############################################################################################################################
#Fix the population names

as.data.frame(unique(cbind(as.numeric(as.factor(c(as.character(pops$V2),as.character(pops$V3)))),c(as.character(pops$V2),as.character(pops$V3)))))->code
as.data.frame(cbind(as.character(pops$V2), as.character(pops$V3)))->popcodes
for(i in 1:n.pops){
as.character(popcodes$V1)->popcodes$V1
as.character(popcodes$V2)->popcodes$V2
as.character(code$V2[i])->p
as.character(code$V1[i])->q
popcodes[popcodes==p]<-q
}
data.frame(pops$V1, popcodes$V1, popcodes$V2)->pops
names(pops)<-c('file','pop1','pop2')

############################################################################################################################
#Get file list and check existence of files

as.character(pops$file)->filelist
length(filelist)->n.files

#Check whether any of the Fst files are missing
cat('\n')
for(f in 1:n.files){
cat('Checking ', filelist[f], '\n')
	if(file.exists(filelist[f])==FALSE){
	paste("Fst file ", filelist[f], " could not be found.", sep="", collapse="")->message
	stop(message)}
}

cat('\nFST VALUES:\n')
paste('\nFST VALUES:\n')->end_text

############################################################################################################################
#Set a list to receive the Fst distributions
	list()->fst.list
	vector()->means
	vector()->devs
	vector()->comps

############################################################################################################################
#Loop through the files
for(f in 1:n.files){

#Set current filename
	filelist[f]->path_current
	
#Set current population names
	pops[f,2:3]->coord
	as.numeric(as.character(coord$pop1))->x
	as.numeric(as.character(coord$pop2))->y
	pops.names[x]->pop1
	pops.names[y]->pop2

#SQLdf file import
	data <- file(eval(path_current))
	sqldf("select * from data", file.format = list(header=FALSE, sep='\t', comment.char=''))->fst.raw
	closeAllConnections()
	names(fst.raw)<-c('A','AB','f','fst','pval')
	
		####################################################################################################################
		#In case we are in list mode, we need to modify the dataframe and parameters.
		if(list_mode==1){
			cbind(fst.raw, fst.raw$pval)->fst.raw
			list.flag->fst.raw$pval
			names(fst.raw)<-c('A','AB','f','fst','pval', 'pval.true')
			pval<-1 #This is not a true pval anymore - maybe this could be made more rigorous at some point.
		}
		####################################################################################################################

############################################################################################################################
#Set up the sliding window
	nrow(fst.raw)->n.sites
	floor(n.sites/window.size)->n.window
	seq(1,window.size*n.window, window.size)->window.start
	seq(window.size,window.size*n.window, window.size)->window.end
	as.data.frame(cbind(window.start, window.end))->window.pos
	names(window.pos)<-c('start','end')

#Set the result vectors
	fst.dist <- vector(length=n.window)
	pval.dist <- vector(length=n.window)
	n.variable <- 0 #set counter for variable sites
	
############################################################################################################################	
#Loop throught the successive windows
	for(w in 1:n.window){

	window.pos[w,1]->start
	window.pos[w,2]->end
	
	fst.raw[start:end,]->fst.chunk
	fst.chunk[fst.chunk$pval>=pval,]->fst.window
	nrow(fst.window)->n.variable.window				#counting the number of retained sites
	n.variable + n.variable.window -> n.variable	#incrementing the total number

	sum(fst.window$A)->As
	sum(fst.window$AB)->AsBs
	(As/AsBs)->fst.window.mean
	fst.window.mean->fst.dist[w]
	
		####################################################################################################################
		if(list_mode==1){
		mean(fst.window$pval.true, na.rm=T)->pval.window.mean
		pval.window.mean->pval.dist[w]
		}
		####################################################################################################################
		
} #end of the sliding window loop
############################################################################################################################

fst.list[[f]]<-fst.dist #add the current loop's output to the result list for all files
signif(mean(fst.dist, na.rm=T),6)->fst.mean		#current file's mean over all windows
signif(sqrt(var(fst.dist, na.rm=T)),6)->fst.dev	#current file's dev over all windows

		####################################################################################################################
		if(list_mode==1){
		signif(mean(pval.dist, na.rm=T),3)->pval.mean
		} #current file's pval on selected sites over all windows
		####################################################################################################################

c(means, fst.mean)->means	#add to the pairwise statistics
c(devs, fst.dev)->devs		#add to the pairwise statistics
paste(pop1, '~', pop2, sep='', collpase='')->comp.name
c(comps, comp.name)->comps

	####################################################################################################################	
	#Print these statistics to the stdout and to a table
	if(list_mode==1){
		#Print the values in the frame-shaped block
		paste(n.variable, ' sites (', signif((n.variable/n.window),2), ' per window).', sep='', collapse='')-> variable.str
		cat(pop1, '::', pop2, '\t|  ', fst.mean, '\t± ', fst.dev, '\t|  mean site prob. = ', pval.mean,'\t|  ', variable.str, '\n', sep='') #print to standard out
		paste(pop1, '::', pop2, '\t|\t', fst.mean, '\t± ', fst.dev, '\t|\tmean site probability = ', pval.mean, sep='')->running_line	#print to frame
		rbind(end_text, running_line)->end_text
		
	} else if (list_mode==0) {
		#Print the values in the frame-shaped block
		paste(n.variable, ' sites (', signif((n.variable/n.window),2), ' per window).', sep='', collapse='')-> variable.str
		cat(pop1, '::', pop2, '\t|\t', fst.mean, '\t± ', fst.dev,'\t|  ', variable.str, '\n', sep='') #print to standard out
		paste(pop1, '::', pop2, '\t|\t', fst.mean, '\t± ', fst.dev, sep='')->running_line	#print to frame
		rbind(end_text, running_line)->end_text
	}
	####################################################################################################################

#Print the values in the matrix-shaped block
fst.mean->res[y,x]	#upper triangle
fst.dev->res[x,y] 	#lower triangle

############################################################################################################################
#PRINT_SITES option has two modes: [1] print sites with a pval threshold, [2] print sites in the list.

if(print_sites==1){
	
	eval(paste(path_out, '_', pop1, '_', pop2, '.pos', sep='', collapse=''))->pos_out
	
	nrow(sites)->n.pos
	nrow(fst.raw)->n.sites
	
		####################################################################################################################		
		#Abort loop if the length of the two files do not match.
		if(n.pos != n.sites){stop('\tNumber of sites in Fst files and in Sites files does not match. Will not print positions.')
		}
		####################################################################################################################	

	names(sites)<-c('CHROM','POS')
	data.frame(sites$CHROM, sites$POS, fst.raw$A, fst.raw$AB, fst.raw$f, fst.raw$A, fst.raw$pval)->data.raw
	names(data.raw)<-c('CHROM','POS','A','AB','f','FST','PVAL')
	data.raw[data.raw$PVAL>=pval,]->data.pval
	nrow(data.pval)->n.kept
	write.table(data.pval, eval(pos_out), sep='\t', row.names=F, col.names=T, quote=F)

	########################################################################################################################

} else if (print_sites==2){
	
	#list.flag was defined above as a vector with 1 for sites in the list, and 0 elsewhere

	eval(paste(path_out, '_', pop1, '_', pop2, '.sites', sep='', collapse=''))->pos_out
	
	nrow(sites)->n.pos
	nrow(fst.raw)->n.sites
	
		####################################################################################################################		
		#Abort loop if the length of the two files do not match.
		if(n.pos != n.sites){stop('\tNumber of sites in Fst files and in Sites files does not match. Will not print positions.')
		}
		####################################################################################################################	

	names(sites)<-c('CHROM','POS')
	data.frame(sites$CHROM, sites$POS, fst.raw$A, fst.raw$AB, fst.raw$f, fst.raw$A, fst.raw$pval)->data.raw
	names(data.raw)<-c('CHROM','POS','A','AB','f','FST','PVAL')
	
	data.raw[list.flag==1,]->data.list
	nrow(data.list)->n.kept
	write.table(data.list, eval(pos_out), sep='\t', row.names=F, col.names=T, quote=F)

}#end of the print_sites task.

############################################################################################################################
} #end of the file loop

names(fst.list)<-comps	#We can make the result list a bit nicer with names.

############################################################################################################################
#PRINT LOG FILE

	paste('PAIRWISE FST COMPARISONS, SITES WITH A PROBABILITY OF BEING VARIABLE >= ', pval, '\nWINDOW SIZE ', window.size, 'bp.\nUpper triangle: mean values, Lower triangle: standard deviation.\n', sep='', collapse='')->header
	eval(paste(path_out, '.log', sep='', collapse=''))->log_out
	write.table(header, eval(log_out), quote=F, row.names=F, col.names=F, append=F)
	write.table(res, eval(log_out), sep='\t', na='-', quote=F, col.names=NA, append=T)
	write.table(end_text, eval(log_out), quote=F, row.names=F, col.names=F, append=T)


############################################################################################################################
#PRINT_MEANS (Print the distribution of means)

	if(print_means==1){
	
	eval(paste(path_out, '.means.pdf', sep='', collapse=''))->pdf_out
	eval(paste(path_out, '.means.summary', sep='', collapse=''))->table_out
	
	#Build a dataframe with the means and deviations
	data.frame(comps, means, devs)->summary
	
	#Write to ggplot
	ggplot(summary) -> s
		s + geom_point(aes(x=comps, y=means), size=5, alpha=0.75, col='firebrick') -> s
		s + geom_errorbar(aes(x=comps, ymax=means+devs, ymin=means-devs), size=0.25, width=0.25) -> s
		s + theme_bw() -> s
		s + labs(list(title='Pairwise Fst comparisons (mean and std dev.)', x='', y='Mean pairwise Fst'), vjust=1) -> s
		s + theme(plot.title=element_text(size=16,lineheight=.8,vjust=1.2,family="Courier"), axis.text.x=element_text(size=12,family="Courier",angle=90), axis.text.y=element_text(size=12,family="Courier"), axis.title.y=element_text(size=12,lineheight=.8,vjust=0.8, family="Courier"))
	ggsave(pdf_out)
	unlink("Rplots.pdf", force=TRUE)
	
	#Write to text file
	names(summary)<-c('PAIR','MEAN','STDEV')
	write.table(summary, eval(table_out), quote=F, sep='\t', row.names=F, col.names=T)
	
	}#end ot the print_means option

############################################################################################################################
#PRINT_DIST (Print the full distributions)

#The total number of comparisons is n.files, and the distributions of Fst values are in fst.list
	if(print_dist==1){
	require(gridExtra)
	
	eval(paste(path_out, '.dist.pdf', sep='', collapse=''))->pdf_out
	eval(paste(path_out, '.dist.summary', sep='', collapse=''))->table_out

	
	max(unlist(lapply(fst.list, function(x) x[which.max(abs(x))])))->xmax
	list()->plot.list
	mean(means, na.rm=T)->overall_mean
	
	#Produce a plot per pairwise comparison
	for(l in 1:n.files){
		fst.list[[l]]->dist
		mean(dist, na.rm=T)->local_mean
		title<-comps[l]
		data.frame(dist)->dist.df
	
		#Write to ggplot
		ggplot(dist.df) -> d
			d + geom_histogram(aes(x=dist,y=..density..), binwidth=0.0025, alpha=0.25, col='black') -> d
			d + scale_x_continuous(limits = c(0, xmax)) -> d
			d + geom_vline(xintercept=overall_mean, col='black', linetype='dashed', size=0.75) -> d
			d + geom_vline(xintercept=local_mean, col='firebrick', linetype='dashed', size=0.75) -> d
			d + theme_bw() -> d
			d + labs(list(title=title, x='', y=''), vjust=1) -> d
			d + theme(plot.title=element_text(size=16,lineheight=.8,vjust=1.2,family="Courier"), axis.text.x=element_text(size=12,family="Courier"), axis.text.y=element_text(size=12,family="Courier"), axis.title.y=element_text(size=12,lineheight=.8,vjust=0.8, family="Courier"))

		plot.list[[l]]<-d
	}

	#Plot the result in pdf
	pdf(eval(pdf_out), width=21, height=27)
	
		i<-1
		plot<-list()
		
		for (n in 1:n.files){
	  	plot[[i]]<-plot.list[[n]]
	  		if(i %% 6 == 0) {
	    		print(do.call(grid.arrange, plot))
	    		plot<-list()
	    		i<-0}
	  	i<-i+1}
	  	
	if (length(plot) != 0) {print(do.call(grid.arrange, plot))}
	dev.off()
	
	#Write a table with the corresponding data
	write.table(fst.list, eval(table_out), quote=F, row.names=F, sep='\t')
	
	}#end of the print_dist option

############################################################################################################################
} # END OF MATRIX MODE
############################################################################################################################



############################################################################################################################
##	IF RUNNING IN SINGLE-FILE MODE
############################################################################################################################
	if(run_mode==2){
############################################################################################################################
#Check whether the list file exist

if(file.exists(eval(path_fst))==FALSE){
	stop("Fst file not found.")}
	
############################################################################################################################
#SQLdf file import

	cat('\n\tReading per-site-Fst from ', basename(path_fst), '...', sep='')
	data <- file(eval(path_fst))
	sqldf("select * from data", file.format = list(header=FALSE, sep='\t', comment.char=''))->fst.raw
	closeAllConnections()
	cat('done.\n')
	names(fst.raw)<-c('A','AB','f','fst','pval')

############################################################################################################################
#PRINT_SITES option has two modes: [1] print sites with a pval threshold, [2] print sites in the list.

if(print_sites==1){
	
	eval(paste(path_out, '.pos', sep='', collapse=''))->pos_out
	nrow(sites)->n.pos
	nrow(fst.raw)->n.sites
	
		####################################################################################################################		
		#Abort loop if the length of the two files do not match.
		if(n.pos != n.sites){stop('\tNumber of sites in Fst files and in Sites files does not match. Will not print positions.')
		}
		####################################################################################################################	

	names(sites)<-c('CHROM','POS')
	data.frame(sites$CHROM, sites$POS, fst.raw$A, fst.raw$AB, fst.raw$f, fst.raw$A, fst.raw$pval)->data.raw
	names(data.raw)<-c('CHROM','POS','A','AB','f','FST','PVAL')
	data.raw[data.raw$PVAL>=pval,]->data.pval
	nrow(data.pval)->n.kept
	cat('\tKept ', n.kept, ' sites out of ', n.pos, '.\n', sep='')
	
	cat('\tWriting filtered positions...')
	write.table(data.pval, eval(pos_out), sep='\t', row.names=F, col.names=T, quote=F)
	cat('done.\n')

	########################################################################################################################

} else if (print_sites==2){
	
	#list.flag was defined above as a vector with 1 for sites in the list, and 0 elsewhere

	eval(paste(path_out, '.sites', sep='', collapse=''))->pos_out
	nrow(sites)->n.pos
	nrow(fst.raw)->n.sites
	
		####################################################################################################################		
		#Abort loop if the length of the two files do not match.
		if(n.pos != n.sites){stop('\tNumber of sites in Fst files and in Sites files does not match. Will not print positions.')
		}
		####################################################################################################################	

	names(sites)<-c('CHROM','POS')
	data.frame(sites$CHROM, sites$POS, fst.raw$A, fst.raw$AB, fst.raw$f, fst.raw$A, fst.raw$pval)->data.raw
	names(data.raw)<-c('CHROM','POS','A','AB','f','FST','PVAL')
	
	data.raw[list.flag==1,]->data.list
	nrow(data.list)->n.kept
	cat('\tKept ', n.kept, ' sites out of ', n.pos, '.\n', sep='')	
	cat('\tWriting filtered positions...')
	write.table(data.list, eval(pos_out), sep='\t', row.names=F, col.names=T, quote=F)
	cat('done.\n')

}#end of the print_sites task.

############################################################################################################################
#Resume the main task (sliding window calculation)

		####################################################################################################################
		#In case we are in list mode, we need to modify the dataframe and parameters.
		if(list_mode==1){
			cbind(fst.raw, fst.raw$pval)->fst.raw
			list.flag->fst.raw$pval
			names(fst.raw)<-c('A','AB','f','fst','pval', 'pval.true')
			pval<-1 #This is not a true pval anymore - maybe this could be made more rigorous at some point.
		}
		####################################################################################################################

#Set up the sliding window
	nrow(fst.raw)->n.sites
	floor(n.sites/window.size)->n.window
	seq(1,window.size*n.window, window.size)->window.start
	seq(window.size,window.size*n.window, window.size)->window.end
	as.data.frame(cbind(window.start, window.end))->window.pos
	names(window.pos)<-c('start','end')

if(choose_p==0){ #This conditions spans over log file writing, and print_dist option, as these don't make sense in choose_p mode.

#Set the result vectors
	fst.dist <- vector(length=n.window)
	pval.dist <- vector(length=n.window)
	
#Loop throught the successive windows
	cat('\tPerforming sliding window averaging...')
	for(w in 1:n.window){

	window.pos[w,1]->start
	window.pos[w,2]->end
	
	fst.raw[start:end,]->fst.chunk
	fst.chunk[fst.chunk$pval>=pval,]->fst.window

	sum(fst.window$A)->As
	sum(fst.window$AB)->AsBs
	(As/AsBs)->fst.window.mean
	fst.window.mean->fst.dist[w]

		####################################################################################################################
		if(list_mode==1){
		mean(fst.window$pval.true, na.rm=T)->pval.window.mean
		pval.window.mean->pval.dist[w]
		}
		####################################################################################################################

	}#end of the window loop
	cat('done.\n')
############################################################################################################################

		####################################################################################################################
		if(list_mode==1){
		signif(mean(pval.dist, na.rm=T),3)->pval.mean
		} #current file's pval on selected sites over all windows
		####################################################################################################################

signif(mean(fst.dist, na.rm=T),6)->fst.mean
signif(sqrt(var(fst.dist, na.rm=T)),6)->fst.dev

	########################################################################################################################	
	#Print these statistics to the stdout and to a table
	if(list_mode==1){
		#Print the values in the frame-shaped block
		cat('\t___________________________________________________________________________________________\n')
		cat('\n\tMEAN FST:\t', fst.mean, '\t± ', fst.dev, '\t|   mean site probability = ', pval.mean,'\n', sep='') #print to standard out
		cat('\t___________________________________________________________________________________________\n\n')
		paste('MEAN FST:\t', fst.mean, '\t± ', fst.dev, '\t|   mean site probability = ', pval.mean,'\n', sep='')->res #print to string
	} else if (list_mode==0) {
		#Print the values in the frame-shaped block
		cat('\t____________________________________________________\n')
		cat('\n\tMEAN FST:\t', fst.mean, '\t± ', fst.dev, '\n', sep='') #print to standard out
		cat('\t____________________________________________________\n\n')
		paste('MEAN FST:\t', fst.mean, '\t± ', fst.dev, '\n', sep='')->res #print to string
	}
	########################################################################################################################

############################################################################################################################
#PRINT LOG FILE

	paste('MEAN WINDOWED FST FOR ', path_fst ,'\nPROBABILITY OF SITES BEING VARIABLE >= ', pval, '\nWINDOW SIZE ', window.size, 'bp\n', sep='', collapse='')->header
	eval(paste(path_out, '.log', sep='', collapse=''))->log_out
	write.table(header, eval(log_out), quote=F, row.names=F, col.names=F, append=F)
	write.table(res, eval(log_out), sep='\t', quote=F, row.names=F, col.names=F, append=T)

############################################################################################################################
#PRINT_DIST OPTION [single-file mode]

	if(print_dist==1){

	eval(paste(path_out, '.dist.pdf', sep='', collapse=''))->pdf_out
	eval(paste(path_out, '.dist.summary', sep='', collapse=''))->table_out

	title <- paste('Windowed Fst for ', basename(path_fst), ', ', window.size, 'bp windows', sep='', collapse='')
	data.frame(dist=fst.dist)->dist.df
	ggplot(dist.df) -> d
		d + geom_histogram(aes(x=dist,y=..density..), binwidth=0.0024, alpha=0.25, col='black') -> d
		d + geom_vline(xintercept=fst.mean, col='firebrick', linetype='dashed', size=0.75) -> d
		d + theme_bw() -> d
		d + labs(list(title=title, x='', y=''), vjust=1) -> d
		d + theme(plot.title=element_text(size=12,lineheight=.8,vjust=1.2,family="Courier"), axis.text.x=element_text(size=12,family="Courier"), axis.text.y=element_text(size=12,family="Courier"), axis.title.y=element_text(size=12,lineheight=.8,vjust=0.8, family="Courier"))

	ggsave(pdf_out)
	unlink("Rplots.pdf", force=TRUE)
	write.table(fst.dist, eval(table_out), quote=F, row.names=F, sep='\t')

	}#end of the print_dist option

############################################################################################################################

}#end ot the choose_p exception

############################################################################################################################
#CHOOSE_P OPTION : Calculates the distribution on Fst ~ pval and print outs the result

	if(choose_p==1){
	
	eval(paste(path_out, '.choose_p.pdf', sep='', collapse=''))->pdf_out
	eval(paste(path_out, '.choose_p.summary', sep='', collapse=''))->table_out

	p_thres <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,1)
	length(p_thres)->n.p
	vector()->means
	vector()->devs
	vector()->kept_sites

	for(p in 1:n.p){
		pval <- p_thres[p]
			cat('\tPerforming sliding window averaging for pval=', pval, '...', sep='')
			
		#Set the result vector	
		fst.dist <- vector(length=n.window)
	
		#Loop throught the successive windows
		for(w in 1:n.window){
			window.pos[w,1]->start
			window.pos[w,2]->end
			fst.raw[start:end,]->fst.chunk
			fst.chunk[fst.chunk$pval>=pval,]->fst.window
			sum(fst.window$A)->As
			sum(fst.window$AB)->AsBs
			(As/AsBs)->fst.window.mean
			fst.window.mean->fst.dist[w]
		}


signif(mean(fst.dist, na.rm=T),6)->fst.mean
signif(sqrt(var(fst.dist, na.rm=T)),6)->fst.dev
c(means, fst.mean)->means	#add to the pairwise statistics
c(devs, fst.dev)->devs		#add to the pairwise statistics

nrow(fst.raw[fst.raw$pval>=pval,])->n.kept
c(kept_sites, n.kept)->kept_sites

cat('retained ', n.kept, ' sites.\n', sep='')

} #End ot the P-threshold loop

	data.frame(p_thres, means, devs, kept_sites)->choose_p.summary
	names(choose_p.summary) <- c('PVAL','MEAN','STDEV','N_SITES')

	#Plot the distribution
	ggplot(choose_p.summary, aes(x=PVAL))->p
		p + geom_errorbar(aes(ymax=MEAN+STDEV, ymin=MEAN-STDEV), size=0.5, width=0) -> p
		p + geom_point(aes(y=MEAN), size=3, alpha=0.75, col='firebrick') -> p
		p + geom_line(aes(y=MEAN), size=0.25, col='firebrick', linetype='dashed') -> p
		p + geom_line(aes(y=MEAN+STDEV), size=0.25, linetype='dashed') -> p
		p + geom_line(aes(y=MEAN-STDEV), size=0.25, linetype='dashed') -> p
		p + labs(list(title='Windowed Fst ~ probability of a site being variable (pval)', x='Probability threshold', y='Mean pairwise Fst'), vjust=1) -> p
		p + theme(plot.title=element_text(size=12,lineheight=.8,vjust=1.2,family="Courier"), axis.text.x=element_text(size=12,family="Courier"), axis.text.y=element_text(size=12,family="Courier"), axis.title.x=element_text(size=12,lineheight=.8,vjust=0, family="Courier"), axis.title.y=element_text(size=12,lineheight=.8,vjust=0.8, family="Courier"))

ggsave(pdf_out)
unlink("Rplots.pdf", force=TRUE)
write.table(choose_p.summary, eval(table_out), quote=F, sep='\t', row.names=F, col.names=T)


} #End of the choose_p option


############################################################################################################################
} # END OF SINGLE-FILE MODE
############################################################################################################################

cat('\n\tΤΑ ΔΗ ΝΥΝ ΠΑΝΤΑ ΤΕΛΕΙΤΑΙ\n\n')
q(save="no")
	
