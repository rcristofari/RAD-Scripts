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
	_____________________________________________________________________________________________________

	Fst Tools
	_____________________________________________________________________________________________________

	Calculates different pairwise-Fst estimators based on a VCF file (--vcf mode). Alternatively, 
	it can simulate HWE genotype data for finite populations undergoing random drift (--sim mode),
	taking into account the uncertainty due to sequencing coverage.
	The simulator only handles equal population sizes for the moment. You need to specify the number
	of samples per population, but only one population size.
	_____________________________________________________________________________________________________

	Arguments:
	--vcf=				VCF file containing empirical genotypes
		--popmap=		a population map
	--sim				simulate genotypes for Fst calculation.
		--nPop=			number of populations to simulate [2]
		--nInd=			number of individuals in these populations [10,000]
		--nSamples=		a list of sample sizes per pop, coma-separated [24,24]
		--nSites=		number of sites to simulate [1000]
		--nGen=			number of generations of drift [1000]
		--depth=		mean sequencing depth [3]
		--missing=		proportion of missing genotypes [0.2]
	
	--out=				path and filename to write the output.
	--window=			Window size for Fst averaging [50]
	
	--help      			print this text
 
	Example:
	./FST_Tools.R
 	_____________________________________________________________________________________________________\n\n")
q(save="no")}

###########################################################################################
#	#STANDARD ARGUMENT PARSER
###########################################################################################

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

#Mode selection:
if(is.null(argsL$vcf) & is.null(argsL$sim)){
	stop('You have to provide a VCF file, or simulation parameters')
} else if (is.null(argsL$vcf)==F & is.null(argsL$sim)==F){
	stop('You have to choose between simulation and VCF mode')
} else if (is.null(argsL$vcf) & is.null(argsL$sim)==F){
	run_mode <- 1	#simulation mode
} else if (is.null(argsL$vcf)==F & is.null(argsL$sim)){
	run_mode <- 2	#vcf mode
	argsL$vcf -> path_vcf
}

#--popmap setting:
if(run_mode==2 & is.null(argsL$popmap)){
	stop('You must provide a population map')
} else if (run_mode==2 & is.null(argsL$popmap)==F){
	argsL$popmap ->path_popmap
}

#--out default:
if(is.null(argsL$out)){
	paste(getwd(), 'out.fst', sep='', collapse='')->path_out
} else {
	argsL$out -> path_out
}

#Fst calculation default:
#if(is.null(argsL$reich) & is.null(argsL$reynolds) & is.null(argsL$weir)){
#	fst_mode <- 1 #defaults to Reich's Fst
#} else if (is.null(argsL$reich)==F & is.null(argsL$reynolds) & is.null(argsL$weir)){
#	fst_mode <- 1 #defaults to Reich's Fst
#} else if (is.null(argsL$reich) & is.null(argsL$reynolds)==F & is.null(argsL$weir)){
#	fst_mode <- 2 #defaults to Reynolds' Fst
#} else if (is.null(argsL$reich) & is.null(argsL$reynolds) & is.null(argsL$weir)==F){
#	fst_mode <- 3 #defaults to Weir's Fst
#} else {
#	stop('You can not calculate several estimators in a single run')
#}

#Have to tweak that to the program calculates and outputs all three.

#--nPop	default:
if(is.null(argsL$nPop)){
	n.pop <- 2	
} else {
	as.numeric(argsL$nPop) -> n.pop
}

#--nInd	default:
if(is.null(argsL$nInd)){
	n.ind <- 10000
} else {
	as.numeric(argsL$nInd) -> n.ind
}

#--nSamples default:
if(is.null(argsL$nSamples)){
	n.samples <- c(24,24)
} else {
	as.numeric(unlist(strsplit(argsL$nSamples, ","))) -> n.samples
}

#--nSites default:
if(is.null(argsL$nSites)){
	n.sites <- 1000
} else {
	as.numeric(argsL$nSites) -> n.sites
}

#--nGen	default:
if(is.null(argsL$nGen)){
	n.generations <- 1000
} else {
	as.numeric(argsL$nGen) -> n.generations
}

#--depth default:
if(is.null(argsL$depth)){
	d <- 3
} else {
	as.numeric(argsL$depth) -> d
}

#--missing default:
if(is.null(argsL$missing)){
	missing <- 0.2
} else {
	as.numeric(argsL$missing) -> missing
}

#Window size default:
if(is.null(argsL$window)){
	window.size <- 50
} else { 
	as.numeric(argsL$window) -> window.size
}

#########################################################################################################################
#SIMULATION BLOCK (neutral drift between n pops, streams out a geno.list object)
#########################################################################################################################
if(run_mode==1){

#Simulate the ancestral population allele frequencies, drawn from a (log)normal distribution.
#Create n populations of n1 and n2 diploids, and S sites
#vector of maf in pop1 and in pop2 (2*n1+1) and (2*n2+1)
#rbinom from each saf, converted to a maf, and sent to the cell of the res vector.
#Iterated over as many generations as necessary. Inside loop to treat both pops.
#At the end of the process, convert the safs to genotypes: sample(c(0,1), 2*n, replace=T, prob=c((1-saf),saf)) to get a genotype realization of that saf.
#Maybe rather needs to make HWE genotypes ? saf^2 2*saf*(1-saf) and (1-saf)^2 ?
#012 genotype matrix of nind*nsites
#Create a mtarix of sequencing depth per pop. cols=ind, lines=sites
#Each col averages up to DEPTH, but randomly drawn around it (allowing it to be 0)

###########################################################################
#Required arguments:
#n.pop			<-	4
#n.ind			<-	100
#n.sites			<-	50000
#n.generations 	<-	3
#d				<-	3
#Need to implement unequal population sizes (for asymmetrical Fst), and unequal sample sizes.
#n.pops <- c(A,B,C,D) not yet. Might slow the whole thing down too much, unless I find a clever way to store it in the list for drift.step.
#n.samples <- c(25,25,15,4)
###########################################################################

#A starting distribution for minor allele frequencies could be drawn from the normal distribution (or other).
freq.weights <- dnorm(seq(0,1,0.001), 0.1, 0.1)
starting.freqs <- sample(seq(0,1,0.001), n.sites, replace=T, prob=freq.weights)
1-starting.freqs -> starting.freqs

#Each population is represented by a vector of length n.sites, stored as an element of the pop.list list
pop.list <- list()

#At the starting state, all populations share the same allele frequencies:
for(pop in 1:n.pop){
	starting.freqs -> pop.list[[pop]]
}

#Define functions to apply the drift process to the whole population list:
#One drift step: re-evaluation of the allele frequency as a redraw from the binomial distribution.
	drift.step <- function(f){ rbinom(1, size=n.ind, prob=f) / n.ind }
	apply.drift <- function(p){ sapply(p, drift.step) }

#We will store all generations in one superlist:
pop.time.list <- list()
pop.list->pop.time.list[[1]]

cat('\n\tSimulating drift over ', n.generations, ' generations...\n', sep='', collapse='')
pb <- txtProgressBar(min = 0, max = n.generations, style = 3)		#set the progress bar

#Iterate over the number of generations:
for(t in 1:n.generations){
	lapply(pop.time.list[[t]], apply.drift)->pop.time.list[[t+1]]
	setTxtProgressBar(pb,t)
}

#The final states are in pop.time.list[[n.generations]]. They will be used to populate genotypes matrices.



#Generate genotypes from the drifted allele frequencies, and populate a matrix per pop, in geno.list (individuals in rows, sites in columns):
geno.list <- list()
for(p in 1:n.pop){

	#First we define a couple functions to assign genotypes in 012 format given a known minor allele frequency:
	#A function to draw hardy-weinberg weigths from a minor allele frequency p:
	hwe.weights <- function(x){ c((x^2),(2*x*(1-x)),((1-x)^2)) }
	#A function to draw n genotypes from that probability distribution:
	draw.genotype <- function(x){ sample(c(0,1,2), n.samples[p], replace=T, prob=hwe.weights(x))}
	sapply(pop.time.list[[n.generations]][[p]], draw.genotype)->geno.list[[p]]	
}

#And the same for the sequencing depth: here the depth is drawn from a gamma distribution of rate 1 and mean d
depth.list <- list()
for(p in 1:n.pop){
	#Sequencing depth has a 'gamma+invariant' distribution
	sequencing <- sample(c(ceiling(rgamma(((1-missing)*(n.sites*n.samples[p])),d,1)), rep(0,(missing*(n.sites*n.samples[p])))))
	d/mean(sequencing)->correction
	sequencing*correction->sequencing # Now the depth is exactly d, including invariant sites
	matrix(nrow=n.samples[p], ncol=n.sites, data=sequencing)->depth.list[[p]]	
}

#Now we define a function to simulate the sampling error due to low sequencing coverage:
#g is the 012 genotype, and d the sequencing depth.
seq.prob <- function(g,d){ rbind(c(1,0,0),c((2^(-d)),(1-2^(1-d)),(2^(-d))), c(0,0,1))[(g+1),] }
seq.sample <- function(g,d){ if(d>0){sample(c(0,1,2), 1, prob=seq.prob(g,d))} else {return(NA)} }

#Now we have to go into full-blown looping until something better turns up...
#Create a list to gather the results:
seq.list <- list()

cat('\n\tSimulating sequencing error:')

for(p in 1:n.pop){
	cat('\n\tPopulation ', p, '/', n.pop, '\n', sep='') 
	seq.list[[p]]<-matrix(nrow=n.samples[p], ncol=n.sites, data=-999)
	pb <- txtProgressBar(min = 0, max = n.samples[p], style = 3)		#set the progress bar
	for(i in 1:n.samples[p]){
		for(s in 1:n.sites){
			g<-geno.list[[p]][i,s]
			d<-depth.list[[p]][i,s]
			seq.sample(g,d)->seq.list[[p]][i,s]
		}	
	setTxtProgressBar(pb,i)
	}
}

names(seq.list)<-sprintf("POP%02d",1:n.pop)
geno.list <- list()
seq.list -> geno.list


} #End of simulation block
#########################################################################################################################
#VCF FILE IMPORT AND FORMAT BLOCK (streaming out a geno.list object)
#########################################################################################################################
if(run_mode==2){

#required arguments:
#path_vcf <- '~/Desktop/emperor_75.recode.vcf'
#path_popmap <- '~/Desktop/emperor.popmap'

cat("\n\tReading VCF file...")
paste("cat ", path_vcf, " | grep -v '##' > ", path_vcf, ".tmp", sep='', collapse='')->cmd
system(cmd)
paste(path_vcf, '.tmp', sep='', collapse='') -> path_tmp
read.table(eval(path_tmp), sep='\t', header=T, comment.char='')->vcf # en sqldf
paste("rm ", path_vcf, ".tmp", sep='', collapse='')->cmd
system(cmd)
cat("done.\n")

#Reset the number of sites as the observed value, except if a ceiling value is specified.
if(is.null(argsL$nSites)){
	n.sites <- nrow(vcf)
} else {
	as.numeric(argsL$nSites) -> n.sites
}

cat("\tReading population map...")
read.table(path_popmap, header=F, sep='\t')->popmap
ind.list<-split(popmap[,1], popmap[,2])
lapply(ind.list, as.character) -> ind.list
length(names(ind.list)) ->n.pop
cat("done.\n")

geno.list<-list()
#Create a list of the vectors of individuals in each pop as ind.list
#Then ind.list[[p]] is the list to retain for that pop

cat("\tParsing VCF and population map:\n")
for(p in 1:n.pop){
	
cat("\tPopulation ", p, "/", n.pop,":\n", sep='', collapse="")
#Retrieve vcf fields corresponding to that population
ind.list[[p]]->pop.ind
vcf[,pop.ind]->pop

#A complicated way to keep only the genotypes for each individual
as.list(pop)->pop.list
pop.split <- function(x){strsplit(as.character(x), ':')[[1]][1]}
pop.apply <- function(x){unlist(lapply(x, pop.split))}
lapply(pop.list, pop.apply)->gen.list

#Then convert them to 012 format
sub.<-function(x){gsub('./.', 'NA', x)}
sub0<-function(x){gsub('0/0','0',x)}
sub1<-function(x){gsub('0/1','1',x)}
sub2<-function(x){gsub('1/1','2',x)}

cat('\t|##        |\t20%')
lapply(gen.list, sub0)->gen.list
cat('\r\t|####      |\t40%')
lapply(gen.list, sub1)->gen.list
cat('\r\t|######    |\t60%')
lapply(gen.list, sub2)->gen.list
cat('\r\t|########  |\t80%')
lapply(gen.list, sub.)->gen.list
cat('\r\t|##########|\t100%\tLoading data...')
lapply(gen.list, as.numeric)->gen.list
#And transform that into a genotype matrix
do.call(cbind.data.frame, gen.list)->gen.df
as.matrix(gen.df)->geno
cat('done.\n')

#Populate the genotype list with these matrices
t(geno)->geno.list[[p]]

}
names(ind.list)->names(geno.list)


} #End of the VCF reading block
#########################################################################################################################
#FST CALCULATION BLOCK (takes a geno.list object, and performs sliding window Fst calculation)
#########################################################################################################################

cat("\n\tPerforming Fst estimation:\n")
##################################################################################################
#Compute Reich et al. 2009 Fst estimator:
#a is the number of observed minor alleles,
#n is the number of observed alleles (minor+major)
reich_fst <- function(ai,aj,ni,nj){
	hi <- (ai*(ni-ai))/(ni*(ni-1))
	hj <- (aj*(nj-aj))/(nj*(nj-1))
	#Estimator of N:
	N <- (((ai/ni)-(aj/nj))^2 - (hi/ni) - (hj/nj))
	#Estimator for D:
	D <- N + hi + hj
	#Estimator of Fst : S(N) / S(D) over all markers (or windowed)
	#Fst <- N/D
	return(c(N,D))
}
	#minor allele count a1 is the sum of the site row,
	#total allele count n1 is the 2*number of non-NA genotypes
	
##################################################################################################
#Compute Reynolds, Weir and Cockerham 1983 Fst estimator:
reynolds_fst <- function(ai,aj,ni,nj){
	#in this formula, n is the number of individuals, not of chromosomes, so me must divide by two the result of get.n()
	ni <- ni/2
	nj <- nj/2
	pi <- ai/(2*ni)
	pj <- aj/(2*nj)
	ps <- (ai+aj)/(2*(ni+nj))
	alphai <- 2*pi*(1-pi)
	alphaj <- 2*pj*(1-pj)
	bs <- (ni*alphai + nj*alphaj) / (ni+nj-1)
	as <- ((4*ni*(pi-ps)^2 + 4*nj*(pj-ps)^2) - bs) / (2*((2*ni*nj)/(ni+nj)))
	N <- as
	D <- as+bs
	return(c(N,D))
}

##################################################################################################
#Compute Weir & Cockerham 1984 Fst estimator:     #there's a problem in weir and reynolds leading to underestimation, check formulas !!!
weir_fst <- function(ai,aj,ni,nj){
	ni <- ni/2
	nj <- nj/2
	pi <- ai/(ni*2)
	pj <- aj/(nj*2)
	r<-2
	n <- (ni+nj)/2
	p <- (ai+aj)/(2*((ni+nj)))
	s2 <- ((ni*((pi-p)^2))+(nj*((pj-p)^2)))/n
	hi <- (ai*(ni-ai))/(ni*(ni-1))
	hj <- (aj*(nj-aj))/(nj*(nj-1))
 	h <- ((ni*hi)+(nj*hj))/(2*n)
	nc <- 2*n - ((ni^2 + nj^2)/(2*n))
	A <- (n/nc) * (s2 - (1/(n-1)) * ((p*(1-p)) - (s2/2) - (h/4)))
	B <- (n/(n-1)) * ((p*(1-p)) - (s2/2) - (((2*n-1)/(4*n))*h))
	C <- h/2
	N<-A
	D<-(A+B+C)
	return(c(N,D))
}


##################################################################################################
#Compute Wright 1951 Fst estimator:
#wright.fst <-function(ai,aj,ni,nj){
#	pi <- ai/ni
#	pj <- aj/nj
#	ps <- (ai+aj)/(ni+nj)
#	s2 <- ((((pi-ps)^2)*ni) + (((pj-ps)^2)*pj))/((ni+nj)/2)
#	Fst <- s^2 / (ps*(1-ps))
#	N<-Fst
#	D<-1
#	return(c(N,D))
#}

get.a <- function(x){ sum(x, na.rm=T) }
get.n <- function(x){ 2*sum(!is.na(x)) }
#This may be useful to use empirical heterozygosity instead of a HWE estimator..?
get.h <- function(x){ sum(x==1, na.rm=T) / (sum(!is.na(x))*2) }

##################################################################################################
#Create a result matrix
reich.fst <- matrix(nrow=n.pop, ncol=n.pop, data=-9)
row.names(reich.fst) <- names(geno.list)
colnames(reich.fst) <- names(geno.list)

reynolds.fst <- matrix(nrow=n.pop, ncol=n.pop, data=-9)
row.names(reynolds.fst) <- names(geno.list)
colnames(reynolds.fst) <- names(geno.list)

weir.fst <- matrix(nrow=n.pop, ncol=n.pop, data=-9)
row.names(weir.fst) <- names(geno.list)
colnames(weir.fst) <- names(geno.list)

reich.res<-list()
reynolds.res<-list()
weir.res<-list()

k<-0
pb <- txtProgressBar(min = 0, max = sum(seq(1,n.pop-1, 1)), style = 3)		#set the progress bar

for(p in 1:n.pop){
		geno.list[[p]]->pop1

	r<-p+1
	if(r>n.pop){break}

	for(q in r:n.pop){
		geno.list[[q]]->pop2
		
		reich.ND <- vector()
		reynolds.ND <- vector()
		weir.ND <- vector()					
		
		for(s in 1:n.sites){
		
				reich_fst(get.a(pop1[,s]), get.a(pop2[,s]), get.n(pop1[,s]), get.n(pop2[,s]))->reich
				reynolds_fst(get.a(pop1[,s]), get.a(pop2[,s]), get.n(pop1[,s]), get.n(pop2[,s]))->reynolds
				weir_fst(get.a(pop1[,s]), get.a(pop2[,s]), get.n(pop1[,s]), get.n(pop2[,s]))->weir

				rbind(reich.ND,reich)->reich.ND
				rbind(reynolds.ND,reynolds)->reynolds.ND
				rbind(weir.ND,weir)->weir.ND
				

		}

	k+1->k
	reich.ND->reich.res[[k]]
	reynolds.ND->reynolds.res[[k]]
	weir.ND->weir.res[[k]]
	setTxtProgressBar(pb,k)
	}

}	#The resulting Fst components are now stored in 3 lists, flattened by k index.


##################################################################################################
#Now, we need a sliding window algorithm to extract the Fst averages (it's the same as in ngsFST_Harvester.R)

	
	
	reich.list <- list()	#List to store the windowed fst means	
	reynolds.list <- list()
	weir.list <- list()
	
	reich.means <- vector()
	reynolds.means <- vector()
	weir.means <- vector()
	
	reich.devs <- vector()
	reynolds.devs <- vector()
	weir.devs <- vector()


for(f in 1:k){	#Looping through the k pairwise comparisons
	as.data.frame(reich.res[[f]], row.names=F)->reich.raw
	as.data.frame(reynolds.res[[f]], row.names=F)->reynolds.raw
	as.data.frame(weir.res[[f]], row.names=F)->weir.raw
	
	names(reich.raw)<-c('N','D')
	names(reynolds.raw)<-c('N','D')
	names(weir.raw)<-c('N','D')
	
	#Set up the sliding window
	nrow(reich.raw)->n.sites
	floor(n.sites/window.size)->n.window
	seq(1,window.size*n.window, window.size)->window.start
	seq(window.size,window.size*n.window, window.size)->window.end
	as.data.frame(cbind(window.start, window.end))->window.pos
	names(window.pos)<-c('start','end')
	
	#Set the result vector	
	reich.dist <- vector(length=n.window)
	reynolds.dist <- vector(length=n.window)
	weir.dist <- vector(length=n.window)
	
	#Loop throught the successive windows
	for(w in 1:n.window){
	window.pos[w,1]->start
	window.pos[w,2]->end
	
	reich.raw[start:end,]->reich.window
	sum(reich.window$N, na.rm=T)->N
	sum(reich.window$D, na.rm=T)->D
	(N/D)->reich.window.mean
	reich.window.mean->reich.dist[w]
	
	reynolds.raw[start:end,]->reynolds.window
	sum(reynolds.window$N, na.rm=T)->N
	sum(reynolds.window$D, na.rm=T)->D
	(N/D)->reynolds.window.mean
	reynolds.window.mean->reynolds.dist[w]
	
	weir.raw[start:end,]->weir.window
	sum(weir.window$N, na.rm=T)->N
	sum(weir.window$D, na.rm=T)->D
	(N/D)->weir.window.mean
	weir.window.mean->weir.dist[w]

	}

reich.list[[f]]<-reich.dist
signif(mean(reich.dist, na.rm=T),6)->reich.mean
signif(sqrt(var(reich.dist, na.rm=T)),6)->reich.dev
c(reich.means, reich.mean)->reich.means
c(reich.devs, reich.dev)->reich.devs


reynolds.list[[f]]<-reynolds.dist
signif(mean(reynolds.dist, na.rm=T),6)->reynolds.mean
signif(sqrt(var(reynolds.dist, na.rm=T)),6)->reynolds.dev
c(reynolds.means, reynolds.mean)->reynolds.means
c(reynolds.devs, reynolds.dev)->reynolds.devs

weir.list[[f]]<-weir.dist
signif(mean(weir.dist, na.rm=T),6)->weir.mean
signif(sqrt(var(weir.dist, na.rm=T)),6)->weir.dev
c(weir.means, weir.mean)->weir.means
c(weir.devs, weir.dev)->weir.devs

} #Finished looping through the k pairs


k<-0
for(p in 1:n.pop){
	r <- p+1
	
	NA->reich.fst[p,p]
	NA->reynolds.fst[p,p]
	NA->weir.fst[p,p]
	
	if(r>n.pop){break}
	for(q in r:n.pop){
	k<-k+1
	
	reich.fst[p,q]<-reich.means[k]
	reich.fst[q,p]<-reich.devs[k]

	reynolds.fst[p,q]<-reynolds.means[k]
	reynolds.fst[q,p]<-reynolds.devs[k]

	weir.fst[p,q]<-weir.means[k]
	weir.fst[q,p]<-weir.devs[k]

	}
}


sink(path_out)
cat("\n\nRESULTS: FST estimates in upper triangle, STDEV in lower triangle")
cat("\n\nReich's Fst:\n\n")
print(reich.fst)
cat("\nReynolds' Fst:\n\n")
print(reynolds.fst)
cat("\nWeir's Fst:\n\n")
print(weir.fst)
cat("\n\n")
sink()

cat('\n\tΤΑ ΔΗ ΝΥΝ ΠΑΝΤΑ ΤΕΛΕΙΤΑΙ\n\n')

quit('no')
