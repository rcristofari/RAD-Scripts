import argparse, time, os, random
from Bio import SeqIO, AlignIO
from collections import Counter

parser = argparse.ArgumentParser(description="Create an input map for BEAST EPSB from Stacks output")
parser.add_argument('--fasta', help='Fasta output from Stacks')
parser.add_argument('--popmap', help='Population map (tab-separated')
parser.add_argument('--n_ind', help="coma-separated list of the minimum number of individuals per populations [all]")
parser.add_argument('--exact', help="Number of individuals per population are exact (downsample the locus if necessary)", action='store_true')
parser.add_argument('--jointly', help="Only output a locus if it has the required sample sizes jointly in all populations", action='store_true')
parser.add_argument('--haploid', help='Randomly pick one allele for each individual e.g. to remove low-coverage-related homozygous bias', action='store_true')
parser.add_argument('--skip_fixed', help="Do not output non-polymorphic loci", action='store_true')
parser.add_argument('--out', help="output directory")
args = parser.parse_args()
start_time = time.time()
print("---------------------------------------------------------------------------------")

#----------------------------------------------------------------------------------------------------------------------#
# Parse and handle command line arguments:
if args.fasta:
    fasta_path = args.fasta
    if not os.path.isfile(fasta_path):
        print(f"Fasta file {fasta_path} could not be found")
        quit()
else:
    print("You must specify an input fasta file - exiting.")
    quit()

if args.popmap:
    popmap_path = args.popmap
    if not os.path.isfile(popmap_path):
        print(f"Population map {popmap_path} could not be found")
        quit()
else:
    print("You must specify an input population map file - exiting.")
    quit()

if args.out:
    out_path = str(args.out).rstrip("/") + "/"
    if not os.path.exists(out_path):
        print(f"Output directory {out_path} could not be found")
        quit()
else:
    print("You must specify an output directory - exiting.")
    quit()

do_haploid = args.haploid
strploid = "haploid" if do_haploid else "diploid"

do_skip = args.skip_fixed

if args.n_ind:
    try:
        n_ind = [int(x) for x in args.n_ind.split(",")]
    except ValueError:
        print(f"Invalid number of indidivuals {args.n_ind}")
        quit()
    do_exact = args.exact
    do_jointly = args.jointly
else:
    n_ind = None

#----------------------------------------------------------------------------------------------------------------------#
# Open and parse the population map
samples, pops = [], []
with open(popmap_path, 'r') as ifile:
    for line in ifile:
        if not line.startswith("#"):
            sample, pop = line.strip("\n").split("\t")
            samples.append(sample)
            pops.append(pop)
pop_dict = dict(zip(samples, pops)) # The sample-pop correspondence
n_pop = len(set(pops))
pop_idx = dict(zip(set(pops), [x for x in range(n_pop)])) # The pop-integer ID correspondence
# Make a dict with the number of requested individuals per pop:
if n_ind:
    n_ind_pop = dict(zip([x for x in set(pops)], n_ind))
else:
    n_ind_pop = dict(zip([x for x in set(pops)], [0]*n_pop))

#----------------------------------------------------------------------------------------------------------------------#
# A function to count the number of SNPs at a locus
def n_snps(locus):
    # Check if a position is monomorphic
    def is_monomorphic(x):
        x = [y for y in x if y != "N"]
        return all(y == x[0] for y in x)
    # Apply that to all positions:
    SNPS = 0
    for pos in range(locus.get_alignment_length()):
        SNPS += not is_monomorphic(locus[:, pos])
    return SNPS

#----------------------------------------------------------------------------------------------------------------------#
# A function to process a locus and write out the result
def process_locus(locus):
    # Split the locus in a list of alignment objects, by population:
    locus_pop = [AlignIO.MultipleSeqAlignment(None) for i in range(n_pop)]
    for record in locus:
        sample = record.description.split(" ")[-1][1:-1]
        idx = pop_idx[pop_dict[sample]]
        locus_pop[idx].append(record)

    # Iterate through each population and process it
    for i in range(n_pop):
        # Duplicate the haploid sequences:
        samples = []
        for record in locus_pop[i]:
            samples.append(record.description.split(" ")[-1][1:-1])
        haploid = [i for i, x in enumerate([Counter(samples)[s] for s in samples]) if x == 1]
        for h in haploid:
            locus_pop[i].append(locus_pop[i][h])
        locus_pop[i].sort(reverse=False)
        # Check if the locus passes the sample size controls, and downsample if necessary:
        n_samples = int(len(locus_pop[i])/2)

        if (n_samples >= n_ind_pop[pop_dict[sample]] and n_samples != 0) and not (do_exact is True and n_ind_pop[pop_dict[sample]] == 0):
            # Subsample randomly some individuals if we want exact sample sizes:
            if do_exact:
                random_idx = sorted(random.sample([x for x in range(n_samples)], n_ind_pop[pop_dict[sample]]))
                subsampled_locus = AlignIO.MultipleSeqAlignment(None)
                for x in random_idx:
                    subsampled_locus.append(locus_pop[i][x*2])
                    subsampled_locus.append(locus_pop[i][x*2+1])
                locus_pop[i] = subsampled_locus
                n_samples = n_ind_pop[pop_dict[sample]]

            # Subsample randomly one of the two alleles if we want haploid:
            if do_haploid:
                random_idx = [x*2 + random.randint(0, 1) for x in range(n_samples)]
                haploid_locus = AlignIO.MultipleSeqAlignment(None)
                for x in random_idx:
                    haploid_locus.append(locus_pop[i][x])
                locus_pop[i] = haploid_locus

            # Write out the data to a NEXUS file separately
            if not do_jointly:
                snps = n_snps(locus_pop[i])
                if not do_skip or (do_skip and snps > 0):
                    this_out_path = out_path + f"Locus_{locus_pop[i][0].description.split('_')[1]}_{list(pop_idx.keys())[i].upper()}_{n_samples}_{strploid}_samples_{snps}_snp.nex"
                    AlignIO.write(locus_pop[i], this_out_path, "nexus")

        # Delete altogether if it does not match sample size requirements:
        else:
            locus_pop[i] = None

    # If we call jointly, go back on the whole locus before writing out:
    if do_jointly:
        available_sizes = [len(locus_pop[i])/2 if locus_pop[i] else 0 for i in range(n_pop)]
        if all(available_sizes[i] >= n_ind[i] for i in range(n_pop)):
            for i in range(n_pop):
                snps = n_snps(locus_pop[i])
                if not do_skip or (do_skip and snps > 0):
                    this_out_path = out_path + f"Locus_{locus_pop[i][0].description.split('_')[1]}_{list(pop_idx.keys())[i].upper()}_{n_samples}_{strploid}_samples_{snps}_snp.nex"
                    AlignIO.write(locus_pop[i], this_out_path, "nexus")

#----------------------------------------------------------------------------------------------------------------------#
# MAIN LOOP
# Count lines in the file (for reporting purposes)
n_lines = sum(1 for i in open(fasta_path, 'rb'))/2
# Set the main vars
sample, locus_id, allele = None, None, None
prev_sample, prev_locus_id, prev_allele = None, None, None
this_locus = []

#----------------------------------------------------------------------------------------------------------------------#
start_time = time.time()
with open(fasta_path) as ifile:
    for K, record in enumerate(SeqIO.parse(ifile, 'fasta')):
        locus_id = record.id.split("_")[1]
        record.annotations["molecule_type"] = "DNA"
        if not this_locus or locus_id == prev_locus_id:
            this_locus.append(record)
        else:
            # Process the locus, and initiate the next one:
            process_locus(this_locus)
            this_locus = [record]
        prev_locus_id = locus_id

        # Handle running display:
        if K % 1000 == 0:
            elapsed = int(round(time.time() - start_time))
            pct = round(100*(K/n_lines), 1)
            remaining = "---"
            if pct > 0:
                remaining = round(elapsed * 100 / pct)
            print(f"{time.asctime()} | Processed {K} sequences in {elapsed} seconds | {pct}% done | {remaining} seconds remaining\r", end="")

    # Process the last locus
    process_locus(this_locus)
#----------------------------------------------------------------------------------------------------------------------#

print(f"\nFinished in {time.time() - start_time} seconds")
