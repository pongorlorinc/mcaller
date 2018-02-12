import os
import optparse

import pysam
import sys
import numpy

import scipy
from scipy.stats import poisson
from scipy.stats import fisher_exact

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = optparse.OptionParser()

parser.add_option('-g', '--ref',
    action="store", dest="genome_file",
    help="path to fasta (required)", default="NA")

parser.add_option('-b', '--bed',
    action="store", dest="bed",
    help="bed file (required)", default="NA")

parser.add_option('-i', '--in',
    action="store", dest="bam",
    help="bam files, semicolon \';\' separated (required)\n", default="NA")

parser.add_option('-o', '--out',
    action="store", dest="out_folder",
    help="output folder to store files (required)\n", default="NA")

parser.add_option('-y', '--coord', action='store_true',
    help='treat bed as exact coordiantes (quicker if mutations are known)', dest="coord", default=0)

parser.add_option('-x', '--aln', type=int,
    action="store", dest="min_alignment_score",
    help="minimum alignment score", default=1)

parser.add_option('-s', '--reg', type=int,
    action="store", dest="region_interval",
    help="region size for noise calculation", default=10)

parser.add_option('-n', '--noise', type=int,
    action="store", dest="noise",
    help="noise probability", default=0.05)

parser.add_option('-t', '--ncont', type=int,
    action="store", dest="normal_mut_thresh",
    help="normal sample max contamination", default=1)

parser.add_option('-w', '--mincov', type=int,
    action="store", dest="min_cov",
    help="minimum coverage during analysis", default=5)

parser.add_option('-e', '--edit', type=int,
    action="store", dest="max_read_edit_distance",
    help="max edit distance", default=4)

parser.add_option('-c', '--covfilt', type=int,
    action="store", dest="filter_coverage",
    help="filter mutations where coverage was not reached by samples", default=100)

parser.add_option('-p', '--pval', type=int,
    action="store", dest="background_p_thresh",
    help="p-value cutoff", default=0.05)

parser.add_option('-q', '--qfreq', type=int,
    action="store", dest="freq_filter_high_cov",
    help="mutation frequency filter for high coverage regions", default=2)

parser.add_option('-f', '--freq', type=int,
    action="store", dest="freq_thresh",
    help="mutation frequency filter for non-high coverage regions", default=5)

parser.add_option('-u', '--ncov', type=int,
    action="store", dest="min_norm_cov",
    help="minimum reads in normal sample", default=5)

parser.add_option('-m', '--mut', type=int,
    action="store", dest="min_somatic_mut",
    help="minimum mutant reads in cancer sample", default=4)

parser.add_option('-a', '--indel', type=int,
    action="store", dest="minimum_indel_mut",
    help="minimum mutant reads supporting indel in cancer sample", default=4)

parser.add_option('-d', '--maxprox', type=int,
    action="store", dest="prox_indel_mut_count",
    help="maximum proximal reads with indel before mutation filtering", default=2)

parser.add_option('-j', '--indeldist', type=int,
    action="store", dest="prox_indel_dist",
    help="distance when looking for proximal indels", default=5)

parser.add_option('-k', '--proxfile', action='store_true',
    help='(boolean) filter for proximal indels (0:no / 1:yes)', dest="proximal_indel_filter", default=1)

parser.add_option('-l', '--indnoise', type=int,
    action="store", dest="indel_noise_freq",
    help="default noise when calculating background", default=0.01)

options, args = parser.parse_args()

if(options.genome_file == "NA" or options.bed == "NA" or options.bam == "NA" or options.out_folder == "NA"):
	parser.print_help()
	sys.exit()




out_folder = options.out_folder

if not os.path.exists(out_folder):
	os.makedirs(out_folder)

if options.genome_file.endswith(".fa"):
	genome_file = options.genome_file

bamlist = list()
bamlist = options.bam.rstrip().split(',')

num_of_samples = len(bamlist)
print ("\n---------------------------------")
print ("mcaller: joint genotyping program")
print ("---------------------------------")
print ("Main parameters:\n")
print ("\tBAM files:\t%s" % bamlist)
print ("\tGenome:\t\t%s" % genome_file)
print ("\tOutfold:\t\"%s\"" % out_folder)
print ("\tBEDfile:\t%s" % options.bed)
print ("\n\nAdditional parameters:\n")
print ("\tTreat bed as coordinate: %d\n" %options.coord)
print ("\t(boolean) filter for proximal indels (0:no / 1:yes): %d" %options.proximal_indel_filter)
print ("\tDistance when looking for proximal indels: %d" %options.prox_indel_dist)
print ("\tMaximum proximal reads with indel before mutation filtering: %d\n" %options.prox_indel_mut_count)
print ("\tMinimum alignment score: %d" %options.min_alignment_score)
print ("\tMax edit distance in read: %d" %options.max_read_edit_distance)
print ("\tRegion size for noise calculation: %d" %options.region_interval)
print ("\tSequencing error probability: %d" %options.noise)
print ("\tNormal sample max contamination: %d" %options.normal_mut_thresh)
print ("\tMinimum coverage during analysis: %d\n" %options.min_cov)
print ("\tMutation p-value cutoff: %d" %options.background_p_thresh)
print ("\tMutation frequency filter for high coverage regions: %d" %options.freq_filter_high_cov)
print ("\tMutation frequency filter for low coverage regions: %d" %options.freq_thresh)
print ("\tMinimum mutant reads in (any) cancer sample: %d" %options.min_somatic_mut)
print ("\tMinimum mutant reads supporting indel in cancer sample: %d" %options.minimum_indel_mut)
print ("\tMinimum reads in normal sample to call germline mutation: %d" %options.min_norm_cov)
print ("\tFilter mutations where coverage was not reached by any samples: %d" %options.filter_coverage)
print ("---------------------------------\n")
input_coords = options.coord

bedlogfile = out_folder + "/entries_finished.bed"
bedloghandle = open(bedlogfile, "w+")

mutationfile = out_folder + "/mutations.txt"
mutationhandle = open(mutationfile, "w+")

mutationhandle.write("chr\tpos\tref\talt\tstatus\taccept\t%s\n" % ("\t".join(os.path.basename(bamlist[y]) for y in range(0, len(bamlist)))))

print ("Importing genome:")
records = SeqIO.to_dict(SeqIO.parse(open(genome_file), 'fasta'))
print ("\tfinished importing genome")

###################### VARIABLES ######################
base = dict()
rb = dict()

base["A"] = 0
base["C"] = 1
base["G"] = 2
base["T"] = 3
rb[0] = "A"
rb[1] = "C"
rb[2] = "G"
rb[3] = "T"

min_alignment_score = options.min_alignment_score
region_interval = options.region_interval

noise = options.noise
background_p_thresh = options.background_p_thresh

q_max = 42
q_m_height = 4 * 42

min_cov = options.min_cov

max_read_edit_distance = options.max_read_edit_distance

filter_coverage = options.filter_coverage
normal_mut_thresh = options.normal_mut_thresh
freq_filter_high_cov = options.freq_filter_high_cov
freq_thresh = options.freq_thresh
min_somatic_mut = options.min_somatic_mut

min_norm_cov = options.min_norm_cov

indel_noise_freq = options.indel_noise_freq

proximal_indel_filter = options.proximal_indel_filter
prox_indel_dist = options.prox_indel_dist
prox_indel_mut_count = options.prox_indel_mut_count
minimum_indel_mut = options.minimum_indel_mut

###################### VARIABLES ######################



def get_edit_distance(read, in_chr, start, g_m_start, m_end, mismatch):
	edit_distance = 0

	genome_pos = read.reference_start
	read_pos = 0

	for (cigarType,cigarLength) in read.cigartuples:
		if cigarType == 0: #MATCH bases in read (M cigar code)
			genome_end = genome_pos + cigarLength
			read_end = read_pos + cigarLength
					
			ref_seq = records[in_chr].seq[genome_pos:genome_end].upper()
			read_seq = read.query_sequence[read_pos:read_end]
			frag_len = len(ref_seq)

			for i in range(0, frag_len):
				t_g_pos = genome_pos + i - g_m_start

				if read_seq[i] in base:
					if t_g_pos < m_end:
						if ref_seq[i] != read_seq[i]:
							mismatch[t_g_pos] += 1
							edit_distance += 1

			genome_pos = genome_end
			read_pos = read_end

		elif cigarType == 1: #INS in read (I cigar code)
			t_g_pos = genome_pos - g_m_start
			read_end = read_pos + cigarLength
			read_pos = read_end
			edit_distance += 1

			for x in range(-2, 2):
				temp_pos = t_g_pos + x

				if temp_pos < m_end:
					mismatch[temp_pos] += 1

		elif cigarType == 2: #DELETION in read (D cigar code)
			t_g_pos = genome_pos - g_m_start
			genome_end = genome_pos + cigarLength
			genome_pos = genome_end
			edit_distance += 1

			for x in range(-2, 2):
				temp_pos = t_g_pos + x

				if temp_pos < m_end:
					mismatch[temp_pos] += 1

		elif cigarType == 3: #SPLICED region in read representing numer of skipped bases on reference (N cigar code)
			genome_pos += cigarLength #jump to next coordinate by skipping bases

		elif cigarType == 4: #SOFTCLIP
			read_pos += cigarLength

		else:
			pass

	return edit_distance, mismatch	

#Input: pysam read type, which contains read information and ref_name (chromosome name where read aligned)
#Output: Boolean, True is read accepted, False if rejected 
def accept_read(in_read):
	if in_read.mapping_quality < min_alignment_score:
		return 0

	if in_read.is_duplicate:
		return 0

	if in_read.reference_id < 0:
		return 0

	if in_read.next_reference_id < 0:
		return 0

	if in_read.reference_name != in_read.next_reference_name:
		return 0

	if in_read.cigartuples is None:
		return 0

	if in_read.cigartuples[0][0] == 5 or in_read.cigartuples[-1][0] == 5:
		return 0

	return 1

def import_reads(in_chr, start, end, samfile, g_m_start, m_end, mismatch, reads):
	for read in samfile.fetch(in_chr, start, end):
		if accept_read(read) == 1:
			edit_distance, mismatch = get_edit_distance(read, in_chr, start, g_m_start, m_end, mismatch)

			if edit_distance < max_read_edit_distance:
				reads.append(read)

	return reads, mismatch

def import_data(reads, coverage, base_array, r_q_matrix, f_q_matrix, m_array, in_chr, g_m_start, m_end, active_pos, indels, s_num, indarray, mismatch):
	for read in reads:
		genome_pos = read.reference_start
		read_pos = 0
		
		for (cigarType,cigarLength) in read.cigartuples:
			if cigarType == 0: #MATCH bases in read (M cigar code)
				genome_end = genome_pos + cigarLength
				read_end = read_pos + cigarLength
					
				ref_seq = records[in_chr].seq[genome_pos:genome_end].upper()
				read_seq = read.query_sequence[read_pos:read_end]
				read_qual = read.query_qualities[read_pos:read_end]

				frag_len = len(ref_seq)
				t_g_pos = genome_pos - g_m_start - 1

				for i in range(0, frag_len):
					t_g_pos += 1

					if t_g_pos < m_end and read_seq[i] in base and mismatch[t_g_pos] > 1:
						
						t_index = q_max*base[read_seq[i]] + read_qual[i]

						if t_index > -1 and t_index < q_m_height:
							if read.is_reverse:
								r_q_matrix[t_index][t_g_pos] += 1

							else:
								f_q_matrix[t_index][t_g_pos] += 1

							if ref_seq[i] != read_seq[i]:
								m_array[t_g_pos] += 1
								active_pos[t_g_pos] = 1

							coverage[t_g_pos] += 1

						base_array[base[read_seq[i]]][t_g_pos] += 1


				genome_pos = genome_end
				read_pos = read_end

			elif cigarType == 1: #INS in read (I cigar code)
				read_end = read_pos + cigarLength
				read_seq = read.query_sequence[read_pos-1:read_end]
				ref_seq = records[in_chr].seq[genome_pos-1]
				t_g_pos = genome_pos - g_m_start - 1
				indel_key = ref_seq + "\t" + read_seq

				indel_strand = 0

				if t_g_pos < m_end:
					if read.is_reverse:
						indel_strand = 1

					if t_g_pos not in indels:
						indels[t_g_pos] = {}

					if indel_key not in indels[t_g_pos]:
						indels[t_g_pos][indel_key] = ([[s_num, min(read.query_qualities[read_pos:read_end]), indel_strand, t_g_pos + 1]])

					else:
						indels[t_g_pos][indel_key].append([s_num, min(read.query_qualities[read_pos:read_end]), indel_strand, t_g_pos + 1])

					indarray[t_g_pos] += 1
					active_pos[t_g_pos] = 1
				read_pos = read_end

			elif cigarType == 2: #DELETION in read (D cigar code)
				genome_end = genome_pos + cigarLength
				ref_seq = "".join(records[in_chr].seq[y] for y in range(genome_pos-1, genome_end))
				read_seq = read.query_sequence[read_pos-1:read_pos]
				t_g_pos = genome_pos - g_m_start

				if t_g_pos < m_end:
					indel_key = ref_seq + "\t" + read_seq
					indel_strand = 0

					if read.is_reverse:
						indel_strand = 1

					if t_g_pos not in indels:
						indels[t_g_pos] = {}

					if indel_key not in indels[t_g_pos]:
						indels[t_g_pos][indel_key] = ([[s_num, min(read.query_qualities[read_pos:read_pos+1]), indel_strand, t_g_pos]])
						#print (indels[t_g_pos])

					else:
						indels[t_g_pos][indel_key].append([s_num, min(read.query_qualities[read_pos:read_pos+1]), indel_strand, t_g_pos])
						#print (indels[t_g_pos])


					indarray[t_g_pos] += 1
					active_pos[t_g_pos] = 1

				genome_pos = genome_end

			elif cigarType == 3: #SPLICED region in read representing numer of skipped bases on reference (N cigar code)
				genome_pos += cigarLength #jump to next coordinate by skipping bases

			elif cigarType == 4: #SOFTCLIP
				read_pos += cigarLength

			else:
				pass

	return coverage, base_array, r_q_matrix, f_q_matrix, m_array, active_pos, indels, indarray

def check_coverage_thresh(cov, pos, num_of_samples):
	for i in range(0, num_of_samples):
		if cov[i][pos] < min_cov:
			return 0

	return 1

def active_region(cov, mut, pos, num_of_samples, verbose):
	for i in range(0, num_of_samples):
		if cov[i][pos] > 0:
			background_noise = cov[i][pos] * noise
			background_test = poisson.cdf(background_noise, mut[i][pos])

			if background_test < background_p_thresh:
				return background_test
			

	return 1

def calculate_region_pvalues(pos, num_of_samples, b_a, ref_segment, m_end, mut, cov, st):
	background_probs = list()

	region_start = pos - region_interval
	region_end = pos + region_interval + 1
	min_p = 1

	if region_start < 0:
		region_end += (-1) * region_start
		region_start = 0

	if region_end > m_end:
		region_start -= region_end - m_end
		region_end = m_end

	for j in range(0, num_of_samples):
		region_noise = 0
		noise_sum = 0
		num_of_bases = 0
		region_mean_cov = 0
		region_mean_mut = 0

		for i in range(region_start, region_end):
			if mut[j][i] > 0 and i != pos and mut[j][i] / cov[j][i] < 0.05:
				region_noise += mut[j][i]
				noise_sum += cov[j][i]
				num_of_bases += 1

		if region_noise > 0 and num_of_bases > 0:
			region_mean_cov = noise_sum / num_of_bases
			region_mean_mut = region_noise / num_of_bases
			region_noise = (region_mean_mut/region_mean_cov) * cov[j][pos] + noise * cov[j][pos]

		else:
			region_noise = noise * cov[j][pos]

		background_probs.append([poisson.cdf(region_noise, b_a[j][0][pos]), poisson.cdf(region_noise, b_a[j][1][pos]), poisson.cdf(region_noise, b_a[j][2][pos]), poisson.cdf(region_noise, b_a[j][3][pos])])
		background_probs[-1][base[ref_segment[pos]]] = 1.0

		for i in range(0, 4):
			if background_probs[-1][i] < min_p:
				min_p = background_probs[-1][i]
 
	return background_probs, min_p

def calculate_strand_and_qual_probs(pos, num_of_samples, f_q_m, r_q_m, ref_segment, background_probs):
	strand_probs = list()
	qual_probs = list()

	for j in range(0, num_of_samples):
		forward_sum, reverse_sum, f_q_sum, r_q_sum = 1, 1, 1, 1
		strand_probs.append([1, 1, 1, 1, 1])
		qual_probs.append([1, 1, 1, 1, 1])
		qual_sums = 0

		for k in range(0, q_m_height):
			forward_sum += f_q_m[j][k][pos]
			reverse_sum += r_q_m[j][k][pos]

			if f_q_m[j][k][pos] > 0:
				qual_sums += f_q_m[j][k][pos] * (k - k/q_max * q_max)

			if r_q_m[j][k][pos] > 0:
				qual_sums += r_q_m[j][k][pos] * (k - k/q_max * q_max)

			
		strand_probs[-1][0] = fisher_exact([[(forward_sum+reverse_sum)/2, (forward_sum+reverse_sum)/2], [forward_sum, reverse_sum]])[1]
		qual_sums = qual_sums / (forward_sum + reverse_sum)

		for k in range(0, 4):
			if background_probs[j][k] < 0.05:
				local_forward_sum = 0
				local_reverse_sum = 0

				for l in range(k*q_max, (k+1)*q_max):
					local_forward_sum += f_q_m[j][l][pos]
					local_reverse_sum += r_q_m[j][l][pos]

					if f_q_m[j][l][pos] > 0:
						qual_probs[-1][k] += f_q_m[j][l][pos] * (l - l/q_max * q_max)

					if r_q_m[j][l][pos] > 0:
						qual_probs[-1][k] += r_q_m[j][l][pos] * (l - l/q_max * q_max)

				if local_forward_sum > 0 or local_reverse_sum > 0:
					strand_probs[-1][k+1] = fisher_exact([[local_forward_sum, local_reverse_sum], [forward_sum, reverse_sum]])[1]
					qual_probs[-1][k] = poisson.cdf(qual_probs[-1][k]/(local_forward_sum + local_reverse_sum), qual_sums)

				else:
					qual_probs[-1][k] = 1

	return strand_probs, qual_probs

def get_region(bed_input):
	in_chr = bed_input[0]
	start = int(bed_input[1])
	end = int(bed_input[2])
	end = end + 1

	local_muts = list()

	if in_chr not in records:
		return

	g_m_start = -1
	g_m_end = -1

	m_start = 0
	m_end = 0

	
	g_m_start = start - 1000
	g_m_end = end

	ref_segment = records[in_chr].seq[g_m_start:g_m_end+100].upper()

	m_end = end - g_m_start + 100
	
	#initialize matrices
	f_q_m = list()
	r_q_m = list()
	b_a = list()
	cov = list()
	mut = list()
	reads = [[] for x in xrange(num_of_samples)]
	mismatch = [0]*m_end
	indarray = list()
	active_pos = dict()
	indels = dict()

	for i in range(0, num_of_samples):
		f_q_m.append(numpy.zeros((q_m_height, m_end), dtype=numpy.int))
		r_q_m.append(numpy.zeros((q_m_height, m_end), dtype=numpy.int))

		b_a.append(numpy.zeros((4, m_end+100), dtype=numpy.int))
		cov.append(numpy.zeros((m_end), dtype=numpy.int))
		indarray.append(numpy.zeros((m_end), dtype=numpy.int))
		mut.append(numpy.zeros((m_end), dtype=numpy.int))

	#matrix initialize over
	for i in range(0, num_of_samples):
		r_bam = pysam.AlignmentFile(bamlist[i], "rb")
		reads[i], mismatch = import_reads(in_chr, start, end, r_bam, g_m_start, m_end, mismatch, reads[i])
		r_bam.close()

	for i in range(0, num_of_samples):
		cov[i], b_a[i], r_q_m[i], f_q_m[i], mut[i], active_pos, indels, indarray[i] = import_data(reads[i], cov[i], b_a[i], r_q_m[i], f_q_m[i], mut[i], in_chr, g_m_start, m_end, active_pos, indels, i, indarray[i], mismatch)

	if input_coords == 0:
		active_coords = active_pos.keys()
		active_coords.sort()

	else:
		active_coords = list()
		active_coords = [start-g_m_start-1]

	for i in active_coords:
		ref_base = ref_segment[i]

		if i in indels.keys():
			for key in indels[i].keys():
				ind_cov = [0] * num_of_samples
				ind_freq = [0] * num_of_samples
				ind_forw = [0] * num_of_samples
				ind_rev = [0] * num_of_samples
				ind_quals = [0] * num_of_samples
				qual_sums = [0] * num_of_samples
				forward_sums = [0] * num_of_samples
				reverse_sums = [0] * num_of_samples

				indel_probs = [1] * num_of_samples

				indel_norm_qual = [0] * num_of_samples
				norm_qual = [0] * num_of_samples

				qual_bias = [1] * num_of_samples
				strand_bias = [1] * num_of_samples

				indel_result = list()

				for els in indels[i][key]:
					if els[0] <= num_of_samples:
						ind_cov[els[0]] += 1

						if els[2] == 0:
							ind_forw[els[0]] += 1

						else:
							ind_rev[els[0]] += 1

						ind_quals[els[0]] += els[1]

				for j in range(0, num_of_samples):
					if ind_cov[j] > 0 and cov[j][i] > 0:
						ind_freq[j] = ind_cov[j] * 100 / cov[j][i]


				for j in range(0, num_of_samples):
					for k in range(0, q_m_height):
						forward_sums[j] += f_q_m[j][k][i]
						reverse_sums[j] += r_q_m[j][k][i]

						if f_q_m[j][k][i] > 0:
							qual_sums[j] += f_q_m[j][k][i] * (k - k/q_max * q_max)

						if r_q_m[j][k][i] > 0:
							qual_sums[j] += r_q_m[j][k][i] * (k - k/q_max * q_max)

				for j in range(0, num_of_samples):
					if cov[j][i] > 0:
						norm_qual[j] = qual_sums[j] / cov[j][i]

					if ind_cov[j] > 1:
						indel_probs[j] = poisson.cdf(int(cov[j][i] * indel_noise_freq), ind_cov[j])
						indel_norm_qual[j] = ind_quals[j] / ind_cov[j]

				for j in range(0, num_of_samples):
						if indel_norm_qual[j] > 0 and norm_qual[j] > 0:
							qual_bias[j] = poisson.cdf(indel_norm_qual[j], norm_qual[j])

						if ind_cov > 0:
							strand_bias[j] = fisher_exact([[ind_forw[j], ind_rev[j]],[forward_sums[j], reverse_sums[j]]])[1]

				for j in range(0, num_of_samples):
					indel_result.append(str(cov[j][i]) + ";" + str(ind_cov[j]) + "," + str(ind_forw[j]) +";" + str(ind_rev[j]) + "," + str(forward_sums[j]) + ";" + str(reverse_sums[j]) + "," + str("%.2g" % indel_probs[j]) + "," + str("%.2g" % qual_bias[j]) + "," + str("%.2g" % strand_bias[j]))

				germline = 0
				somatic = 0
				accept = 0

				if indel_probs[0] < 0.05:
					if ind_freq[0] > 5 and ind_cov[0] > 4:
						accept = 1

					mutationhandle.write("%s\t%d\t%s\tgermline\t%d\t%s\n" % (in_chr, els[3] + g_m_start, key, accept, "\t".join(indel_result)))

				else:
					if min(indel_probs[1:num_of_samples]) < 0.05:
						s_accept = 1

						if max(qual_bias[1:num_of_samples]) < 0.05:
							s_accpet = 0

						if cov[0][i] < filter_coverage:
							if ind_cov[0] > 0:
								s_accept = 0

						else:
							if ind_cov[0] > normal_mut_thresh:
								s_accept = 0

							if max(ind_freq[1:num_of_samples]) < freq_filter_high_cov:
								s_accept = 0

						if max(ind_cov[1:num_of_samples]) < minimum_indel_mut:
							s_accept = 0

						if s_accept == 1:
							accept = 1

						if cov[0][i] < min_norm_cov:
							accept = 0

					mutationhandle.write("%s\t%d\t%s\tsomatic\t%d\t%s\n" % (in_chr, els[3] + g_m_start, key, accept, "\t".join(indel_result)))

		if ref_base in base:
			background_probs, min_p = calculate_region_pvalues(i, num_of_samples, b_a, ref_segment, m_end, mut, cov, g_m_start)
			background_p = active_region(cov, mut, i, num_of_samples, 0)
			
			if background_p < 0.05 or min_p < 0.05:
				strand_probs, qual_probs = calculate_strand_and_qual_probs(i, num_of_samples, f_q_m, r_q_m, ref_segment, background_probs)

				multi_normal = 0
				normal_mut = -1
				multi_allelic = 0
				somatic = -1

				somatic = list()
				germline = list()

				sample_results = list()

				for j in range(0, num_of_samples):
					sample_results.append( str(cov[j][i])+";"+",".join(str(b_a[j][y][i]) for y in range(0, 4)) + ";" + ','.join(str("%.2g" % x) for x in background_probs[j]) + ";" + ':'.join(str("%.2g" % x) for x in strand_probs[j]) + ";" + ','.join(str("%.2g" % x) for x in qual_probs[j]))
				sample_results = '\t'.join(str(x) for x in sample_results)

				for j in range(0, 4):
					if background_probs[0][j] < 0.05 and j != base[ref_base]:
						germline.append(j)


				for j in range(0, 4):
					local_mut = 0

					for k in range(1, num_of_samples):
						if background_probs[k][j] < 0.05 and j != base[ref_base]:
							local_mut = 1

					if local_mut == 1 and j != base[ref_base]:
						in_germline = 0

						for l in range(0,len(germline)):
							if germline[l] == j:
								in_germline = 1

							if in_germline == 0:
								somatic.append(j)

				if len(germline) > 0:
					mut_type = "germline"

					if len(germline) > 1:
						mut_type = "germline_multi"

					for j in range(0, len(germline)):
						if b_a[0][germline[j]][i] > 4 and b_a[0][germline[j]][i] * 100 / cov[0][i] > 5:
							mutationhandle.write("%s\t%d\t%c\t%c\t%s\t1\t%s\n" % (in_chr, i + g_m_start + 1, ref_base, rb[germline[j]], mut_type, sample_results))

						else:
							mutationhandle.write("%s\t%d\t%c\t%c\t%s\t0\t%s\n" % (in_chr, i + g_m_start + 1, ref_base, rb[germline[j]], mut_type, sample_results))

				if len(somatic) > 0:
					mut_type = "somatic"

					if len(somatic) > 1:
						mut_type = "somatic_multi"

					for j in range(0, len(somatic)):
						accept = 0

						for k in range(1, num_of_samples):
							s_accept = 1

							if cov[k][i] > 0:
								if qual_probs[k][somatic[j]] < 0.05:
									s_accept = 0

								if cov[k][i] < filter_coverage:
									if b_a[0][somatic[j]][i] > 0:
										s_accept = 0




								else:
									if b_a[0][somatic[j]][i] > normal_mut_thresh:
										s_accept = 0

									if b_a[k][somatic[j]][i]*100/cov[k][i] < freq_filter_high_cov:
										s_accept = 0


								if b_a[k][somatic[j]][i] < min_somatic_mut:
									if b_a[k][somatic[j]][i]*100/cov[k][i] < freq_thresh:
										s_accept = 0

								if proximal_indel_filter == 1:
									if max(indarray[k][i-prox_indel_dist:i+prox_indel_dist]) >= prox_indel_mut_count:
										s_accept = 0

							else:
								s_accept = 0

							if s_accept == 1:
								accept = 1
						if min_norm_cov > cov[0][i]:
							accept = 0

						mutationhandle.write("%s\t%d\t%c\t%c\tsomatic\t%d\t%s\n" % (in_chr, i + g_m_start + 1, ref_base, rb[somatic[j]], accept, sample_results))

	tbed = "\t".join(bed_input)
	print (tbed)
	bedloghandle.write("%s\n" % (tbed))

	
bedfile = open(options.bed, "r")

genes = [line.rstrip().split('\t') for line in bedfile]

for k in genes:
	get_region(k)
mutationhandle.close()
bedloghandle.close()

