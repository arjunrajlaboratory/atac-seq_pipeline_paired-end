import sys
import logging

from subprocess import Popen, PIPE
import os
import gzip
import random

# note: the call to this script should be preprended with a 'set_atac_pmacs_env' command that sets the paths for the appropriate versions of each tool called in this pipeline. 

# TODO (optional). Rewrite into "functional" style.
# TODO (optional). Process arguments with argparse

sample_name       = sys.argv[1]
input_data_dir    = sys.argv[2]
output_dir        = sys.argv[3]
n_bowtie2_threads = sys.argv[4]
read1_path     = '{0}/{1}'.format(input_data_dir, sample_name + '_R1.fastq.gz')
read2_path     = '{0}/{1}'.format(input_data_dir, sample_name + '_R2.fastq.gz')

#software references
refs_dir = '/home/esanford/refs'
bowtie2_index     = refs_dir + '/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
gencode_tss_areas = refs_dir + '/GENCODEv24_TSS_hg38_radius500.bed'
refseq_tss_areas  = refs_dir + '/UCSC_RefSeq_TSS_hg38_radius500.bed'
chrom_sizes_file  = '/home/esanford/software/bedGraphToBigWig/hg38.chrom.sizes'
hg38_blacklist    = '/home/esanford/refs/hg38.blacklist.sorted.bed'

#setup output filepaths
output_dir_sample  = output_dir + '/' + sample_name
if not os.path.exists(output_dir_sample):
	os.makedirs(output_dir_sample)

bowtie2_stderr_logfile         = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'log_stderr_bowtie2.txt')
bowtie2_aligned_reads          = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'bowtie2_aligned_reads.bam')

tmp_early_filt_output          = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'tmp.lightFilt.bam')
early_filt_fixmate_output      = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'lightFilt.fixmate.bam')
filtered_bam_file              = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'filt.bam')
filtered_bam_file_namesorted   = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'filt.nameSorted.bam')
picard_metrics_file            = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'picardMetrics.txt')
picard_output                  = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'picardMarkedDups.bam')
final_bam_file                 = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.bam')
final_bam_stats                = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.stats.txt')
final_bam_file_namesorted      = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.nameSorted.bam')

early_tagalign_file            = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'early.tagAlign.gz')
final_tagalign_file            = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.tagAlign.gz')
final_Tn5shifted_tagAlign_file = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.tagAlign_tn5_shifted.gz')
final_bedpe_file               = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.bedpe.gz')
macs2_unsorted_tagAlign_file   = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'macs2_formatted.unsorted.tagAlign')
macs2_tagAlign_file            = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'Tn5_insertion_points.tagAlign.gz')


subsampled_tagalign_file       = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'subsampled.tagAlign.gz')

insert_size_histogram_plot     = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'insert_size_histogram.pdf')
insert_size_histogram_data     = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'insert_size_data.txt')
mitochondrial_read_report      = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'mitochondrial_read_report.txt')
refseq_tss_report              = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'refseq_tss_report.txt')
gencode_tss_report             = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'gencode_tss_report.txt')
pbc_qc_file                    = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.pbc.qc.txt')

#files to be deleted after pipeline runs
intermediate_file_list = [tmp_early_filt_output, early_filt_fixmate_output, early_tagalign_file, filtered_bam_file_namesorted, filtered_bam_file, 
						  final_Tn5shifted_tagAlign_file, picard_output, final_bam_file_namesorted, macs2_unsorted_tagAlign_file, macs2_tagAlign_file]

#constants
bt2_max_frag_len    = 5000
n_subsampled_reads  = 10000000
# note: in encode atac specification, for "ENCODE3" the window size is 73 and the shift size is -37
macs2_smoothing_window_sizes = [ 73, 150,  200]
macs2_shift_sizes            = [-37, -75, -100]   #shift sizes should be int, -1/2 * corresponding smoothing window sizes
macs2_cap_num_summits        = 300000
macs2_summit_window_radius   = 50

def main():
	logger=logging.getLogger()

	##################################
	# Step 0: FastQC for each read file
	##################################

	logger.info('Running FastQC for sample: {0}'.format(sample_name))
	os.system(' '.join(['fastqc', read1_path, '--outdir', output_dir_sample]))
	os.system(' '.join(['fastqc', read2_path, '--outdir', output_dir_sample]))

	##################################
	# Align reads
	##################################

	logger.info('Aligning reads for sample: {0}'.format(sample_name))

	cmd = ' '.join(['bowtie2', '--local', 
		           '--maxins', str(bt2_max_frag_len), '-x', bowtie2_index,
		           '--threads', n_bowtie2_threads, 
		           '-1', read1_path, '-2', read2_path, 
		           '2>', bowtie2_stderr_logfile, 
		      '|', 'samtools', 'view', '-u', '/dev/stdin',	
		      '|', 'samtools', 'sort', '/dev/stdin', bowtie2_aligned_reads[:-4]]) # [:-4] is necessary because samtools adds a bam suffix
	print cmd
	os.system(cmd)	
	logger.info('Finished aligning reads for sample: {0}'.format(sample_name))

	##################################
	# Filter alignments
	##################################
	logger.info('Filtering reads for sample: {0}'.format(sample_name))

	logger.debug('Removing unmapped reads, reads with unmapped mates, or reads not passing vendor QC, then sorting by QNAME: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-F', '524', '-f', '2', '-u', bowtie2_aligned_reads,
						'|', 'samtools', 'sort', '-n', '/dev/stdin', tmp_early_filt_output[:-4]]))
	# samtools view: 
	#               -F 524   (-F is "blacklist" for certain flag bits)
	#					bit 512: "not passing filters, such as platform/vendor quality controls"
	#					bit 8: "next segment in the template unmapped"
	#					bit 4: "segment unmapped"
	#               -f 2     (-f is "whitelist" for certain flag bits) 
	#					bit 2: "each segment properly aligned according to the aligner"

	logger.debug('Adding header, filtering out secondary alignments, and fixing mates: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-h', tmp_early_filt_output, 
				   '|', 'samtools', 'fixmate', '-r', '/dev/stdin', early_filt_fixmate_output]))

	logger.debug('Removing low map-quality alignments, optical duplicates, and sorting by position: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-q', '20', '-F', '1804', '-f', '2', '-u', early_filt_fixmate_output,
						'|', 'samtools', 'sort', '/dev/stdin', filtered_bam_file[:-4]]))
	# samtools view: 
	#               -F 1804  (-F is "blacklist" for certain flag bits)
	#					bit 1024: pcr or optical duplicate 
	#					bit 512: "not passing filters, such as platform/vendor quality controls"
	#					bit 256: secondary alignment. "It is typically used to flag alternative mappings when multiple mappings are presented in a SAM"
	#					bit 8: "next segment in the template unmapped"
	#					bit 4: "segment unmapped"
	#               -f 2     (-f is "whitelist" for certain flag bits) 
	#					bit 2: "each segment properly aligned according to the aligner"

	logger.debug('Marking duplicates with Picard for sample: {0}'.format(sample_name))
	os.system(' '.join(['MarkDuplicates', 
						'INPUT={0}'.format(filtered_bam_file), 'OUTPUT={0}'.format(picard_output),
						 'METRICS_FILE={0}'.format(picard_metrics_file), 'VALIDATION_STRINGENCY=LENIENT',
						 'ASSUME_SORTED=true', 'REMOVE_DUPLICATES=false', 'VERBOSITY=WARNING']))

	logger.debug('Writing final BAM file for sample (with PCR duplicates removed): {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-F', '1804', '-f', '2', '-b', picard_output, '>', final_bam_file]))

	logger.debug('Indexing final BAM file for sample: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'index', final_bam_file]))

	logger.info('Finished filtering reads for sample: {0}'.format(sample_name))

	##################################
	# Generate additional QC reports
	##################################
	logger.debug('Computing read stats for sample: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'flagstat', final_bam_file, '>', final_bam_stats]))

	logger.info('Creating PBC QC file for BAM with PCR duplicates for sample: {0}'.format(sample_name))
	# PBC File output format:
		# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
	logger.debug('Creating temporary name-sorted filtered BAM file: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'sort', '-n', filtered_bam_file, filtered_bam_file_namesorted[:-4]]))

	logger.debug('Running PBC computation and writing report: {0}'.format(sample_name))
	cmd = "bamToBed -bedpe -i {0}".format(filtered_bam_file_namesorted) + \
		   " | grep -v 'chrM' | sort | uniq -c" + \
		   " | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'" + \
		   " > {0}".format(pbc_qc_file)
	os.system(cmd)

	logger.debug('Plotting insert size histogram for sample: {0}'.format(sample_name))
	os.system(' '.join(['CollectInsertSizeMetrics', 'INPUT={0}'.format(final_bam_file),
						'OUTPUT={0}'.format(insert_size_histogram_data), 'H={0}'.format(insert_size_histogram_plot)]))

	logger.info('Calculating percent mitochondrial reads for sample: {0}'.format(sample_name))

	logger.debug('Creating tagAlign formatted file for lightly-filtered aligned reads: {0}'.format(sample_name))
	cmd = "bamToBed -i {0}".format(early_filt_fixmate_output) + \
	 r"""| awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' """ + \
		"| gzip -nc > {0}".format(early_tagalign_file)
	os.system(cmd)

	# calculate total read pairs (before MQ filtering)
	p1 = Popen(['zcat', early_tagalign_file], stdout=PIPE)
	p2 = Popen(['wc', '-l'], stdin=p1.stdout, stdout=PIPE, stderr=PIPE)
	p1.stdout.close()
	output, err = p2.communicate()
	return_code = p2.wait()
	n_total_read_pairs = int(output) / 2
	# calculate number of mitochondrial read pairs (before MQ filtering)
	p1 = Popen(['zcat', early_tagalign_file], stdout=PIPE)
	p2 = Popen(['grep', 'chrM'], stdin=p1.stdout, stdout=PIPE)
	p3 = Popen(['wc', '-l'], stdin=p2.stdout, stdout=PIPE, stderr=PIPE)
	p1.stdout.close()
	p2.stdout.close()
	output, err = p3.communicate()
	return_code = p3.wait()
	n_mitochondrial_read_pairs = int(output) / 2
	with open(mitochondrial_read_report, 'w') as f:
		f.write('\t'.join(['TotalReadPairs', 'MitochondrialReadPairs', 'FracMitochondrialReads'])+'\n')
		f.write('\t'.join(map(str, [n_total_read_pairs, n_mitochondrial_read_pairs, 
						   float(n_mitochondrial_read_pairs)/float(n_total_read_pairs)]))+'\n')


	################################################
	# Convert PE BAM to tagAlign and BEDPE files
	################################################
	
	# Create "virtual single-ended file" 
	logger.debug('Creating tagAlign formatted file: {0}'.format(sample_name))
	cmd = "bamToBed -i {0}".format(final_bam_file) + \
	 r"""| awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' """ + \
		"| gzip -nc > {0}".format(final_tagalign_file)
	os.system(cmd)

	# Shift tagAlign file due to Tn5 insertion (original ATAC seq paper authors did this)
	#     "For peak-calling and footprinting, we adjusted the read start sites to represent the center of the transposon binding event. 
	#      Previous descriptions of the Tn5 transposase show that the transposon binds as a dimer and inserts two adaptors separated by 9 bp (ref. 11). 
	#      Therefore, all reads aligning to the + strand were offset by +4 bp, and all reads aligning to the - strand were offset -5 bp." 
	logger.info('Creating Tn5-shifted final tagAlign file: {0}'.format(sample_name))

	cmd = "zcat {0} ".format(final_tagalign_file) + \
	 r"""| awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4; $3 = $3 + 4} else if ($6 == "-") {$2 = $2 - 5; $3 = $3 - 5} print $0}' """ + \
		"| gzip -nc > {0}".format(final_Tn5shifted_tagAlign_file)
	os.system(cmd)

	logger.debug('Creating temporary name-sorted final BAM: {0}'.format(sample_name))
	os.system('samtools sort -n {0} {1}'.format(final_bam_file, final_bam_file_namesorted[:-4]))
	logger.info('Creating BEDPE file from name-sorted final BAM, excluding mitochondrial reads: {0}'.format(sample_name))
	os.system("bamToBed -bedpe -mate1 -i {0} | grep -v 'chrM' | gzip -nc > {1}".format(final_bam_file_namesorted, final_bedpe_file))

	logger.debug('Subsampling tagAlign file to {1} reads, excluding chrM: {0}'.format(sample_name, n_subsampled_reads))
	cmd = 'zcat {0} '.format(final_bedpe_file) + \
	'| shuf -n {0} --random-source={1} '.format(n_subsampled_reads, final_bedpe_file) + \
	r"""| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' """ + \
	'| gzip -nc > {0}'.format(subsampled_tagalign_file)
	os.system(cmd)

	def calc_and_write_frac_reads_close_to_TSS(final_tagalign_file, reference_tss_areas, output_file):
		p1 = Popen(['bedtools', 'intersect', '-u', '-a', final_tagalign_file, '-b', reference_tss_areas], stdout=PIPE)
		p2 = Popen(['grep', '-v', 'chrM'], stdin=p1.stdout, stdout=PIPE)
		p3 = Popen(['wc', '-l'], stdin=p2.stdout, stdout=PIPE, stderr=PIPE)
		p1.stdout.close()
		p2.stdout.close()
		output, err = p3.communicate()
		return_code = p3.wait()
		n_tss_reads = int(output)

		p1 = Popen(['zcat', final_tagalign_file], stdout=PIPE)
		p2 = Popen(['wc', '-l'], stdin=p1.stdout, stdout=PIPE, stderr=PIPE)
		p1.stdout.close()
		output, err = p2.communicate()
		return_code = p2.wait()
		n_total_reads = int(output)

		with open(output_file, 'w') as f:
			f.write('\t'.join(['NumReads(R1+R2)', 'NumReadsCloseToTSS(R1+R2)', 'FracTSS_Reads']) + '\n')
			f.write('\t'.join(map(str, [n_total_reads, n_tss_reads, float(n_tss_reads)/n_total_reads])) + '\n')

	logger.info('Calculating fraction of reads overlapping RefSeq TSS regions: {0}'.format(sample_name))
	calc_and_write_frac_reads_close_to_TSS(final_Tn5shifted_tagAlign_file, refseq_tss_areas, refseq_tss_report)
	logger.info('Calculating fraction of reads overlapping Gencode TSS regions: {0}'.format(sample_name))
	calc_and_write_frac_reads_close_to_TSS(final_Tn5shifted_tagAlign_file, gencode_tss_areas, gencode_tss_report)

	#################################
	# Call Peaks
	#################################

	logger.debug('Generating tagAlign file for MACS2 peak calling: {0}'.format(sample_name))
	fout = open(macs2_unsorted_tagAlign_file, 'w')
	random.seed(0)
	with gzip.open(final_bedpe_file, 'rt') as f:
		for line in f:
			fields    = line.strip().split('\t')
			r1_chr    = fields[0]
			r1_pos1   = int(fields[1])
			r1_pos2   = int(fields[2])
			r1_strand = fields[8]
			r2_chr    = fields[3]
			r2_pos1   = int(fields[4])
			r2_pos2   = int(fields[5])
			r2_strand = fields[9]
			assert(r1_chr == r2_chr)
			if(r1_strand == r2_strand):
				logger.warn('Final BEDPE contains read with unexpected R1-R2 orientation!: ' + line)
			#Tn5 shifting done by original ATAC seq paper authors
			if r1_strand == '+':
				r1_pos1 += 4
				r1_pos2 += 4
				r2_pos1 -= 5
				r2_pos2 -= 5
			elif r1_strand == '-':
				r1_pos1 -= 5
				r1_pos2 -= 5
				r2_pos1 += 4
				r2_pos2 += 4
			else:
				logger.warn('Final BEDPE final contains improperly-formatted line!: ' + line)
			posn_list = [r1_pos1, r1_pos2, r2_pos1, r2_pos2]
			min_posn = min(posn_list)
			max_posn = max(posn_list)
			if min_posn in [r1_pos1, r1_pos2]:
				min_strand = r1_strand
				max_strand = r2_strand
			elif min_posn in [r2_pos1, r2_pos2]:
				min_strand = r2_strand
				max_strand = r1_strand
			else:
				raise()
			# Pick only one of the Tn5 insertion points (ENCODE pipeline only keeps one in paired-end datasets, perhaps due to possibly breaking MACS2 statistical assumptions when including both?)
			if random.randint(0,1):
				fout.write('\t'.join(map(str, [r1_chr, min_posn, min_posn, 'N', 1000, min_strand])) + '\n')
			else:
				fout.write('\t'.join(map(str, [r1_chr, max_posn, max_posn, 'N', 1000, max_strand])) + '\n')
	fout.flush()
	fout.close()
	os.system(' '.join(['sort', '-k1,1', '-k2,2n', macs2_unsorted_tagAlign_file, '>', macs2_tagAlign_file[:-3]]))
	os.system('gzip -f {0}'.format(macs2_tagAlign_file[:-3]))

	# note: peak calling parameter set from Omni-ATAC paper: "macs2 callpeak --nomodel --nolambda --keep-dup all --call-summits"
	logger.info('Calling peaks with MACS2: {0}'.format(sample_name))

	for shiftsize, ext_size in zip(macs2_shift_sizes, macs2_smoothing_window_sizes):
		macs2_cvg_bedgraph_output      = '{0}/{2}/{1}_treat_pileup.bdg'.format(output_dir_sample, sample_name, 'macs2_output_windowSize{0}'.format(ext_size))
		macs2_tmp_sorted_bedgraph_file = '{0}/{2}/{1}.macs2.tmp.sorted.bedGraph'.format(output_dir_sample, sample_name, 'macs2_output_windowSize{0}'.format(ext_size))
		macs2_tmp_clipped_bedgraphFile = '{0}/{2}/{1}.macs2.tmp.sorted.clipped.bedGraph'.format(output_dir_sample, sample_name, 'macs2_output_windowSize{0}'.format(ext_size))
		macs2_smooth_bigWig_file       = '{0}/{2}/{1}.smoothed_Tn5insertionPoints.bigWig'.format(output_dir_sample, sample_name, 'macs2_output_windowSize{0}'.format(ext_size))
		[intermediate_file_list.append(x) for x in [macs2_cvg_bedgraph_output, macs2_tmp_sorted_bedgraph_file, macs2_tmp_clipped_bedgraphFile]]

		os.system(' '.join(['macs2', 'callpeak', '--nomodel', '--nolambda', '--keep-dup', 'all',
							'--call-summits', '-B', '--SPMR', '--format', 'BED',
							'-p', '0.01', '--shift', str(shiftsize), '--extsize', str(ext_size),
							'-t', macs2_tagAlign_file, '--outdir', output_dir_sample + '/macs2_output_windowSize{0}'.format(ext_size), 
							'--name', sample_name]))

		logger.debug('Sorting temporary bedgraph file from MACS2: {0}'.format(sample_name))
		os.system(' '.join(['sort', '-k1,1', '-k2,2n', macs2_cvg_bedgraph_output, '>', macs2_tmp_sorted_bedgraph_file]))

		logger.debug('Clipping peak segments larger than chromosome sizes: {0}'.format(sample_name))
		def make_chr_size_dict(chrom_sizes_file):
			chr_size_dict = {}
			with open(chrom_sizes_file) as f:
				for line in f:
					fields = line.strip().split('\t')
					chr_name = fields[0]
					chr_size = int(fields[1])
					chr_size_dict[chr_name] = chr_size
			return chr_size_dict
		chr_size_dict = make_chr_size_dict(chrom_sizes_file)

		fout = open(macs2_tmp_clipped_bedgraphFile, 'w')
		with open(macs2_tmp_sorted_bedgraph_file) as f:
			for line in f:
				fields    = line.strip().split('\t')
				chr_field = fields[0]
				coord1    = int(fields[1])
				coord2    = int(fields[2])
				chr_size  = chr_size_dict[chr_field]
				if coord1 >= chr_size and coord2 >= chr_size:
					continue
				elif coord2 >= chr_size:
					coord2 = chr_size - 1
				fout.write('\t'.join(map(str, [chr_field, coord1, coord2, fields[3]])) + '\n')
		fout.flush()
		fout.close()

		macs2_summits_file                  = '{0}/{1}/{2}_summits.bed'.format(output_dir_sample, 'macs2_output_windowSize{0}'.format(ext_size), sample_name)
		macs2_summits_blacklistFilteredfile = '{0}/{1}/{2}_summits.blacklistFiltered.bed'.format(output_dir_sample, 'macs2_output_windowSize{0}'.format(ext_size), sample_name)
		macs2_summits_blFilt_sorted         = '{0}/{1}/{2}_summits.blacklistFiltered.sorted.bed'.format(output_dir_sample, 'macs2_output_windowSize{0}'.format(ext_size), sample_name)
		macs2_summits_top300k               = '{0}/{1}/{2}_summits.top300k.bed'.format(output_dir_sample, 'macs2_output_windowSize{0}'.format(ext_size), sample_name)
		macs2_summits_top300k_posnSorted    = '{0}/{1}/{2}_summits.top300k.posnSorted.bed'.format(output_dir_sample, 'macs2_output_windowSize{0}'.format(ext_size), sample_name)
		macs2_summits_top300k_windows       = '{0}/{1}/{2}_summits.top300k.windowedRadius{3}.bed'.format(output_dir_sample, 'macs2_output_windowSize{0}'.format(ext_size), sample_name, macs2_summit_window_radius)

		logger.debug('Filtering blacklisted regions out of MACS2 summits: {0}'.format(sample_name))
		os.system(' '.join(['bedtools', 'intersect', '-v', '-sorted',
			                '-a', macs2_summits_file, '-b', hg38_blacklist,
			                '>', macs2_summits_blacklistFilteredfile]))

		#todo: rename peaks by rank?
		logger.debug('Selecting top {1} summits: {0}'.format(sample_name, macs2_cap_num_summits))
		os.system(' '.join(['sort', '-k', '5gr,5gr', macs2_summits_blacklistFilteredfile,
							'>', macs2_summits_blFilt_sorted]))
		os.system(' '.join(['head', '-n', str(macs2_cap_num_summits), macs2_summits_blFilt_sorted,
							'>', macs2_summits_top300k]))

		logger.debug('Sorting top MACS2 summits by genomic position: {0}'.format(sample_name))
		os.system(' '.join(['sort', '-k1,1', '-k2,2n', macs2_summits_top300k,
							'>', macs2_summits_top300k_posnSorted]))

		logger.debug('Converting summits to small genomic windows: {0}'.format(sample_name))
		os.system(' '.join(['cat', macs2_summits_top300k_posnSorted, '|', 'awk', 
						r"""'BEGIN{OFS="\t"}{""" + '$2 = $2 - {0}; $3 = $3 + {0}; print $0'.format(macs2_summit_window_radius) + r"}'",
						    '>', macs2_summits_top300k_windows]))


		[intermediate_file_list.append(x) for x in [macs2_summits_blacklistFilteredfile, macs2_summits_blFilt_sorted, macs2_summits_top300k]]

		logger.debug('Converting bedgraph to bigWig file (for viewing in IGV): {0}'.format(sample_name))
		os.system(' '.join(['bedGraphToBigWig', macs2_tmp_clipped_bedgraphFile, chrom_sizes_file, macs2_smooth_bigWig_file]))

	##################################
	# Final step: remove temporary intermediate files
	##################################
	logger.info('Deleting intermediate files: {0}'.format(sample_name))
	for f in intermediate_file_list:
		os.system(' '.join(['rm', f]))


if __name__ == '__main__':
    #args = _argparser().parse_args(sys.argv[1:])    
    logging.basicConfig(format='%(levelname)s: [%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S %Z')
    logger=logging.getLogger()
    logger.setLevel(logging.DEBUG)
    main()


#OLD CODE SNIPPETS

	# logger.info('Creating coverage bigWig tracks from final BAM file: {0}'.format(sample_name))
	# os.system(' '.join(['bedtools', 'genomecov', '-bg', '-ibam', final_bam_file, '>', final_tmp_bedgraph_file]))
	# logger.debug('Sorting temporary bedgraph file: {0}'.format(sample_name))
	# os.system(' '.join(['sort', '-k1,1', '-k2,2n', final_tmp_bedgraph_file, '>', final_tmp_sorted_bedgraph_file]))
	# logger.debug('Converting bedgraph to bigWig file: {0}'.format(sample_name))
	# os.system(' '.join(['bedGraphToBigWig', final_tmp_sorted_bedgraph_file, chrom_sizes_file, final_bigWig_file]))
	# logger.debug('Deleting temporary files used to make final coverage bigWig file: {0}'.format(sample_name))
	# os.system(' '.join(['rm', final_tmp_bedgraph_file]))
	# os.system(' '.join(['rm', final_tmp_sorted_bedgraph_file]))
