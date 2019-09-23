import sys
import logging
import argparse
from subprocess import Popen, PIPE
import os
import gzip
import random

#note: the call to this script should be preprended with a 'set_atac_pmacs_env' command that sets the paths for the appropriate versions of each tool called in this pipeline. 

#software references
refs_dir          = '/project/arjunrajlab/refs'
bowtie2_index     = refs_dir + '/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
gencode_tss_areas = refs_dir + '/atac_pipeline_refs/GENCODEv24_TSS_hg38_radius500.bed'
refseq_tss_areas  = refs_dir + '/atac_pipeline_refs/UCSC_RefSeq_TSS_hg38_radius500.bed'
chrom_sizes_file  = refs_dir + '/atac_pipeline_refs/hg38.chrom.sizes'
hg38_blacklist    = refs_dir + '/atac_pipeline_refs/hg38.blacklist.sorted.bed'


class ATAC_Sample:
	"""
	Lightweight class to store filepaths for atac-seq sample processing
	"""
	def __init__(self, sample_name, input_data_dir, output_dir, macs2_window_size, macs2_summit_window_radius):
		self.name            = sample_name
		self.input_data_dir  = input_data_dir
		self.output_data_dir = output_dir + '/' + sample_name
		# input fastq files (paired end)
		self.read1_path     = '{0}/{1}'.format(input_data_dir, sample_name + '_R1.fastq.gz')
		self.read2_path     = '{0}/{1}'.format(input_data_dir, sample_name + '_R2.fastq.gz')
		# bowtie2 alignment files
		self.bowtie2_stderr_logfile = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'log_stderr_bowtie2.txt')
		self.bowtie2_aligned_reads  = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'bowtie2_aligned_reads.bam')
		# alignment filtering files
		self.tmp_early_filt_output        = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'tmp.lightFilt.bam')
		self.early_filt_fixmate_output    = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'lightFilt.fixmate.bam')
		self.filtered_bam_file            = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'filt.bam')
		self.filtered_bam_file_namesorted = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'filt.nameSorted.bam')
		self.picard_metrics_file          = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'picardMetrics.txt')
		self.picard_output                = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'picardMarkedDups.bam')
		self.final_bam_file               = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'final.bam')
		self.final_bam_stats              = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'final.stats.txt')
		self.final_bam_file_namesorted    = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'final.nameSorted.bam')
		# text file representations of filtered aligned reads
		self.early_tagalign_file            = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'early.tagAlign.gz')
		self.final_tagalign_file            = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'final.tagAlign.gz')
		self.final_Tn5shifted_tagAlign_file = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'final.tagAlign_tn5_shifted.gz')
		self.final_bedpe_file               = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'final.bedpe.gz')
		self.macs2_unsorted_tagAlign_file   = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'macs2_formatted.unsorted.tagAlign')
		self.macs2_tagAlign_file            = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'Tn5_insertion_points.tagAlign.gz')
		# subsampled files
		self.subsampled_tagalign_file = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'subsampled.tagAlign.gz')
		# QC and summary report files
		self.insert_size_histogram_plot = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'insert_size_histogram.pdf')
		self.insert_size_histogram_data = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'insert_size_data.txt')
		self.mitochondrial_read_report  = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'mitochondrial_read_report.txt')
		self.refseq_tss_report          = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'refseq_tss_report.txt')
		self.gencode_tss_report         = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'gencode_tss_report.txt')
		self.pbc_qc_file                = '{0}/{1}.{2}'.format(self.output_data_dir, sample_name, 'final.pbc.qc.txt')
		# macs2 peak calling files
		macs2_folder_name = 'macs2_output_windowSize{0}'.format(macs2_window_size)
		self.macs2_output_dir = '{0}/{1}'.format(self.output_data_dir, macs2_folder_name)
		self.macs2_cvg_bedgraph_output      = '{0}/{2}/{1}_treat_pileup.bdg'.format(self.output_data_dir, sample_name, macs2_folder_name)
		self.macs2_tmp_sorted_bedgraph_file = '{0}/{2}/{1}.macs2.tmp.sorted.bedGraph'.format(self.output_data_dir, sample_name, macs2_folder_name)
		self.macs2_tmp_clipped_bedgraphFile = '{0}/{2}/{1}.macs2.tmp.sorted.clipped.bedGraph'.format(self.output_data_dir, sample_name, macs2_folder_name)
		self.macs2_smooth_bigWig_file       = '{0}/{2}/{1}.smoothed_Tn5insertionPoints.bigWig'.format(self.output_data_dir, sample_name, macs2_folder_name)
		#TODO
		self.macs2_narrowPeak_file               = '{0}/{1}/{2}_peaks.narrowPeak'.format(self.output_data_dir, macs2_folder_name, sample_name)
		self.macs2_narrowPeak_fullChrsOnly       = '{0}/{1}/{2}_fullChrsOnly.narrowPeak'.format(self.output_data_dir, macs2_folder_name, sample_name)
		self.macs2_narrowPeak_blacklistFiltered  = '{0}/{1}/{2}_blacklistFiltered.narrowPeak'.format(self.output_data_dir, macs2_folder_name, sample_name)
		self.macs2_summits_file                  = '{0}/{1}/{2}_summits.bed'.format(self.output_data_dir, macs2_folder_name, sample_name)
		#TODO
		self.macs2_summits_file_fullChrsOnly     = '{0}/{1}/{2}_summits.fullChrsOnly.bed'.format(self.output_data_dir, macs2_folder_name, sample_name)
		self.macs2_summits_blacklistFilteredfile = '{0}/{1}/{2}_summits.blacklistFiltered.bed'.format(self.output_data_dir, macs2_folder_name, sample_name)
		#TODO
		self.macs2_summits_blFilt_sorted         = '{0}/{1}/{2}_summits.blacklistFiltered.sorted.bed'.format(self.output_data_dir, macs2_folder_name, sample_name)
		self.macs2_summits_top300k               = '{0}/{1}/{2}_summits.top300k.bed'.format(self.output_data_dir, macs2_folder_name, sample_name)
		self.macs2_summits_top300k_posnSorted    = '{0}/{1}/{2}_summits.top300k.posnSorted.bed'.format(self.output_data_dir, macs2_folder_name, sample_name)
		self.macs2_summits_top300k_windows       = '{0}/{1}/{2}_summits.top300k.windowedRadius{3}.bed'.format(self.output_data_dir, macs2_folder_name, sample_name, macs2_summit_window_radius)

		#files to be deleted after pipeline runs
		self.intermediate_file_list = [self.tmp_early_filt_output, self.early_filt_fixmate_output, self.early_tagalign_file, self.filtered_bam_file_namesorted, 
								       self.filtered_bam_file, self.picard_output, self.final_bam_file_namesorted, self.macs2_unsorted_tagAlign_file,
								       self.macs2_cvg_bedgraph_output, self.macs2_tmp_sorted_bedgraph_file, self.macs2_tmp_clipped_bedgraphFile, self.macs2_summits_file_fullChrsOnly,
								       self.macs2_narrowPeak_fullChrsOnly, self.macs2_summits_blFilt_sorted, self.macs2_summits_top300k]


def run_command(cmd):
	logger = logging.getLogger()
	logger.debug(cmd)
	os.system(cmd)


def setup_output_filepaths(s_obj):
	if not os.path.exists(s_obj.output_data_dir):
		os.makedirs(s_obj.output_data_dir)


def run_fastQC(s_obj):
	run_command(' '.join(['fastqc', s_obj.read1_path, '--outdir', s_obj.output_data_dir]))
	run_command(' '.join(['fastqc', s_obj.read2_path, '--outdir', s_obj.output_data_dir]))


def bowtie2_align_reads(s_obj, bt2_max_frag_len, num_bowtie2_threads):
	cmd = ' '.join(['bowtie2', '--local', 
		           '--maxins', str(bt2_max_frag_len), '-x', bowtie2_index,
		           '--threads', str(num_bowtie2_threads), 
		           '-1', s_obj.read1_path, '-2', s_obj.read2_path, 
		           '2>', s_obj.bowtie2_stderr_logfile, 
		      '|', 'samtools', 'view', '-u', '/dev/stdin',	
		      '|', 'samtools', 'sort', '/dev/stdin', s_obj.bowtie2_aligned_reads[:-4]]) # [:-4] is necessary because samtools adds a bam suffix
	run_command(cmd)	


def filter_alignments(s_obj):
	logger = logging.getLogger()
	logger.info('Removing unmapped reads, reads with unmapped mates, or reads not passing vendor QC, then sorting by QNAME: {0}'.format(s_obj.name))
	run_command(' '.join(['samtools', 'view', '-F', '524', '-f', '2', '-u', s_obj.bowtie2_aligned_reads,
						'|', 'samtools', 'sort', '-n', '/dev/stdin', s_obj.tmp_early_filt_output[:-4]]))
	# samtools view: 
	#               -F 524   (-F is "blacklist" for certain flag bits)
	#					bit 512: "not passing filters, such as platform/vendor quality controls"
	#					bit 8: "next segment in the template unmapped"
	#					bit 4: "segment unmapped"
	#               -f 2     (-f is "whitelist" for certain flag bits) 
	#					bit 2: "each segment properly aligned according to the aligner"

	logger.info('Adding header, filtering out secondary alignments, and fixing mates: {0}'.format(s_obj.name))
	run_command(' '.join(['samtools', 'view', '-h', s_obj.tmp_early_filt_output, 
				   '|', 'samtools', 'fixmate', '-r', '/dev/stdin', s_obj.early_filt_fixmate_output]))
	logger.info('Removing low map-quality alignments, optical duplicates, and sorting by position: {0}'.format(s_obj.name))
	run_command(' '.join(['samtools', 'view', '-q', '20', '-F', '1804', '-f', '2', '-u', s_obj.early_filt_fixmate_output,
						'|', 'samtools', 'sort', '/dev/stdin', s_obj.filtered_bam_file[:-4]]))
	# samtools view: 
	#               -F 1804  (-F is "blacklist" for certain flag bits)
	#					bit 1024: pcr or optical duplicate 
	#					bit 512: "not passing filters, such as platform/vendor quality controls"
	#					bit 256: secondary alignment. "It is typically used to flag alternative mappings when multiple mappings are presented in a SAM"
	#					bit 8: "next segment in the template unmapped"
	#					bit 4: "segment unmapped"
	#               -f 2     (-f is "whitelist" for certain flag bits) 
	#					bit 2: "each segment properly aligned according to the aligner"
	logger.info('Marking duplicates with Picard for sample: {0}'.format(s_obj.name))
	run_command(' '.join(['MarkDuplicates', 
						'INPUT={0}'.format(s_obj.filtered_bam_file), 'OUTPUT={0}'.format(s_obj.picard_output),
						 'METRICS_FILE={0}'.format(s_obj.picard_metrics_file), 'VALIDATION_STRINGENCY=LENIENT',
						 'ASSUME_SORTED=true', 'REMOVE_DUPLICATES=false', 'VERBOSITY=WARNING']))
	logger.info('Writing final BAM file for sample (with PCR duplicates removed): {0}'.format(s_obj.name))
	run_command(' '.join(['samtools', 'view', '-F', '1804', '-f', '2', '-b', s_obj.picard_output, '>', s_obj.final_bam_file]))
	logger.info('Indexing final BAM file for sample: {0}'.format(s_obj.name))
	run_command(' '.join(['samtools', 'index', s_obj.final_bam_file]))


def generate_summary_reports(s_obj):
	logger = logging.getLogger()

	logger.info('Computing samtools stats for sample: {0}'.format(s_obj.name))
	run_command(' '.join(['samtools', 'flagstat', s_obj.final_bam_file, '>', s_obj.final_bam_stats]))

	logger.info('Creating PBC QC file for BAM with PCR duplicates for sample: {0}'.format(s_obj.name))
	# PBC File output format:
		# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
	logger.info('Creating temporary name-sorted filtered BAM file: {0}'.format(s_obj.name))
	run_command(' '.join(['samtools', 'sort', '-n', s_obj.filtered_bam_file, s_obj.filtered_bam_file_namesorted[:-4]]))

	logger.info('Running PBC computation and writing report: {0}'.format(s_obj.name))
	cmd = "bamToBed -bedpe -i {0}".format(s_obj.filtered_bam_file_namesorted) + \
		   " | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$6,$9,$10}'" + \
		   " | grep -v 'chrM' | sort | uniq -c" + \
		   " | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'" + \
		   " > {0}".format(s_obj.pbc_qc_file)
	run_command(cmd)

	logger.info('Plotting insert size histogram for sample: {0}'.format(s_obj.name))
	run_command(' '.join(['CollectInsertSizeMetrics', 'INPUT={0}'.format(s_obj.final_bam_file),
						'OUTPUT={0}'.format(s_obj.insert_size_histogram_data), 'H={0}'.format(s_obj.insert_size_histogram_plot)]))

	logger.info('Calculating percent mitochondrial reads for sample: {0}'.format(s_obj.name))

	logger.info('Creating tagAlign formatted file for lightly-filtered aligned reads: {0}'.format(s_obj.name))
	cmd = "bamToBed -i {0}".format(s_obj.early_filt_fixmate_output) + \
	 r"""| awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' """ + \
		"| gzip -nc > {0}".format(s_obj.early_tagalign_file)
	run_command(cmd)

	# calculate total read pairs (before MQ filtering)
	p1 = Popen(['zcat', s_obj.early_tagalign_file], stdout=PIPE)
	p2 = Popen(['wc', '-l'], stdin=p1.stdout, stdout=PIPE, stderr=PIPE)
	p1.stdout.close()
	output, err = p2.communicate()
	return_code = p2.wait()
	n_total_read_pairs = int(output) / 2
	# calculate number of mitochondrial read pairs (before MQ filtering)
	p1 = Popen(['zcat', s_obj.early_tagalign_file], stdout=PIPE)
	p2 = Popen(['grep', 'chrM'], stdin=p1.stdout, stdout=PIPE)
	p3 = Popen(['wc', '-l'], stdin=p2.stdout, stdout=PIPE, stderr=PIPE)
	p1.stdout.close()
	p2.stdout.close()
	output, err = p3.communicate()
	return_code = p3.wait()
	n_mitochondrial_read_pairs = int(output) / 2
	with open(s_obj.mitochondrial_read_report, 'w') as f:
		f.write('\t'.join(['TotalReadPairs', 'MitochondrialReadPairs', 'FracMitochondrialReads'])+'\n')
		f.write('\t'.join(map(str, [n_total_read_pairs, n_mitochondrial_read_pairs, 
						   float(n_mitochondrial_read_pairs)/float(n_total_read_pairs)]))+'\n')


def generate_text_alignment_formats(s_obj):
	logger = logging.getLogger()
	# Create "virtual single-ended file" 
	logger.info('Creating tagAlign formatted file: {0}'.format(s_obj.name))
	cmd = "bamToBed -i {0}".format(s_obj.final_bam_file) + \
	 r"""| awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' """ + \
		"| gzip -nc > {0}".format(s_obj.final_tagalign_file)
	run_command(cmd)
	# Shift tagAlign file due to Tn5 insertion (original ATAC seq paper authors did this)
	#     "For peak-calling and footprinting, we adjusted the read start sites to represent the center of the transposon binding event. 
	#      Previous descriptions of the Tn5 transposase show that the transposon binds as a dimer and inserts two adaptors separated by 9 bp (ref. 11). 
	#      Therefore, all reads aligning to the + strand were offset by +4 bp, and all reads aligning to the - strand were offset -5 bp." 
	logger.info('Creating Tn5-shifted final tagAlign file: {0}'.format(s_obj.name))
	cmd = "zcat {0} ".format(s_obj.final_tagalign_file) + \
	 r"""| awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4; $3 = $3 + 4} else if ($6 == "-") {$2 = $2 - 5; $3 = $3 - 5} print $0}' """ + \
		"| gzip -nc > {0}".format(s_obj.final_Tn5shifted_tagAlign_file)
	run_command(cmd)
	logger.info('Creating temporary name-sorted final BAM: {0}'.format(s_obj.name))
	run_command('samtools sort -n {0} {1}'.format(s_obj.final_bam_file, s_obj.final_bam_file_namesorted[:-4]))
	logger.info('Creating BEDPE file from name-sorted final BAM, excluding mitochondrial reads: {0}'.format(s_obj.name))
	run_command("bamToBed -bedpe -mate1 -i {0} | grep -v 'chrM' | gzip -nc > {1}".format(s_obj.final_bam_file_namesorted, s_obj.final_bedpe_file))

	logger.info('Generating tagAlign file for MACS2 peak calling: {0}'.format(s_obj.name))
	fout = open(s_obj.macs2_unsorted_tagAlign_file, 'w')
	random.seed(0)
	with gzip.open(s_obj.final_bedpe_file, 'rt') as f:
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
	run_command(' '.join(['sort', '-k1,1', '-k2,2n', s_obj.macs2_unsorted_tagAlign_file, '>', s_obj.macs2_tagAlign_file[:-3]]))
	run_command('gzip -f {0}'.format(s_obj.macs2_tagAlign_file[:-3]))


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


def macs2_call_peaks(s_obj, shiftsize, ext_size, macs2_cap_num_summits, macs2_summit_window_radius):
	logger = logging.getLogger()

	run_command(' '.join(['macs2', 'callpeak', '--nomodel', '--nolambda', '--keep-dup', 'all',
						'--call-summits', '-B', '--SPMR', '--format', 'BED',
						'-q', '0.05', '--shift', str(shiftsize), '--extsize', str(ext_size),
						'-t', s_obj.macs2_tagAlign_file, '--outdir', s_obj.macs2_output_dir, 
						'--name', s_obj.name]))

	logger.info('Sorting temporary bedgraph file from MACS2: {0}'.format(s_obj.name))
	run_command(' '.join(['sort', '-k1,1', '-k2,2n', s_obj.macs2_cvg_bedgraph_output, '>', s_obj.macs2_tmp_sorted_bedgraph_file]))

	logger.info('Clipping peak segments larger than chromosome sizes: {0}'.format(s_obj.name))
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

	fout = open(s_obj.macs2_tmp_clipped_bedgraphFile, 'w')
	with open(s_obj.macs2_tmp_sorted_bedgraph_file) as f:
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

	fullChrOnly_grep_command = r'grep -P "^(chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY)\t"'
	logger.info('Filtering non-full-chr contigs out of MACS2 narrowPeak file: {0}'.format(s_obj.name))
	run_command(' '.join([fullChrOnly_grep_command, s_obj.macs2_narrowPeak_file, ">", s_obj.macs2_narrowPeak_fullChrsOnly]))
	logger.info('Filtering blacklisted regions out of MACS2 narrowPeak file: {0}'.format(s_obj.name))
	run_command(' '.join(['bedtools', 'intersect', '-v', '-sorted',
		                '-a', s_obj.macs2_narrowPeak_fullChrsOnly, '-b', hg38_blacklist,
		                '>', s_obj.macs2_narrowPeak_blacklistFiltered]))
	logger.info('Filtering non-full-chr contigs out of MACS2 summits: {0}'.format(s_obj.name))
	run_command(' '.join([fullChrOnly_grep_command, s_obj.macs2_summits_file, ">", s_obj.macs2_summits_file_fullChrsOnly]))
	logger.info('Filtering blacklisted regions out of MACS2 summits: {0}'.format(s_obj.name))
	run_command(' '.join(['bedtools', 'intersect', '-v', '-sorted',
		                '-a', s_obj.macs2_summits_file_fullChrsOnly, '-b', hg38_blacklist,
		                '>', s_obj.macs2_summits_blacklistFilteredfile]))
	#optional todo: rename peaks by rank?
	logger.info('Selecting top {1} summits: {0}'.format(s_obj.name, macs2_cap_num_summits))
	run_command(' '.join(['sort', '-k', '5gr,5gr', s_obj.macs2_summits_blacklistFilteredfile,
						'>', s_obj.macs2_summits_blFilt_sorted]))
	run_command(' '.join(['head', '-n', str(macs2_cap_num_summits), s_obj.macs2_summits_blFilt_sorted,
						'>', s_obj.macs2_summits_top300k]))
	logger.info('Sorting top MACS2 summits by genomic position: {0}'.format(s_obj.name))
	run_command(' '.join(['sort', '-k1,1', '-k2,2n', s_obj.macs2_summits_top300k,
						'>', s_obj.macs2_summits_top300k_posnSorted]))
	logger.info('Converting summits to small genomic windows: {0}'.format(s_obj.name))
	run_command(' '.join(['cat', s_obj.macs2_summits_top300k_posnSorted, '|', 'awk', 
					r"""'BEGIN{OFS="\t"}{""" + '$2 = $2 - {0}; $3 = $3 + {0}; print $0'.format(macs2_summit_window_radius) + r"}'",
					    '>', s_obj.macs2_summits_top300k_windows]))
	logger.info('Converting bedgraph to bigWig file (for viewing in IGV): {0}'.format(s_obj.name))
	run_command(' '.join(['bedGraphToBigWig', s_obj.macs2_tmp_clipped_bedgraphFile, chrom_sizes_file, s_obj.macs2_smooth_bigWig_file]))


def _arg_parser():
	script_usage = script_usage = '''%(prog)s [options] sample_name input_data_directory output_data_directory'''
	parser = argparse.ArgumentParser(usage=script_usage)

	g = parser.add_argument_group("Required Arguments")
	g.add_argument('sample_name', type=str, help='The sample name')
	g.add_argument('input_data_directory', type=str,
	                help='Folder containing input fastq files that correspond to the name of the sample')
	g.add_argument('output_data_directory', type=str,
                    help='Folder in which to store pipeline outputs')

	g = parser.add_argument_group("Parameters")
	g.add_argument('--num_bowtie2_threads', type=int, default=1,
    	           help='Number of threads for running the bowtie2 aligner in parallel')
	g.add_argument('--bowtie2_max_fragment_length', type=int, default=2000,
				   help="Maximum fragment size to consider for aligning paired end reads")
	# note: in encode atac specification, for "ENCODE3" the window size is 73 and the shift size is -37. However 150 is the default in Kundaje lab pipeline.
	g.add_argument('--macs2_smoothing_window_size', type=int, default=150,
			       help="The width of the window for smoothing inferred Tn5 insertion points, used by MACS2 peak caller for calling peaks and creating bigWig files")
	g.add_argument('--macs2_summit_window_radius', type=int, default=50,
				   help="Hard-coded peak size for summit-based analysis")
	g.add_argument('--macs2_cap_num_summits', type=int, default=300000,
    			   help="Upper bound for the number of MACS2 summits used in downstream analysis")
	return parser


def main(args):
	logging.basicConfig(format='%(levelname)s: [%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S %Z')
	logger=logging.getLogger()
	logger.setLevel(logging.DEBUG)

	sample_name         = args.pop('sample_name')
	input_data_dir      = args.pop('input_data_directory')
	output_dir          = args.pop('output_data_directory')
	num_bowtie2_threads = args.pop('num_bowtie2_threads')
	bt2_max_frag_len    = args.pop('bowtie2_max_fragment_length')
	macs2_smoothing_window_size = args.pop('macs2_smoothing_window_size')
	macs2_shift_size            = int(round(macs2_smoothing_window_size/2) * -1) #shift size should be int, -1/2 * corresponding smoothing window sizes
	macs2_cap_num_summits       = args.pop('macs2_cap_num_summits')
	macs2_summit_window_radius  = args.pop('macs2_summit_window_radius')

	s_obj = ATAC_Sample(sample_name, input_data_dir, output_dir, macs2_smoothing_window_size, macs2_summit_window_radius)
	setup_output_filepaths(s_obj)

	logger.info('Running FastQC for sample: {0}'.format(s_obj.name))
	run_fastQC(s_obj)
	logger.info('Aligning reads for sample: {0}'.format(s_obj.name))
	bowtie2_align_reads(s_obj, bt2_max_frag_len, num_bowtie2_threads)
	logger.info('Filtering alignments for sample: {0}'.format(s_obj.name))
	filter_alignments(s_obj)
	logger.info('Generating summary and QC reports for sample: {0}'.format(s_obj.name))
	generate_summary_reports(s_obj)
	logger.info('Generating text-based alignment file formats: {0}'.format(s_obj.name))
	generate_text_alignment_formats(s_obj)
	logger.info('Calculating fraction of reads overlapping RefSeq TSS regions: {0}'.format(s_obj.name))
	calc_and_write_frac_reads_close_to_TSS(s_obj.final_Tn5shifted_tagAlign_file, refseq_tss_areas, s_obj.refseq_tss_report)
	logger.info('Calculating fraction of reads overlapping Gencode TSS regions: {0}'.format(s_obj.name))
	calc_and_write_frac_reads_close_to_TSS(s_obj.final_Tn5shifted_tagAlign_file, gencode_tss_areas, s_obj.gencode_tss_report)
	logger.info('Calling peaks with MACS2: {0}'.format(s_obj.name))
	## note: peak calling parameter set from Omni-ATAC paper is: "macs2 callpeak --nomodel --nolambda --keep-dup all --call-summits"
	macs2_call_peaks(s_obj, macs2_shift_size, macs2_smoothing_window_size, macs2_cap_num_summits, macs2_summit_window_radius)
	logger.info('Deleting intermediate files: {0}'.format(s_obj.name))
	for f in s_obj.intermediate_file_list:
		run_command(' '.join(['rm', f]))


if __name__ == '__main__':
    args = vars(_arg_parser().parse_args(sys.argv[1:]))
    main(args)


#OLD CODE SNIPPETS

# #TODO: find purpose for subsampled files... do pseudoreplicate peak calling and IDR! (use half of reads)
# n_subsampled_reads  = 10000000
# logger.debug('Subsampling tagAlign file to {1} reads, excluding chrM: {0}'.format(sample_name, n_subsampled_reads))
# cmd = 'zcat {0} '.format(final_bedpe_file) + \
# '| shuf -n {0} --random-source={1} '.format(n_subsampled_reads, final_bedpe_file) + \
# r"""| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' """ + \
# '| gzip -nc > {0}'.format(subsampled_tagalign_file)
# os.system(cmd)

	# logger.info('Creating coverage bigWig tracks from final BAM file: {0}'.format(sample_name))
	# os.system(' '.join(['bedtools', 'genomecov', '-bg', '-ibam', final_bam_file, '>', final_tmp_bedgraph_file]))
	# logger.debug('Sorting temporary bedgraph file: {0}'.format(sample_name))
	# os.system(' '.join(['sort', '-k1,1', '-k2,2n', final_tmp_bedgraph_file, '>', final_tmp_sorted_bedgraph_file]))
	# logger.debug('Converting bedgraph to bigWig file: {0}'.format(sample_name))
	# os.system(' '.join(['bedGraphToBigWig', final_tmp_sorted_bedgraph_file, chrom_sizes_file, final_bigWig_file]))
	# logger.debug('Deleting temporary files used to make final coverage bigWig file: {0}'.format(sample_name))
	# os.system(' '.join(['rm', final_tmp_bedgraph_file]))
	# os.system(' '.join(['rm', final_tmp_sorted_bedgraph_file]))
