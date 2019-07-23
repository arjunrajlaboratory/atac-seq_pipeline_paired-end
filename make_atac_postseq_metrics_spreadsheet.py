import glob
import re
import sys

#'/home/esanford/dev_atac_seq_pipeline/pipeline_outputs/HDD1'
input_data_dir = sys.argv[1]
output_file = sys.argv[2]

class Sample:
#lightweight class for sample paths, heavily borrowed/pasted from atac_pipeline.py
	def __init__(self, name, input_data_dir):
		self.name            = sample_name
		self.output_data_dir = input_data_dir + '/' + sample_name
		# fastqc files
		self.fastqc_r1 = '{0}/{1}_{2}'.format(self.output_data_dir, sample_name, 'R1_fastqc.html')
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

#just for making sample name list
sample_final_bams = glob.glob(input_data_dir + '/*/*.final.bam')
sample_objs = []
for s in sample_final_bams:

	sample_name_regex = '(.*/)(.*).final.bam'
	re_match_obj = re.match(sample_name_regex, s)
	sample_name = re_match_obj.group(2)

	sample_objs.append(Sample(sample_name, input_data_dir))


with open(output_file, 'w') as f:
#write header
	f.write('\t'.join(r'SampleName|initial # read pairs (pre-alignment and filtering)|final # read pairs (aligned, filtered, dups removed)|percent PCR/optical duplicates|est. library size (from Picard MarkDuplicates report)|PCR bottleneck coeff 1|mitochondrial read fraction|TSS read fraction (RefSeq)'.split('|')))
	f.write('\n')
	for s in sample_objs:

		total_num_reads = -1
		with open(s.fastqc_r1) as fr:
			re_num_total_reads = '.*<td>Total Sequences</td><td>([0-9]+)</td>.*'
			for line in fr:
				if re.match(re_num_total_reads, line):
					match_obj = re.match(re_num_total_reads, line)
					total_num_reads = int(match_obj.group(1))

		final_num_reads = -1
		with open(s.final_bam_stats) as fr:
			for line in fr.readlines():
				if 'read1' in line:
					final_num_reads = int(line.split()[0])

		frac_dups = -1
		est_libsize = -1
		with open(s.picard_metrics_file) as fr:
			for i in range(50): #sloppy code to avoid using readlines() method
				line = fr.readline()
				if line.startswith('## METRICS CLASS'):
					fr.readline() #discard header
					relevant_line = fr.readline()
					frac_dups = float(relevant_line.split('\t')[-2])
					est_libsize = int(relevant_line.split('\t')[-1])
					break

		# PBC File output format:
			# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
		pbc1 = -1
		with open(s.pbc_qc_file) as fr:
			pbc1 = float(fr.readline().split()[-2])

		frac_mito_reads = -1
		with open(s.mitochondrial_read_report) as fr:
			fr.readline() #discard header
			frac_mito_reads = float(fr.readline().split()[-1])

		frac_TSS_reads = -1
		with open(s.refseq_tss_report) as fr:
			fr.readline() #discard header
			frac_TSS_reads = float(fr.readline().split()[-1])


		f.write('\t'.join(map(str, [s.name, total_num_reads, final_num_reads, frac_dups, est_libsize, pbc1, frac_mito_reads, frac_TSS_reads]))+'\n')


