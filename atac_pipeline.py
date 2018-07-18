import sys
import logging

from subprocess import call
import os


# note: the call to this script should be preprended with a "wrapper" command that sets the paths for the appropriate versions of each tool called in this pipeline

#TODO remove unnecessary intermediate files
#TODO rename intermediate files with more descriptive names
#TODO add number of threads for bowtie?

#input data
# input_data_dir = '/home/esanford/data/history_dependence_pilot_18Jul2018_concatenated_lanes'
# sample_name = '10-RA_to_E2andRA'
# output_dir     = '/home/esanford/dev_atac_seq_pipeline/test_pipeline_outputs'
sample_name       = sys.argv[1]
input_data_dir    = sys.argv[2]
output_dir        = sys.argv[3]
n_bowtie2_threads = sys.argv[4]
read1_path     = '{0}/{1}'.format(input_data_dir, sample_name + '_R1.fastq.gz')
read2_path     = '{0}/{1}'.format(input_data_dir, sample_name + '_R2.fastq.gz')

#software references
refs_dir = '/home/esanford/refs'
# hg38 index
bowtie2_index = refs_dir + '/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'



#setup output filepaths
output_dir_sample  = output_dir + '/' + sample_name
if not os.path.exists(output_dir_sample):
	os.makedirs(output_dir_sample)

bowtie2_stderr_logfile         = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'log_stderr_bowtie2.txt')
bowtie2_aligned_reads          = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'bowtie2_aligned_reads.bam')

step1b_tmp_filt_output         = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'testStep1b.tmp.filt.bam')
step1b_tmp_filt_fixmate_output = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'testStep1b.tmp.filt.fixmate.bam')
step1b_filtered_bam_file       = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'testStep1b.filt.bam')
picard_metrics_file            = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'picardMetrics.txt')
picard_output                  = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'picardMarkedDups.bam')
final_bam_file                 = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.bam')
final_bam_stats                = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.stats.txt')
final_bam_file_namesorted      = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.nameSorted.bam')

final_tagalign_file            = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.tagAlign.gz')
final_bedpe_file               = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.bedpe.gz')

subsampled_tagalign_file       = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'subsampled.tagAlign.gz')

pbc_qc_file                    = '{0}/{1}.{2}'.format(output_dir_sample, sample_name, 'final.pbc.qc.txt')

#constants
bt2_multimapping    = 4
bt2_max_frag_len    = 5000
n_subsampled_reads  = 10000000


def main():
	logger=logging.getLogger()

	##################################
	# Step 1a: Align reads
	##################################

	logger.info('Aligning reads for sample: {0}'.format(sample_name))

	cmd = ' '.join(['bowtie2', '-k', str(bt2_multimapping), '--local', 
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
	# Step 1b: Filter alignments
	##################################
	logger.info('Filtering reads for sample: {0}'.format(sample_name))

	#samtools view -F 524 -f 2 -u ${RAW_BAM_FILE} | sambamba -n /dev/stdin -o ${TMP_FILT_BAM_FILE}
	logger.debug('Removing unmapped reads, reads with unmapped mates, or reads not passing vendor QC, then sorting by QNAME: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-F', '524', '-f', '2', '-u', bowtie2_aligned_reads,
						'|', 'samtools', 'sort', '-n', '/dev/stdin', step1b_tmp_filt_output[:-4]]))
	# samtools view: 
	#               -F 524   (-F is "blacklist" for certain flag bits)
	#					bit 512: "not passing filters, such as platform/vendor quality controls"
	#					bit 8: "next segment in the template unmapped"
	#					bit 4: "segment unmapped"
	#               -f 2     (-f is "whitelist" for certain flag bits) 
	#					bit 2: "each segment properly aligned according to the aligner"

	#samtools view -h ${TMP_FILT_BAM_FILE} | assign_multimappers.py -k $multimapping --paired-end | samtools fixmate -r /dev/stdin ${TMP_FILT_FIXMATE_BAM_FILE}
	#./set_atac_paths samtools view -h /Users/emsanford/Desktop/atac_dev_temp/testStep1b.tmp.filt.bam | python ./encode_utils/assign_multimappers.py -k 4 --paired-end | ./set_atac_paths samtools fixmate -r /dev/stdin assignMultimapTest_noChange
	logger.debug('Adding header, filtering out secondary alignments, and fixing mates: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-h', step1b_tmp_filt_output, 
				   '|', 'samtools', 'fixmate', '-r', '/dev/stdin', step1b_tmp_filt_fixmate_output]))

	#samtools view -F 1804 -f 2 -u ${TMP_FILT_FIXMATE_BAM_FILE} | sambamba sort /dev/stdin -o ${FILT_BAM_FILE}
	logger.debug('Removing low map-quality alignments, optical duplicates, and sorting by position: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-q', '20', '-F', '1804', '-f', '2', '-u', step1b_tmp_filt_fixmate_output,
						'|', 'samtools', 'sort', '/dev/stdin', step1b_filtered_bam_file[:-4]]))
	# samtools view: 
	#               -F 1804  (-F is "blacklist" for certain flag bits)
	#					bit 1024: pcr or optical duplicate 
	#					bit 512: "not passing filters, such as platform/vendor quality controls"
	#					bit 256: secondary alignment. "It is typically used to flag alternative mappings when multiple mappings are presented in a SAM"
	#					bit 8: "next segment in the template unmapped"
	#					bit 4: "segment unmapped"
	#               -f 2     (-f is "whitelist" for certain flag bits) 
	#					bit 2: "each segment properly aligned according to the aligner"


	# java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
	logger.debug('Marking duplicates with Picard for sample: {0}'.format(sample_name))
	os.system(' '.join(['MarkDuplicates', 
						'INPUT={0}'.format(step1b_filtered_bam_file), 'OUTPUT={0}'.format(picard_output),
						 'METRICS_FILE={0}'.format(picard_metrics_file), 'VALIDATION_STRINGENCY=LENIENT',
						 'ASSUME_SORTED=true', 'REMOVE_DUPLICATES=false', 'VERBOSITY=WARNING']))

	# samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}
	logger.debug('Writing final BAM file for sample (with PCR duplicates removed): {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-F', '1804', '-f', '2', '-b', picard_output, '>', final_bam_file]))

	# # sambamba index ${FINAL_BAM_FILE}
	logger.debug('Indexing final BAM file for sample: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'index', final_bam_file]))

	# samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}
	logger.debug('Counting read stats for sample: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'flagstat', final_bam_file, '>', final_bam_stats]))

	# TODO: replace this PBC QC file instructions with the instructions for paired end reads (this part I accidentally did for the single-end read pipeline)
	# bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
	logger.debug('Creating PBC QC file for BAM with PCR duplicates for sample: {0}'.format(sample_name))
	cmd = "bamToBed -i {0}".format(step1b_filtered_bam_file) + \
		   " | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$6}'" + \
		   " | grep -v 'chrM' | sort | uniq -c" + \
		   " | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'" + \
		   " > {0}".format(pbc_qc_file)
	os.system(cmd)
	# PBC File output format:
		# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
	logger.info('Finished filtering reads for sample: {0}'.format(sample_name))

	##################################
	# Step 2a: Convert PE BAM to tagAlign
	##################################
	
	# Create "virtual single-ended file" 
	logger.debug('Creating tagAlign formatted file: {0}'.format(sample_name))
	cmd = "bamToBed -i {0}".format(final_bam_file) + \
	 r"""| awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' """ + \
		"| gzip -nc > {0}".format(final_tagalign_file)
	os.system(cmd)

	#bedtools bamtobed -bedpe -mate1 -i ${FINAL_NMSRT_BAM_FILE} | gzip -nc > ${FINAL_BEDPE_FILE}
	logger.debug('Creating temporary name-sorted final BAM: {0}'.format(sample_name))
	os.system('samtools sort -n {0} {1}'.format(final_bam_file, final_bam_file_namesorted[:-4]))
	logger.debug('Creating BEDPE file from name-sorted final BAM: {0}'.format(sample_name))
	os.system("bamToBed -bedpe -mate1 -i {0} | gzip -nc > {1}".format(final_bam_file_namesorted, final_bedpe_file))
	logger.debug('Removing temporary name-sorted final BAM file: {0}'.format(sample_name))
	os.system('rm {0}'.format(final_bam_file_namesorted))

	# zcat ${FINAL_BEDPE_FILE} | grep -v "chrM" | shuf -n ${NREADS} --random-source=${FINAL_BEDPE_FILE}  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' | gzip -nc > ${SUBSAMPLED_TA_FILE}
	logger.debug('Subsampling tagAlign file to {1} reads, excluding chrM: {0}'.format(sample_name, n_subsampled_reads))
	cmd = 'zcat {0} '.format(final_bedpe_file) + \
	'| grep -v "chrM" ' + \
	'| shuf -n {0} --random-source={1} '.format(n_subsampled_reads, final_bedpe_file) + \
	r"""| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' """ + \
	'| gzip -nc > {0}'.format(subsampled_tagalign_file)
	print cmd
	os.system(cmd)

	logger.debug('Subsampling tagAlign file to {1} reads, excluding chrM: {0}'.format(sample_name, n_subsampled_reads))

#def _argparser():

if __name__ == '__main__':
    #args = _argparser().parse_args(sys.argv[1:])    
    logging.basicConfig(format='%(levelname)s: [%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S %Z')
    logger=logging.getLogger()
    logger.setLevel(logging.DEBUG)
    main()


#OLD CODE SNIPPETS

	# ###### BEGIN SUBPROCESS TESTING
	# from subprocess import Popen, PIPE
	# bt2_temp_aligned_reads = '/home/esanford/dev_atac_seq_pipeline/test_pipeline_outputs/ENCLB034CND/ENCLB034CND.bowtie2_alignments.sam'
	# # with open(bowtie2_aligned_reads, 'w') as fout:
	# p1 = Popen(['cat', bt2_temp_aligned_reads], stdout=PIPE)
	# p2 = Popen(['samtools', 'view', '-u', '/dev/stdin'], stdin=p1.stdout, stdout=PIPE)
	# p3 = Popen(['samtools', 'sort', '/dev/stdin', bowtie2_aligned_reads[:-4]], stdin=p2.stdout)
	# p1.stdout.close()
	# p2.stdout.close()
	# return_code = p3.wait()
	# print('Return code = {0}'.format(return_code))

	# # bt2_args = ['bowtie2', '-k', str(bt2_multimapping), '--local', '--maxins', str(bt2_max_frag_len), '-x', bowtie2_index, 
	# # 	        '-1', read1_path, '-2', read2_path, '2>', bowtie2_stderr_logfile]

	# ###### END SUBPROCESS TESTING