import sys
import logging

from subprocess import call
import os


# note: the call to this script should be preprended with a "wrapper" command that sets the paths for the appropriate versions of each tool called in this pipelin

#temporary paths for testing
bowtie2_index = '/Users/emsanford/Desktop/syd_atac_temp/hg19_downloaded_index/hg19'
picard_path = '/Users/emsanford/Desktop/Software/picard.jar'

sample_name = 'ENCLB034CND'
read1_path = '/Users/emsanford/Desktop/atac_dev_temp/ENCLB034CND/ENCFF581SDL.fastq.gz'
read2_path = '/Users/emsanford/Desktop/atac_dev_temp/ENCLB034CND/ENCFF157FMZ.fastq.gz'
sterr_log_file_step1a = '/Users/emsanford/Desktop/atac_dev_temp/stderrStep1a.txt'
step1a_output = '/Users/emsanford/Desktop/atac_dev_temp/testStep1_single_pipeline.bam'
step1b_input = step1a_output
step1b_tmp_filt_output = '/Users/emsanford/Desktop/atac_dev_temp/testStep1b.tmp.filt.bam'
step1b_tmp_filt_fixmate_output = '/Users/emsanford/Desktop/atac_dev_temp/testStep1b.tmp.filt.fixmate.bam'
step1b_filtered_bam_file = '/Users/emsanford/Desktop/atac_dev_temp/testStep1b.filt.bam'
picard_metrics_file = '/Users/emsanford/Desktop/atac_dev_temp/testStep1b.picardMetrics.txt'
picard_input = step1b_filtered_bam_file
picard_output = '/Users/emsanford/Desktop/atac_dev_temp/testStep1b.picardMarkedDups.bam'
final_bam_file = '/Users/emsanford/Desktop/atac_dev_temp/testStep1b.final.bam'
final_bam_stats = '/Users/emsanford/Desktop/atac_dev_temp/testStep1b.final.stats.txt'

#constants
bt2_multimapping = 4


def main():
	logger=logging.getLogger()

	# Step 1a: Align reads
	logger.info('Aligning reads for sample: {0}'.format(sample_name))

	os.system(' '.join(['bowtie2', '-k', str(bt2_multimapping), '--local', '-x', bowtie2_index, '-1', read1_path, '-2', read2_path, 
		  '2>', sterr_log_file_step1a, 
		  '|', 'samtools', 'view', '-Su', '/dev/stdin'
		  '|', 'samtools', 'sort', '-o', step1a_output]))

	logger.info('Finished aligning reads for sample: {0}'.format(sample_name))

	# Step 1b: Filter alignments
	logger.info('Filtering reads for sample: {0}'.format(sample_name))

	#samtools view -F 524 -f 2 -u ${RAW_BAM_FILE} | sambamba -n /dev/stdin -o ${TMP_FILT_BAM_FILE}
	logger.debug('Removing unmapped reads, reads with unmapped mates, or reads not passing vendor QC, then sorting by QNAME: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-F', '524', '-f', '2', '-u', step1b_input,
						'|', 'samtools', 'sort', '-n', '/dev/stdin', '-o', step1b_tmp_filt_output]))
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
	# It appears to unnecessarily discard reads, as "secondary alignments" are discarded by downstream steps and it discards all reads
	# when there are "k" candidate alignments.

	#samtools view -F 1804 -f 2 -u ${TMP_FILT_FIXMATE_BAM_FILE} | sambamba sort /dev/stdin -o ${FILT_BAM_FILE}
	logger.debug('Removing secondary alignments, optical duplicates, and sorting by position: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-F', '1804', '-f', '2', '-u', step1b_tmp_filt_fixmate_output,
						'|', 'samtools', 'sort', '/dev/stdin', '-o', step1b_filtered_bam_file]))
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
	os.system(' '.join(['java', '-Xmx4G', '-jar', picard_path, 'MarkDuplicates', 
						'INPUT={0}'.format(picard_input), 'OUTPUT={0}'.format(picard_output),
						 'METRICS_FILE={0}'.format(picard_metrics_file), 'VALIDATION_STRINGENCY=LENIENT',
						 'ASSUME_SORTED=true', 'REMOVE_DUPLICATES=false']))

	# samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}
	logger.debug('Writing final BAM file for sample (with PCR duplicatees removed): {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'view', '-F', '1804', '-f', '2', '-b', picard_output, '>', final_bam_file]))

	# # sambamba index ${FINAL_BAM_FILE}
	logger.debug('Indexing final BAM file for sample: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'index', final_bam_file]))

	# # samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}
	logger.debug('Counting read stats for sample: {0}'.format(sample_name))
	os.system(' '.join(['samtools', 'flagstat', final_bam_file, '>' ])
	# logger.info('Finished filtering reads for sample: {0}'.format(sample_name))

	#TODO remove unnecessary intermediate files


# # =============================
# # Compute library complexity
# # =============================
# # Sort by name
# # convert to bedPE and obtain fragment coordinates
# # sort by position and strand
# # Obtain unique count statistics

# module add bedtools/2.26

# PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

# # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair


# sambamba sort -n ${FILT_BAM_FILE} -o ${OFPREFIX}.srt.tmp.bam
# bedtools bamtobed -bedpe -i ${OFPREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
# rm ${OFPREFIX}.srt.tmp.bam

# rm ${FILT_BAM_FILE}
# # =============================

#def _argparser():

if __name__ == '__main__':
    #args = _argparser().parse_args(sys.argv[1:])    
    logging.basicConfig(format='%(levelname)s: [%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S %Z')
    logger=logging.getLogger()
    logger.setLevel(logging.DEBUG)
    main()