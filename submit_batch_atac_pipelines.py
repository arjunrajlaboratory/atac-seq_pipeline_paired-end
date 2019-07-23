import sys
import os
import glob
import re

# do not include terminal slash in directories
path_to_set_atac_pmacs_env = sys.argv[1]  # e.g., '/home/esanford/dev_atac_seq_pipeline/set_atac_pmacs_env'
path_to_atac_pipeline = sys.argv[2]       # e.g., '/home/esanford/dev_atac_seq_pipeline/atac_pipeline.py'
input_data_dir = sys.argv[3]
output_dir     = sys.argv[4]
num_bowtie2_threads = 4
memory_per_job_in_mb = 12 * 1024 #default on PMACS is 6GB. 6GB will cause most jobs to terminate if there are ~50-75M reads per sample.

sample_names = set() #start with set to disallow duplicates
for sample_filepath in glob.glob(input_data_dir + '/*'):
	try:
		sample_name_regex = '(.*/)(.*)(_R[1-2]).fastq.gz'
		re_match_obj = re.match(sample_name_regex, sample_filepath)
		sample_name = re_match_obj.group(2)
		sample_names.add(sample_name)
	except AttributeError:
		print('warning: {0} does not match expected filename for input fastq files! Expected filenames conform to this regular expression: {1}'.format(sample_filepath, sample_name_regex))

sample_name_list = list(sample_names)
print('about to submit jobs for these samples: {0}'.format(sample_name_list))
for s in sample_name_list:
	job_cmd = 'bsub -J ' + 'pipeline_' + s + ' -n ' + str(num_bowtie2_threads) + ' ' + \
				   '-M ' + str(memory_per_job_in_mb) + ' ' + \
				   '-o ' + '{0}/{1}/{1}.pipeline_stdout.txt'.format(output_dir, s) + ' ' + \
				   '-e ' + '{0}/{1}/{1}.pipeline_stderr.txt'.format(output_dir, s) + ' ' + \
				   path_to_set_atac_pmacs_env + ' ' + \
				   'python' + ' ' + \
				   path_to_atac_pipeline + ' ' + \
				   '{0} {1} {2} --num_bowtie2_threads {3}'.format(s, input_data_dir, output_dir, num_bowtie2_threads)
	print(job_cmd)
	os.system(job_cmd)
