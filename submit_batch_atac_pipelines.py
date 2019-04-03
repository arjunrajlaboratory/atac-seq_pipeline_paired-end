import sys
import os
import glob
import re

# do not include terminal slash in directories
input_data_dir = sys.argv[1]
output_dir     = sys.argv[2]
num_bowtie2_threads = 4

sample_names = set() #start with set to disallow duplicates
for sample_filepath in glob.glob(input_data_dir + '/*'):
	sample_name_regex = '(.*/)(.*)(_R[1-2]).fastq.gz'
	re_match_obj = re.match(sample_name_regex, sample_filepath)
	sample_name = re_match_obj.group(2)
	sample_names.add(sample_name)
sample_name_list = list(sample_names)
print(sample_name_list)

for s in sample_name_list:
	job_cmd = 'bsub -J ' + 'pipeline_' + s + ' -n ' + str(num_bowtie2_threads) + ' ' \
				   '-o ' + '{0}/{1}/{1}.pipeline_stdout.txt'.format(output_dir, s) + ' ' \
				   '-e ' + '{0}/{1}/{1}.pipeline_stderr.txt'.format(output_dir, s) + ' ' \
				   '/home/esanford/dev_atac_seq_pipeline/set_atac_pmacs_env ' \
			       'python /home/esanford/dev_atac_seq_pipeline/atac_pipeline.py ' \
			       '{0} {1} {2} --num_bowtie2_threads {3}'.format(s, input_data_dir, output_dir, num_bowtie2_threads)
	print job_cmd
	os.system(job_cmd)
