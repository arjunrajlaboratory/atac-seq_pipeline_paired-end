import sys
import os
import glob
import re

input_data_dir = '/home/esanford/data/HDD1/concatenated_data'
output_dir     = '/home/esanford/dev_atac_seq_pipeline/pipeline_outputs/HDD1'
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
	# print job_cmd
	os.system(job_cmd)



#Old code

# sample_names = """
# 01-EtOH
# 02-E2
# 03-RA
# 04-E2andRA
# 05-EtOH_to_EtOH
# 06-E2_to_E2
# 07-RA_to_RA
# 08-EtOH_to_E2andRA
# 09-E2_to_E2andRA
# 10-RA_to_E2andRA
# 11-E2andRA_to_E2andRA
# """

# sample_name_list = sample_names.split()
# remove already processed samples (subtract one for 0-based indexing):
# sample_name_list.pop(10 - 1)
# sample_name_list = [sample_name_list[6-1]]
