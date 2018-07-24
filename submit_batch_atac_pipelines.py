import sys
import os

input_data_dir = '/home/esanford/data/history_dependence_pilot_18Jul2018_concatenated_lanes'
output_dir     = '/home/esanford/dev_atac_seq_pipeline/pipeline_outputs'
num_bowtie2_threads = 2

sample_names = """
01-EtOH
02-E2
03-RA
04-E2andRA
05-EtOH_to_EtOH
06-E2_to_E2
07-RA_to_RA
08-EtOH_to_E2andRA
09-E2_to_E2andRA
10-RA_to_E2andRA
11-E2andRA_to_E2andRA
"""
sample_name_list = sample_names.split()
# remove already processed samples (subtract one for 0-based indexing):
# sample_name_list.pop(10 - 1)
# sample_name_list = [sample_name_list[6-1]]

for s in sample_name_list:
	job_cmd = 'bsub -J ' + 'pipeline_' + s + ' -n ' + str(num_bowtie2_threads) + ' ' \
				   '-o ' + '{0}/{1}/{1}.pipeline_stdout.txt'.format(output_dir, s) + ' ' \
				   '-e ' + '{0}/{1}/{1}.pipeline_stderr.txt'.format(output_dir, s) + ' ' \
				   '/home/esanford/dev_atac_seq_pipeline/set_atac_pmacs_env ' \
			       'python /home/esanford/dev_atac_seq_pipeline/atac_pipeline.py ' \
			       '{0} {1} {2} {3}'.format(s, input_data_dir, output_dir, num_bowtie2_threads)
	# print job_cmd
	os.system(job_cmd)

