import glob
import sys
import os
import re
import time 

base_directory = '/home/esanford/data/HD3_ATAC-seq/data_from_illumina/FASTQ_Generation_2018-10-31_08_01_22Z-134203162'
output_directory = '/home/esanford/data/HD3_ATAC-seq/concatenated_data'
replicate = 'rep1'
fastq_files = glob.glob(base_directory + '/*/*.fastq*')
sample_dict = {}
number_of_lanes = 4
print fastq_files

for f in fastq_files:
	sample_dir_string = f.split('/')[-2]
	fastq_file_string = f.split('/')[-1]

	sample_name_regex = '(.*)(_S)([0-9]+)(_L00.*)'
	re_match_obj = re.match(sample_name_regex, fastq_file_string)
	sample_name = re_match_obj.group(1)
	sample_number = int(re_match_obj.group(3))

	print "{0} : {1}".format(sample_number, sample_name)
	#sample_key = "{0:02d}-{1}".format(sample_number, sample_name)
	sample_key = sample_name + '-' + replicate

	if sample_key not in sample_dict:
		sample_dict[sample_key] = []

	sample_dict[sample_key].append(f)


for sample, filepaths in sample_dict.items():
	read1_files = filter(lambda x: '_R1_' in x, filepaths)
	read2_files = filter(lambda x: '_R2_' in x, filepaths)

	print read1_files
	print read2_files

	assert(len(read1_files) == number_of_lanes)
	assert(len(read2_files) == number_of_lanes)

	r1_output_file = output_directory + '/' + sample + '_R1.fastq'
	cmd1 = 'zcat {0} > {1}'.format(' '.join(read1_files), r1_output_file)
	print cmd1
	os.system(cmd1)
	r2_output_file = output_directory + '/' + sample + '_R2.fastq'
	cmd2 = 'zcat {0} > {1}'.format(' '.join(read2_files), r2_output_file)
	print cmd2
	os.system(cmd2)

	# delay five seconds in case previous step isn't completely finished before proceeding to next step
	time.sleep(5)

	#compress files
	os.system('bsub gzip {0}'.format(r1_output_file))
	os.system('bsub gzip {0}'.format(r2_output_file))




