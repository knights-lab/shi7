from __future__ import print_function
import argparse
import subprocess
import os
import re
import shutil
import multiprocessing
from datetime import datetime

#gcc -m64 -O3 ninja_shi7.c -o ninja_shi7_mac
#gcc -m64 -O3 -msse4.1 ninja_shi7.c -o ninja_shi7_mac
#subprocess.call("for var in 0 1 2; do ls ..; done;", shell=True)
#subprocess.call("cd sample_fastq; ./runTM.cmd", shell=True)
startTime = datetime.now()

parser = argparse.ArgumentParser(description='This is the commandline interface for NINJA-SHI7', usage='python ninja_shi7.py -i <input> -o <output> -t_trim <threads>...')

parser.add_argument('-adaptor', '--adaptor_type', help='Set the type of the adaptor (default: %(default)s)', choices=['None','Nextera','TruSeq3'], default='None')
parser.add_argument('-f', '--flash', help='Enable or disable the FLASH (default: %(default)s)', choices=['enabled','disabled'], default='enabled')
parser.add_argument('-trim','--trimmer', help='Enable or disable the TRIMMER (default: %(default)s)', choices=['enabled','disabled'], default='enabled')
parser.add_argument('-fasta', '--make_fastas', help='Enable or disable the MAKE_FASTAS (default: %(default)s)', choices=['enabled','disabled'], default='enabled')
parser.add_argument('-qiime', '--QIIME', help='Enable or disable the QIIME (default: %(default)s) NOTE: MAKE_FASTAS needs to be enabled in order to run QIIME.', choices=['enabled','disabled'], default='enabled')
parser.add_argument('-append', '--append_fna', help='Enable or disable the fna append mode (default: %(default)s)', choices=['enabled','disabled'], default='disabled')

parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory', required=True)
parser.add_argument('-o', '--output', help='Set the directory path of the destination (default: the directory path where NINJA-SHI7 is located)')
parser.add_argument('-t', '--threads', help='Set the number of threads (default: %(default)s)', default=str(multiprocessing.cpu_count()))
parser.add_argument('-O', '--allow_outies', help='Enable "outie" orientation Choose one option: enabled/disabled (deafult: %(default)s)', default='enabled')
parser.add_argument('-m', '--min_overlap', help='Set the minimum overlap length between two reads (default: %(default)s)', default='20')
parser.add_argument('-M', '--max_overlap', help='Set the maximum overlap length between two reads (default: %(default)s)', default='700')
parser.add_argument('-trim_l', '--trim_length', help='Set the trim length (default: %(default)s)', default='150')
parser.add_argument('-trim_q', '--trim_qual', help='Set the trim qual (default: %(default)s)', default='20')
args = parser.parse_args()
#print args.threads_trimmomatic, args.min_overlap, args.max_overlap, args.input, args.allow_outies

script_path = os.path.dirname(os.path.realpath(__file__))
print('script_path:', script_path)
fastq_path = args.input
print('fastq_path:', fastq_path)
#cd_fastq_path = 'cd /Users/KaiweiAng/Documents/Knights_Lab/NINJA-SHI7/sample_fastq; ' #input from user
#cd_fastq_path = 'cd ' + fastq_path + '; '

#FIRST CHECK IF THE INPUT AND OUTPUT PATH EXIST. IF DO NOT, RAISE EXCEPTION AND EXIT
if os.path.exists(fastq_path) == False:
	print('Error:', fastq_path, 'doesn\'t exist!')
	exit()

if args.output != None and os.path.exists(os.path.dirname(args.output)) == False: 
	print('Error:', os.path.dirname(args.output), 'doesn\'t exist!')
	exit()

if args.output != None and args.output.endswith('.fna') == False:
	print('Error: output file must be a .fna file!')
	exit() 

#copy the fastq samples to a temp folder in the script directory and use it instead of the original fastq folder
if os.path.exists(fastq_path + '/temp'):
	shutil.rmtree(fastq_path + '/temp')
	print('Existing temp directory deleted.')
#subprocess.call('cp -r ' + fastq_path + ' ' + fastq_path + '/temp', shell=True)
shutil.copytree(fastq_path, fastq_path + '/temp')

#check if MAKE_FASTAS is enabled when QIMME is enabled
if args.QIIME == 'enabled' and args.make_fastas == 'disabled':
	print('MAKE_FASTAS needs to be ENABLED in order to run QIIME!')
	exit()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#AXE_ADAPTORS
FASTQ = [f for f in os.listdir(fastq_path + '/temp') if f.endswith('fastq')]
print(FASTQ)
R1_fq = [f for f in FASTQ if re.search('R1',f) != None]
R2_fq = [f for f in FASTQ if re.search('R2',f) != None]
print(R1_fq)
print(R2_fq)
if len(R1_fq) != len(R2_fq) or len(R1_fq) < 1:
	print('Error:', fastq_path + ':', 'The input directory must contain at least one pair of R1 & R2 fastq file!')
	exit()

if args.adaptor_type == 'None':
	for f in R1_fq:
		#subprocess.call('cp temp/' + f + ' ' + 'temp/' + re.sub('R1','fwdp',f), shell=True)
		shutil.copy(fastq_path + '/temp/' + f, fastq_path + '/temp/' + re.sub('R1','fwdp',f))
	for f in R2_fq:
		#subprocess.call('cp temp/' + f + ' ' + 'temp/' + re.sub('R2','revp',f), shell=True)
		shutil.copy(fastq_path + '/temp/' + f, fastq_path + '/temp/' + re.sub('R2','revp',f))
elif args.adaptor_type == 'Nextera':
	for f in R1_fq:
		subprocess.call('cd ' + fastq_path + '/temp; java -jar ' + script_path + '/Trimmomatic-0.36/trimmomatic-0.36.jar PE ' + f + ' ' + re.sub('R1','R2',f) + ' ' + re.sub('R1','fwdp',f) + ' /dev/null ' + re.sub('R1','revp',f) + ' /dev/null ILLUMINACLIP:' + script_path + '/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:2:true -threads ' + args.threads, shell=True)
else:
	for f in R1_fq:
		subprocess.call('cd ' + fastq_path + '/temp; java -jar ' + script_path + '/Trimmomatic-0.36/trimmomatic-0.36.jar PE ' + f + ' ' + re.sub('R1','R2',f) + ' ' + re.sub('R1','fwdp',f) + ' /dev/null ' + re.sub('R1','revp',f) + ' /dev/null ILLUMINACLIP:' + script_path + '/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:2:true -threads ' + args.threads, shell=True)

print('AXE_ADAPTORS done!')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#FLASH
FWDP = [f for f in os.listdir(fastq_path + '/temp') if re.search('_fwdp.+.fastq',f) != None]
print(FWDP)
if args.flash == 'enabled':
	for f in FWDP:
		s = 'cd ' + fastq_path + '/temp; ' + script_path + '/FLASH-1.2.11/flash ' + f + ' ' + re.sub("fwdp","revp",f) + ' -o ' + re.sub('_fwdp.+.fastq','',f) + ' -M ' + args.max_overlap + ' -m ' + args.min_overlap
		if args.allow_outies == 'enabled':
			s = s + ' -O'
		subprocess.call(s, shell=True)
else:
	for f in FWDP:
		#subprocess.call('cp temp/' + f + ' temp/' + re.sub('_fwdp.+.fastq','.extendedFrags.fastq',f), shell=True)
		shutil.copy(fastq_path + '/temp/' + f, fastq_path + '/temp/' + re.sub('_fwdp.+.fastq','.extendedFrags.fastq',f))

print('FLASH done!')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#CREATE_TRIMMER_GENERAL
TRIMMER = [f for f in os.listdir(fastq_path + '/temp') if f.endswith('.extendedFrags.fastq')]
print(TRIMMER)
if args.trimmer == 'enabled':
	if os.path.exists(fastq_path + '/temp/ninja_shi7_report.log'):
		os.remove(fastq_path + '/temp/ninja_shi7_report.log')
	for f in TRIMMER:
		s1 = 'cd ' + fastq_path + '/temp; echo ' + f + ' >> ninja_shi7_report.log'
		s2 = script_path + '/ninja_shi7_mac ' + f + ' ' + re.sub('.extendedFrags.fastq','.trimmed.fastq',f) + ' ' + args.trim_length + ' ' + args.trim_qual + ' ' + args.trim_qual + 'FLOOR 5 ASS_QUALITY 30 >> ninja_shi7_report.log'
		subprocess.call(s1 + ' && ' + s2, shell=True)
else:
	for f in TRIMMER:
		#subprocess.call('cd temp; cp ' + f + ' ' + re.sub('.extendedFrags.fastq','.trimmed.fastq',f), shell=True)
		shutil.copy(fastq_path + '/temp/' + f, fastq_path + '/temp/' + re.sub('.extendedFrags.fastq','.trimmed.fastq',f))

print('CREATE_TRIMMER_GENERAL done!')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#MAKE_FASTAS
if os.path.exists(fastq_path + 'temp/fasta'):
	shutil.rmtree(fastq_path + 'temp/fasta')
os.mkdir(fastq_path + 'temp/fasta')
FASTA = [f for f in os.listdir(fastq_path + '/temp') if f.endswith('.trimmed.fastq')]
print(FASTA)
if args.make_fastas == 'enabled':
	iter = 0
	while iter<len(FASTA):
		fo = open(fastq_path + '/temp/' + FASTA[iter])
		record = fo.readlines()

		trim = []
		i=0
		while i<len(record):
			trim.append(record[i])
			trim.append(record[i+1])
			i=i+4

		trim_record = [re.sub('^@','>',line) for line in trim]
		fo.close()

		trim_string = ''.join(trim_record)
		fo = open(fastq_path + 'temp/fasta/'+re.sub('fastq','fasta',FASTA[iter]),'w') 
		fo.write(trim_string)
		fo.close()
		iter+=1

print('MAKE_FASTAS done!')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#QIIME to .fna
File = [f for f in os.listdir(fastq_path + 'temp/fasta') if f.endswith('.fasta')]
print(File)
if args.output != None:
		fna_path = args.output
else:
	fna_path = os.path.dirname(os.path.dirname(fastq_path)) + '/seqs.fna'

print('fna path:',fna_path)

if args.QIIME == 'enabled':
	if args.append_fna == 'enabled' and os.path.exists(fna_path): #to append
		fo = open(fna_path) #open to read the last sample index
		fo_app = open(fna_path,'a') #open the fna
		last_sample_name = fo.readlines()[-2]
		print('last sample name: ', last_sample_name)
		s = re.split(' ',last_sample_name)
		s = re.split('_',s[0])
		last_idx = int(s[1])
		print('last_idx:',last_idx)
		fo.close()
		sample_idx = last_idx + 1
	else:
		fo_app = open(fna_path,'w') #open the fna
		sample_idx = 0
		print('idx:', sample_idx)

	iter2 = 0
	fastq = []
	while iter2 < len(File):
		fo = open(fastq_path + 'temp/fasta/' + File[iter2])
		record = fo.readlines()
		record = [re.sub('^>','',line) for line in record]
		sample_name = re.sub('.trimmed.fasta','',File[iter2])
		sample_name = '>' + re.sub('[-_]','.',sample_name) + '_'
		for idx, line in enumerate(record):
			if idx%2 == 0:
				line = sample_name + str(sample_idx) + ' ' + line
				sample_idx += 1
			fastq.append(line)
		fo.close()
		iter2+=1

	fastq_string = ''.join(fastq)		
	fo_app.write(fastq_string)
	fo_app.close()

print('QIIME done!')
print('Execution time:',datetime.now() - startTime)
