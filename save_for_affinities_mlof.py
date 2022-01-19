import pandas as pd
import pymysql
import string
import subprocess
import time
from multiprocessing import Pool
import multiprocessing
from tqdm import tqdm
import sys

#mysql database connection

rand_str = lambda n: ''.join([random.choice(string.ascii_lowercase) for i in range(n)])
	
alleles=[i.strip() for i in open('all_alleles.txt').readlines()]
print(alleles)


def get_promiscuity(peptide,alleles=alleles):
	db2=pymysql.connect(host='127.0.0.1',user='jmcgill',passwd='Mysqlpass1984!',db='netmhciipan31')
	sql_comm="select promiscuity from deimmunizecas9_donors where (peptide = \'%s\') and (threshold = %d);"%(peptide,threshold)
	cur=db2.cursor()
	cur.execute(sql_comm)
	results=cur.fetchall()
	if len(results)==0:
		s = rand_str(10)
		outfile=open('./tmp/'+s+'.fasta','w')
		outfile.write('>temp\n')
		outfile.write(peptide)
		outfile.close()
		bash_command = 'for allele in ' + ' '.join(list(alleles.index)) + '; do netmhciipan -f ./tmp/'+ s+'.fasta -a $allele -length '+str(len(peptide))+'; done'
		process = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE)
		process_results=str(process.communicate()).split('\\n')
		#print(process_results[:5])
		# terminal output of netmhc run as a list
		line_list = list([1 if float(i.split()[9])<threshold else 0 for i in process_results if len(i.split()) in [11,12] and i.split()[0]=='1'])
		#line_list = [i+['NA']  if len(i)==11 else i for i in line_list ]
		#print(line_list)
		prom=sum([float(line_list[i])*((alleles[i])) for i in range(len(line_list))])

		sql_comm="insert into deimmunizecas9b (peptide,changes,blosum_scores,promiscuity,threshold) values (\'%s\',\'%s\',\'%s\',%f,%f);"%(peptide,changes,blosum_scores,prom,threshold)		
		cur=db2.cursor()
		cur.execute(sql_comm)
		db2.commit()
		try:
			subprocess.call('rm ./tmp/%s.fasta'%(s))
		except:
			pass
		return prom
	else:
		return list(results)[0][0]
	db2.close()



'''
for pool in range(1):
	for number_of_changes in [10]:
		for run in range(25):
			print('%s Pool: %s Run: %d Number: %d Last Run: %.2fm'%(time.asctime(),pool,run,number_of_changes,time_to_print/60))
			start_total=time.time()
			simulated_annealing(number_of_changes,20,starting_sequence,pool,number_of_iterations=10000,temperature_decrease_alpha=.001)
			end_total=time.time()
			time_to_print=end_total-start_total
'''