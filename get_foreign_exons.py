from utility_scripts_for_get_foreign_exons import *
from collections import Counter
import re 


infile=open('./MLOF_final_drb_dpb_dqb.csv')
header=infile.readline().strip()+',New_type,length_of_foreign,foreign_seq,foreign_perc_rank_1,foreign_perc_rank_2\n'
lines=infile.readlines()
lines=[i.strip().split(',') for i in lines]# if i.strip().split(',')[0]=='WH3YQM8ZP']
f8_aa=open('./f8_coding.fa').readlines()[1]
f8_seq=nuc_to_protein(f8_aa)

#for l in lines[:10]:
#	print(l[:5])
#	print (','.join(l[40:46]))
ct=0
new_lines=[]
for line in lines:
	i=line[:]
	#print(i)
	if i[1] in ['No Record of Subject','Not enough DNA','Excluded','Clinical Only','Clinical Only - Discarded']:
		continue
		
	
	if i[45]=='c.6429+?_6430-?inv':
		i+=get_results_i22_inversion(i[45],'I22-Inversion',i[1])
	elif i[45]=='c.143+?_144-?inv' or i[45]=='c.143+?_144-?Inv' :
		i+=get_results_i1_inversion(i[45],'I1-Inversion',i[1])
	
	elif i[42]=='Frameshift' and i[43]=='Deletion':
		i+=get_results_frameshift_deletion(i[45],'Frameshift-Deletion',i[1])
	elif i[42]=='Frameshift' and i[43]=='Duplication':
		i+=get_results_frameshift_duplication(i[45],'Frameshift-Duplication',i[1])
	elif i[42]=='Frameshift' and i[43]=='Insertion/Deletion':
		i+=get_results_frameshift_indel(i[45],'Frameshift-Insertion-Deletion',i[1])
	elif i[42]=='Frameshift' and i[43]=='Substitution':
		i+=get_results_missense(i[45],'Frameshift-Substituion',i[1])	
	elif i[42]=='Frameshift' and i[43]=='Insertion':
		i+=get_results_frameshift_insertion(i[45],'Frameshift-Insertion',i[1])	
		
	elif i[42]=='Missense':
			i+=get_results_missense(i[45],'Missense',i[1])
		
	elif i[42]=='Nonsense':
		i+=get_results_nonsense(i[45],'Nonsense',i[1])
			
	
	elif i[42].startswith('Small') and i[43]=='Deletion':
		i+=get_results_frameshift_deletion(i[45],'Small-Deletion',i[1])
	elif i[42].startswith('Small') and i[43]=='Insertion/Deletion':
		i+=get_results_frameshift_indel(i[45],'Small-Indel',i[1])
		
	elif i[42].startswith('Large') and i[43]=='Deletion':
		i+=get_results_large_deletion(i[45],'Large-Deletion',i[1])
	
	elif i[42]=='Synonymous':
		i+=get_results_synonymous(i[1],'Synonymous')
	elif i[42]=='NA' and i[43]=='NA' and i[45]=='NA' and i[44]=='NA':
		i+=get_results_synonymous(i[1],'No Clinical Mutation')
	elif i[42]=='Upstream':
		i+=get_results_synonymous(i[1],'Upstream Mutation')
	
	else:
		i+=just_check_other_muts(i[1])
	
	new_lines.append(i)
	
all_for_15mers=[]
all_alleles=[]
for i in new_lines:
	#print((i))
	#print(i[-6])
	for p in i[-2]:
		if p not in all_for_15mers: all_for_15mers.append(p)
	for j in [i[46],i[47]]:
		if j not in all_alleles: all_alleles.append(j)

'''
outfile=open('all_for_15mers.txt','w')	
for i in all_for_15mers:
	outfile.write('%s\n'%i)

outfile.close()

outfile=open('all_alleles.txt','w')	
for i in all_alleles:
	outfile.write('DRB1_%s%s\n'%(i[:2],i[-2:]))

outfile.close()
'''

ct=0
outfile=open('./MLOF_final_drb_dpb_dqb_and_foreign.csv','w')
outfile.write(header)
for i in new_lines:
	print('%d of %d'%(ct,len(new_lines)))
	ct+=1
	i=i[:-1]
	hla1='DRB1_%s%s'%(i[46][:2],i[46][3:])
	hla2='DRB1_%s%s'%(i[47][:2],i[47][3:])
	pct_rank_1=[]
	pct_rank_2=[]
	
	if(len(i[50]))>0:
		for j in i[50]:
			pct_rank_1.append(get_affinity(j,hla1))
			pct_rank_2.append(get_affinity(j,hla2))
	
	i[50]=';'.join(i[50])
	i.append(';'.join(pct_rank_1))
	i.append(';'.join(pct_rank_2))
	
	outfile.write(','.join(map(str,i))+'\n')
outfile.close()
