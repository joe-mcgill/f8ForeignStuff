
import re

def nuc_to_protein(nuc_string):
	table = { 
	        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
	        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
	        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
	        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
	        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
	        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
	        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
	        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
	        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
	        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
	        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
	        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
	        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
	        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
	        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
	        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	    } 
	protein ="" 
	while len(nuc_string)%3 != 0:
		nuc_string=nuc_string[:-1]
	for i in range(0, len(nuc_string), 3): 
		codon = nuc_string[i:i + 3] 
		protein+= table[codon] 
	return protein 


f8_aa=open('./f8_coding.fa').readlines()[1]
f8_seq=nuc_to_protein(f8_aa)

other_mutations=open('./VariantCalls_5-2-2019.csv').readlines()
other_mutations=[i for i in other_mutations if i.split(',')[3]=='F8' and i.split(',')[4][2]!='*' and re.match('p\.([a-zA-Z]{3})([0-9]*)([a-zA-Z]*)',i.split(',')[5])]
lst=['\t'.join(j) for j in sorted([i.split(',') for i in other_mutations],key=lambda x: x[0])]

all_affinities=open('./affinities_15mers.csv').readlines()
all_affinities=[i.strip()[1:-1].split('","') for i in all_affinities]

def get_affinity(peptide,allele):
	print(peptide)
	print(allele)
	temp_affinities=[i for i in all_affinities if i[0]==peptide and i[2]==allele]
	print(temp_affinities)
	return(temp_affinities[0][4])


def add_non_clinical_mutaion(for_aa,pat_id):
	#print(pat_id)
	pat_mutation=[i for i in other_mutations if i.split(',')[0]==pat_id]
	#print(pat_mutation)
	temp_muts=[]
	for i in pat_mutation:
		#print(i)
		pattern=re.compile('^p.([a-zA-Z]*)([0-9]*)([a-zA-Z]*)')
		mtchs=re.search(pattern,i.split(',')[5])
		if mtchs is not None:

			rslts=mtchs.group(1,2,3)
			#print('%s\t%s  \t %s'%(rslts[0],rslts[1],rslts[2]))
			temp_muts+=[f8_seq[(int(rslts[1])+x):(int(rslts[1])+15+x)] for x in range(-15,0,1)]
			#sprint(temp_muts)
	return for_aa+temp_muts
	
len_f8=len(f8_seq)




def turn_f8_string_to_foreign_15mers(pat_f8seq,base_seq=f8_seq):
	for_15mers=[]
	for i in range(len_f8-14):
		startpos=pat_f8seq.find(f8_seq[i:(i+15)])
		if startpos<0:
			#print(f8_seq[i:(i+15)])
			for_15mers.append(f8_seq[i:(i+15)])
	return for_15mers

def get_results_i22_inversion(input_mutation,description,pat_id):
	seq='KWQTYRGNSTGTLMVFFGNVDSSGIKHN'
	for_aa=[seq[i:i+15] for i in range(15)]
	#print(for_aa)
	or_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )


def get_results_i1_inversion(input_mutation,description,pat_id):
	seq='RFPPRVPKSFPFNTSVVYKKTLFVEFTDHLFNIAKPRPPWMGLLGPTIQAEVYDTVVITLKNMASHPVSLHAVGVSYWKASEGAEYDDQTSQREKEDDKVFPGGSHTYVWQVLKENGPMASDPLCLTYSYLSHVDLVKDLNSGLIGALLVCREGSLAKEKTQTLHKFILLFAVFDEGKSWHSETKNSLMQDRDAASARAWPKMHTVNGYVNRSLPGLIGCHRKSVYWHVIGMGTTPEVHSIFLEGHTFLVRNHRQASLEISPITFLTAQTLLMDLGQFLLFCHISSHQHDGMEAYVKVDSCPEEPQLRMKNNEEAEDYDDDLTDSEMDVVRFDDDNSPSFIQIRSVAKKHPKTWVHYIAAEEEDWDYAPLVLAPDDRSYKSQYLNNGPQRIGRKYKKVRFMAYTDETFKTREAIQHESGILGPLLYGEVGDTLLIIFKNQASRPYNIYPHGITDVRPLYSRRLPKGVKHLKDFPILPGEIFKYKWTVTVEDGPTKSDPRCLTRYYSSFVNMERDLASGLIGPLLICYKESVDQRGNQIMSDKRNVILFSVFDENRSWYLTENIQRFLPNPAGVQLEDPEFQASNIMHSINGYVFDSLQLSVCLHEVAYWYILSIGAQTDFLSVFFSGYTFKHKMVYEDTLTLFPFSGETVFMSMENPGLWILGCHNSDFRNRGMTALLKVSSCDKNTGDYYEDSYEDISAYLLSKNNAIEPRSFSQNSRHPSTRQKQFNATTIPENDIEKTDPWFAHRTPMPKIQNVSSSDLLMLLRQSPTPHGLSLSDLQEAKYETFSDDPSPGAIDSNNSLSEMTHFRPQLHHSGDMVFTPESGLQLRLNEKLGTTAATELKKLDFKVSSTSNNLISTIPSDNLAAGTDNTSSLGPPSMPVHYDSQLDTTLFGKKSSPLTESGGPLSLSEENNDSKLLESGLMNSQESSWGKNVSSTESGRLFKGKRAHGPALLTKDNALFKVSISLLKTNKTSNNSATNRKTHIDGPSLLIENSPSVWQNILESDTEFKKVTPLIHDRMLMDKNATALRLNHMSNKTTSSKNMEMVQQKKEGPIPPDAQNPDMSFFKMLFLPESARWIQRTHGKNSLNSGQGPSPKQLVSLGPEKSVEGQNFLSEKNKVVVGKGEFTKDVGLKEMVFPSSRNLFLTNLDNLHENNTHNQEKKIQEEIEKKETLIQENVVLPQIHTVTGTKNFMKNLFLLSTRQNVEGSYDGAYAPVLQDFRSLNDSTNRTKKHTAHFSKKGEEENLEGLGNQTKQIVEKYACTTRISPNTSQQNFVTQRSKRALKQFRLPLEETELEKRIIVDDTSTQWSKNMKHLTPSTLTQIDYNEKEKGAITQSPLSDCLTRSHSIPQANRSPLPIAKVSSFPSIRPIYLTRVLFQDNSSHLPAASYRKKDSGVQESSHFLQGAKKNNLSLAILTLEMTGDQREVGSLGTSATNSVTYKKVENTVLPKPDLPKTSGKVELLPKVHIYQKDLFPTETSNGSPGHLDLVEGSLLQGTEGAIKWNEANRPGKVPFLRVATESSAKTPSKLLDPLAWDNHYGTQIPKEEWKSQEKSPEKTAFKKKDTILSLNACESNHAIAAINEGQNKPEIEVTWAKQGRTERLCSQNPPVLKRHQREITRTTLQSDQEEIDYDDTISVEMKKEDFDIYDEDENQSPRSFQKKTRHYFIAAVERLWDYGMSSSPHVLRNRAQSGSVPQFKKVVFQEFTDGSFTQPLYRGELNEHLGLLGPYIRAEVEDNIMVTFRNQASRPYSFYSSLISYEEDQRQGAEPRKNFVKPNETKTYFWKVQHHMAPTKDEFDCKAWAYFSDVDLEKDVHSGLIGPLLVCHTNTLNPAHGRQVTVQEFALFFTIFDETKSWYFTENMERNCRAPCNIQMEDPTFKENYRFHAINGYIMDTLPGLVMAQDQRIRWYLLSMGSNENIHSIHFSGHVFTVRKKEEYKMALYNLYPGVFETVEMLPSKAGIWRVECLIGEHLHAGMSTLFLVYSNKCQTPLGMASGHIRDFQITASGQYGQWAPKLARLHYSGSINAWSTKEPFSWIKVDLLAPMIIHGIKTQGARQKFSSLYISQFIIMYSLDGKKWQTYRGNSTGTLMVFFGNVDSSGIKHNIFNPPIIARYIRLHPTHYSIRSTLRMELMGCDLNSCSMPLGMESKAISDAQITASSYFTNMFATWSPSKARLHLQGRSNAWRPQVNNPKEWLQVDFQKTMKVTGVTTQGVKSLLTSMYVKEFLISSSQDGHQWTLFFQNGKVKVFQGNQDSFTPVVNSLDPPLLTRYLRIHPQSWVHQIALRMEVLGCEAQDLY'
	for_aa=[seq[i:i+15] for i in range(15)]
	#print(for_aa)
	or_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )

def get_results_synonymous(pat_id,desc):
	for_aa=[]
	#print(for_aa)
	or_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [desc,len(for_aa),for_aa,'placeholder for prom'] )
def get_results_frameshift_deletion(input_mutation,description,pat_id):
	pattern=re.compile('c\.[-]*([0-9]+)[-_\?\+]*([0-9]*)[-_\?\+]*del([ACGT]*)')
	match=re.search(pattern,input_mutation)
	re_results=match.group(1,2,3)
	if len(re_results[1])==0:
		new_seq=f8_aa[:int(re_results[0])-1]+f8_aa[int(re_results[0]):]
	else:
		new_seq=f8_aa[:int(re_results[0])-1]+f8_aa[int(re_results[1]):]
	while new_seq[:3]!='ATG':
		new_seq=new_seq[1:]
	#if new_seq[:3]!='ATG':
		#print(new_seq[:3])
	
	new_aa=nuc_to_protein(new_seq).split('_')[0]
	for_aa=turn_f8_string_to_foreign_15mers(new_aa)
	for_aa=add_non_clinical_mutaion(for_aa,pat_id)
	
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )

def get_results_frameshift_duplication(input_mutation,description,pat_id):
	pattern=re.compile('c\.([0-9]+)[_]*([0-9]*)dup([ACGT]*)')
	match=re.search(pattern,input_mutation)
	re_results=match.group(1,2,3)
	if len(re_results[2])==0:
		new_seq=f8_aa[:(int(re_results[0])-1)] + f8_aa[(int(re_results[0])-1):int(re_results[1])] +f8_aa[(int(re_results[0])-1):]
	else:
		new_seq=f8_aa[:int(re_results[0])] + re_results[2] +f8_aa[int(re_results[0]):]
	new_aa=nuc_to_protein(new_seq).split('_')[0]
	for_aa=turn_f8_string_to_foreign_15mers(new_aa)
	for_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )

def get_results_frameshift_insertion(input_mutation,description,pat_id):
	pattern=re.compile('c\.([0-9]+)[_]*([0-9]*)ins([ACGT]*)')
	match=re.search(pattern,input_mutation)
	re_results=match.group(1,2,3)
	if len(re_results[2])==1:
		new_seq=f8_aa[:(int(re_results[0]))] + re_results[2] +f8_aa[(int(re_results[0])):]
	else:
		new_seq=f8_aa[:int(re_results[0])] + re_results[2] +f8_aa[int(re_results[0]):]
	
	new_aa=nuc_to_protein(new_seq).split('_')[0]
	for_aa=turn_f8_string_to_foreign_15mers(new_aa)
	for_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )


def get_results_missense(input_mutation,description,pat_id):
	pattern=re.compile('c\.([0-9]+)([ACGT])>([ACGT])')
	match=re.search(pattern,input_mutation)
	re_results=match.group(1,2,3)
	new_seq=f8_aa[:(int(re_results[0])-1)] + re_results[2] +f8_aa[(int(re_results[0])):]
	new_aa=nuc_to_protein(new_seq)
	for_aa=turn_f8_string_to_foreign_15mers(new_aa)
	#print(for_aa)
	for_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )

def get_results_nonsense(input_mutation,description,pat_id):
	pattern=re.compile('c\.([0-9]+)([ACGT])>([ACGT])')
	match=re.search(pattern,input_mutation)
	re_results=match.group(1,2,3)
	#print(re_results[0])
	new_seq=f8_aa[:int(re_results[0])]
	new_aa=nuc_to_protein(new_seq)
	for_aa=turn_f8_string_to_foreign_15mers(new_aa)
	for_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )	

def get_results_frameshift_indel(input_mutation,description,pat_id):
	pattern=re.compile('c\.([0-9]+)[_]*([0-9]*)del[ACGT]*ins([ACGT]+)')
	match=re.search(pattern,input_mutation)
	re_results=match.group(1,2,3)
	
	if len(re_results[1])==0:
		new_seq=f8_aa[:int(re_results[0])-1]+re_results[2]+f8_aa[int(re_results[0]):]
	else:
		new_seq=f8_aa[:int(re_results[0])-1]+re_results[2]+f8_aa[int(re_results[1]):]
	new_aa=nuc_to_protein(new_seq).split('_')[0]
	for_aa=turn_f8_string_to_foreign_15mers(new_aa)
	for_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )

def get_results_large_deletion(input_mutation,description,pat_id):
	#c.1-?_2113+?del
	#print(input_mutation)

	pattern=re.compile('c\.[-]*([0-9]+)[_]*[+-]*\?*_*[c.]*([0-9]+).+')
	match=re.search(pattern,input_mutation)
	re_results=match.group(1,2)
	new_seq=f8_aa[:int(re_results[0])-1]+f8_aa[int(re_results[1]):]
	new_aa=nuc_to_protein(new_seq).split('_')[0]
	#print(new_aa)
	for_aa=turn_f8_string_to_foreign_15mers(new_aa)
	for_aa=add_non_clinical_mutaion(for_aa,pat_id)
	#print(for_aa)
	return( [description,len(for_aa),for_aa,'placeholder for prom'] )

def just_check_other_muts(pat_id):
	for_aa=add_non_clinical_mutaion([],pat_id)
	return( ['Ambiguous Clinical Mutations',len(for_aa),for_aa,'placeholder for prom'] )
#def get_affinities_foreign(mer_list):
	
#print turn_f8_string_to_foreign_15mers(test_pat_f8)
