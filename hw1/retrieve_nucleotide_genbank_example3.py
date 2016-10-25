'''

Retrieve Genbank entries from the nucleotide database at NCBI.

-----------------------------------------------------------
(c) 2013 Allegra Via and Kristian Rother
    Licensed under the conditions of the Python License

    This code appears in section 20.4.3 of the book
    "Managing Biological Data with Python".
-----------------------------------------------------------
'''

from Bio import Entrez
from Bio import SeqIO
from Bio import Seq	
from Bio.Alphabet import IUPAC

Entrez.email = "Aclass@drexel.edu" 
# search sequences by a combination of keywords
handle = Entrez.esearch(db="nucleotide", term="16S rRNA[gene] AND streptococcus[ORGN] NOT partial AND complete genome")
records = Entrez.read(handle)
print records['Count']

top3_records = records['IdList'][0:3]

handle = Entrez.efetch(db="nucleotide", id=records['IdList'][0], rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

handle = Entrez.efetch(db="nucleotide", id="FN568063", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

sixteen_s=[]
seqs=[]
locations=[]
for feature in record.features:
	if feature.type=='gene' or feature.type == 'rRNA':
		if 'gene' in feature.qualifiers:
			if feature.qualifiers['gene'][0]=='16S rRNA':
				#sixteen_s.append(seq)
				#sixteen_s.append(feature)

				if str(feature.location) not in locations:
					print feature.location
					locations.append(str(feature.location))
					print locations
					sixteen_s.append(feature)
					seqs.append(feature.extract(record.seq))

print len(sixteen_s)



# rRNAs=[];

output_handle=open("rRNAs.fa","w")

#SeqIO.write(final, output_handle, "fasta")
for i in range(len(seqs)):
	output_handle.write(">%s %s %s\n%s\n" % (record.id,record.description,sixteen_s[i].location,str(seqs[i])))
output_handle.close()