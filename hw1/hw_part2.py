#!usr/bin/python

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein

Entrez.email = "A.N.Other@example.com"

# search sequences by a combination of keywords
# search sequences by a combination of keywords
handle = Entrez.esearch(db="nucleotide", term="Galdieria sulphuraria[Organism] NOT partial AND complete genome ")
records = Entrez.read(handle)
print records['Count']

top_record = records['IdList'][0]

handle = Entrez.efetch(db="nucleotide", id=top_record, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
#print record

atpase=[]
seqs=[]
locations=[]
for feature in record.features:
	if feature.type=='CDS': # feature type to be CDS
		#print "hello"
		if 'product' in feature.qualifiers: # products that contain feature qualifiers
			#print "austin sux"
			if 'ATP synthase' in feature.qualifiers['product'][0]:
				#print "lol jk"
				if str(feature.location) not in locations:
					print feature.location
					locations.append(str(feature.location))
					print locations
					atpase.append(feature)
					seqs.append(feature.qualifiers['protein_id'][0])

print len(atpase)

# rRNAs=[];

output_handle=open("G_sulphuraria_atpase_ids.fa","w")

#SeqIO.write(final, output_handle, "fasta")
for i in range(len(seqs)):
	output_handle.write(">%s %s %s\n%s\n" % (record.id,record.description,atpase[i].location,str(seqs[i])))
output_handle.close()