
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein

Entrez.email = "A.N.Other@example.com"
top_records = []

#search nematodes
handle = Entrez.esearch(db="nucleotide",term="nematode[Organism] AND 28S rRNA[Gene] AND 701:10000000[Sequence Length]")
records = Entrez.read(handle)
print records['Count']

#search greeen algae
handle = Entrez.esearch(db="nucleotide",retmax=304,term="Chlorophyta[Organism] AND 28S rRNA[Gene] AND 701:10000000[Sequence Length] ")
records = Entrez.read(handle)
print records['Count']

for i in range(1,304):
	top_records.append(records['IdList'][i])

#search Ascomycete fungi
handle = Entrez.esearch(db="nucleotide",term="ascomycete AND 28S rRNA[Gene] AND 701:10000000[Sequence Length]")
records = Entrez.read(handle)
print records['Count']

# retrieve the sequence by their GI numbers
gi_list = ','.join(top_records)
#it works so i dont need to see this anymore
#print "These are the GI Numbers chosen: ",gi_list
handle = Entrez.efetch(db="nucleotide",id=gi_list,rettype="gb",retmode="xml")
my_genbank_records = Entrez.read(handle)
handle.close()

# convert the Genbank to FASTA
my_fasta_records = []
for i in range(len(my_genbank_records)):
	my_fasta_records.append(SeqRecord(Seq(my_genbank_records[i]['GBSeq_sequence']),id=my_genbank_records[i]['GBSeq_primary-accession'],description=my_genbank_records[i]['GBSeq_definition']))
	
# output file 
out_file = open("long_28rrna_greenalgae.fa","w")
SeqIO.write(my_fasta_records,out_file,"fasta")
out_file.close()

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
	if feature.type=='CDS':
		#print "hello"
		if 'product' in feature.qualifiers:
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

