#!/usr/bin/python

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein

Entrez.email = "A.N.Other@example.com"
top_records = []

#search sequence 
handle = Entrez.esearch(db="nucleotide",term="nematode[Organism] AND 28S rRNA[Gene] AND 701:10000000[Sequence Length]")
records = Entrez.read(handle)
print records['Count']
#top_records.append(records['IdList'][0])


handle = Entrez.esearch(db="nucleotide",retmax=304,term="Chlorophyta[Organism] AND 28S rRNA[Gene] AND 701:10000000[Sequence Length] ")
records = Entrez.read(handle)
print records['Count']

for i in range(1,304):
	top_records.append(records['IdList'][i])
#print records

handle = Entrez.esearch(db="nucleotide",term="ascomycete AND 28S rRNA[Gene] AND 701:10000000[Sequence Length]")
records = Entrez.read(handle)
print records['Count']







#top_records.append(records['IdList'][0])


# retrieve the sequence by their GI numbers
gi_list = ','.join(top_records)
print "These are the GI Numbers chosen: ",gi_list
handle = Entrez.efetch(db="nucleotide",id=gi_list,rettype="gb",retmode="xml")
my_genbank_records = Entrez.read(handle)
handle.close()


# convert the Genbank to FASTA
my_fasta_records = []
for i in range(len(my_genbank_records)):
	my_fasta_records.append(SeqRecord(Seq(my_genbank_records[i]['GBSeq_sequence']),id=my_genbank_records[i]['GBSeq_primary-accession'],description=my_genbank_records[i]['GBSeq_definition']))
	
# output file 
one_file = open("long_28rrna_greenalgae.fa","w")
SeqIO.write(my_fasta_records,one_file,"fasta")
one_file.close()

