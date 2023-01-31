from Bio import SeqIO

count=[]
for record in SeqIO.parse('A_baumannii_phage_refseq.fasta',"fasta"):
    count.append(record.id) # collect the protein ID which are separated by the start point and end point

for record in SeqIO.parse('A_baumannii_phage_refseq.fasta',"fasta"):
    SeqIO.write(record, "A_baumannii_genomic_A_baumannii_{}.fasta".format(record.id), "fasta")
# above code can separate one sequence file to multiple protein names file. Fasta files are not translated into aminoacids they are wrote as nucleic acids

# open the protein fasta files which were translated into aminoacids
example = [seq_record for seq_record in SeqIO.parse("A_baumannii_refseq_proteins.fasta", "fasta")]

# separate each protein sequences according to start point and end point. Files are wrote into aminoacids
final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "A_baumannii_protein_A_baumannii_{}.fasta".format(name), "fasta")
    final=[]

