from Bio import SeqIO

count=[]
for record in SeqIO.parse('A_baumannii_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('A_baumannii_phage_refseq.fasta',"fasta"):
   # SeqIO.write(record, "A_baumannii_genomic_A_baumannii_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("A_baumannii_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "A_baumannii_protein_A_baumannii_{}.fasta".format(name), "fasta")
    final=[]

