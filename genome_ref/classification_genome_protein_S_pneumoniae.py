from Bio import SeqIO

count=[]
for record in SeqIO.parse('S_pneumoniae_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('S_pneumoniae_phage_refseq.fasta',"fasta"):
 #   SeqIO.write(record, "S_pneumoniae_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("S_pneumoniae_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "S_pneumoniae_protein_S_pneumoniae_{}.fasta".format(name), "fasta")
    final=[]

