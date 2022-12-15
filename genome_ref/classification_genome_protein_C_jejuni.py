from Bio import SeqIO

count=[]
for record in SeqIO.parse('C_jejuni_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('C_jejuni_phage_refseq.fasta',"fasta"):
 #   SeqIO.write(record, "C_jejuni_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("C_jejuni_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "C_jejuni_protein_C_jejuni_{}.fasta".format(name), "fasta")
    final=[]

