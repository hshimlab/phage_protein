from Bio import SeqIO

count=[]
for record in SeqIO.parse('Shigella_flexneri_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('Shigella_flexneri_phage_refseq.fasta',"fasta"):
 #   SeqIO.write(record, "Shigella_flexneri_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("Shigella_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "Shigella_protein_{}.fasta".format(name), "fasta")
    final=[]

