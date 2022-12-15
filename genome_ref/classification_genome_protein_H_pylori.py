from Bio import SeqIO

count=[]
for record in SeqIO.parse('H_pylori_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('H_pylori_phage_refseq.fasta',"fasta"):
 #   SeqIO.write(record, "H_pylori_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("H_pylori_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "H_pylori_protein_H_pylori_{}.fasta".format(name), "fasta")
    final=[]

