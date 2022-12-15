from Bio import SeqIO

count=[]
for record in SeqIO.parse('H_influenzae_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('H_influenzae_phage_refseq.fasta',"fasta"):
   # SeqIO.write(record, "H_influenzae_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("H_influenzae_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "H_influenzae_protein_H_influenzae_{}.fasta".format(name), "fasta")
    final=[]

