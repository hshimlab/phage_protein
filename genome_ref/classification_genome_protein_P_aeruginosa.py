from Bio import SeqIO

count=[]
for record in SeqIO.parse('P_aeruginosa_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('P_aeruginosa_phage_refseq.fasta',"fasta"):
 #   SeqIO.write(record, "P_aeruginosa_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("P_aeruginosa_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "P_aeruginosa_protein_P_aeruginosa_{}.fasta".format(name), "fasta")
    final=[]

