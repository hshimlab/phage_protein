from Bio import SeqIO

count=[]
for record in SeqIO.parse('E_coli_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('E_coli_phage_refseq.fasta',"fasta"):
  #  SeqIO.write(record, "E_coli_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("E_coli_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "E_coli_protein_E_coli_{}.fasta".format(name), "fasta")
    final=[]

