from Bio import SeqIO

count=[]
for record in SeqIO.parse('E_faecium_phage_refseq.fasta',"fasta"):
    count.append(record.id)

#for record in SeqIO.parse('E_faecium_phage_refseq.fasta',"fasta"):
 #   SeqIO.write(record, "E_faecium_genomic_{}.fasta".format(record.id), "fasta")

example = [seq_record for seq_record in SeqIO.parse("E_faecium_refseq_proteins.fasta", "fasta")]


final=[]
for name in count:
    for record in example:
        if name in record.id:
            final.append(record)
    SeqIO.write(final, "E_faecium_protein_E_faecium_{}.fasta".format(name), "fasta")
    final=[]

