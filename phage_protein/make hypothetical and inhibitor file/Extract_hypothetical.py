from Bio import SeqIO

sequences = SeqIO.parse('E_faecium_phage_protein.fasta', 'fasta') # This code is the example of E_faecium
filtered = [seq for seq in sequences if 'yphtothetical' in seq.description.split('|')[1]]

with open('E_faecium_revised.fasta', 'wt') as output:
    SeqIO.write(filtered, output, 'fasta')
