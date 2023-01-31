from Bio import SeqIO


sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta') # this code is the example of A_baumannii

# below codes separates the fasta files by the sequence length
# the range of length are 1 to 500, if there are sequences length over 500, combine the whole example
len_1_20_text = [seq_len for seq_len in sequences if 1<=len(seq_len.seq)<=20]

with open('A_baumannii 1 to 20 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_1_20_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_21_100_text = [seq_len for seq_len in sequences if 21<=len(seq_len.seq)<=100]

with open('A_baumannii 21 to 100 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_21_100_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_101_150_text = [seq_len for seq_len in sequences if 101<=len(seq_len.seq)<=150]

with open('A_baumannii 101 to 150 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_101_150_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_151_200_text = [seq_len for seq_len in sequences if 151<=len(seq_len.seq)<=200]

with open('A_baumannii 151 to 200 inhibitor.fasta', 'wt') as output:
   SeqIO.write(len_151_200_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_201_250_text = [seq_len for seq_len in sequences if 201<=len(seq_len.seq)<=250]

with open('A_baumannii 201 to 250 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_201_250_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_251_300_text = [seq_len for seq_len in sequences if 251<=len(seq_len.seq)<=300]

with open('A_baumannii 251 to 300 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_251_300_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_301_350_text = [seq_len for seq_len in sequences if 301<=len(seq_len.seq)<=350]

with open('A_baumannii 301 to 350 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_301_350_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_351_400_text = [seq_len for seq_len in sequences if 351<=len(seq_len.seq)<=400]

with open('A_baumannii 351 to 400 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_351_400_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_401_450_text = [seq_len for seq_len in sequences if 401<=len(seq_len.seq)<=450]

with open('A_baumannii 401 to 450 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_401_450_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_451_500_text = [seq_len for seq_len in sequences if 451<=len(seq_len.seq)<=500]

with open('A_baumannii 451 to 500 inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_451_500_text, output, 'fasta')

sequences = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')
len_500_over_text = [seq_len for seq_len in sequences if len(seq_len.seq)>500]

with open('A_baumannii 500_over inhibitor.fasta', 'wt') as output:
    SeqIO.write(len_500_over_text, output, 'fasta')

sequences_check = SeqIO.parse('A_baumannii_inhibitor.fasta', 'fasta')










