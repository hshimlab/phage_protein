from Bio import SeqIO
# file_name = "Ecoli_phage_all.fasta"    #fasta file name
def fasta_seq(file_name):
    all_species=[]    # remove the '|' in the description part
    protein_name=[]   # make the list for classifying protein Name
    Number_of_ID=[]   # make the list for classifying sequence ID
    for seq_record in SeqIO.parse(file_name, "fasta"):   # open the fasta file
        all_species.append(seq_record.description.split('|')) # remove '|' in the description part
        Number_of_ID.append(seq_record.id) # add the sequence IDs to the list

    for protein in all_species:
        protein_name.append(protein[2])     # extract the protein name from the description part

    protein_name_list= list(dict.fromkeys(protein_name))   # this code is used for checking the number of duplicate protein names

    protein_dic={}     # make the dictionary for format like protein_name: sequence IDs (Some sequences can have the same protein)

    for key in protein_name:   # Protein name
        for value in Number_of_ID:   # Sequence ID
            if key in protein_dic:   # If there are not protein name in the dictionary, add the new protein name
                protein_dic[key]+=[value]
            else:
                protein_dic[key]=[value]  # If there are protein name in the dictionary, add the new Sequence ID in the value part
            Number_of_ID.remove(value)

            break    # stop the loop if there are no Sequence ID in the Protein name list

    #filt_keys=['host cell division inhibitor Icd-like protein [Enterobacteria phage P4]']  # This part is example for validating the dictionary

    #search=[protein_dic[key] for key in filt_keys]     # It shows the number of sequence ID for protein name example

    #for key, value in sorted(protein_dic.items()):    # It shows the number of values in the keys
        #print(key, len(value))

    return protein_dic

