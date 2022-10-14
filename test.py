from Bio import SeqIO
import re
import itertools
#from collections import Counter
# file_name = "Ecoli_phage_all.fasta"    #fasta file name
# whole data number = 52261
# removed data number = 20959
# length of dic = 2942
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


    protein_revised = []
    for candidate in protein_name: # This code is used to remove symbols which disturb interpreting
        re_candidate= re.sub("[\(\[].*?[\)\]]", "", candidate) # This code remove bracket [] which contains the kinds of bacteria
        protein_revised.append(re_candidate)


    lower_protein_revised=[]   # Python distinguishes alphabet with lower case and upper case
    for lower_alphabet in protein_revised: # This code is used to change the alphabet form with lower cases
        lower_alphabet=lower_alphabet.lower()
        lower_protein_revised.append(lower_alphabet)



    removed_protein = [elements for elements in lower_protein_revised if 'hypothetical' not in elements]
    removed_protein_final = [elements for elements in removed_protein if 'putative' not in elements]
    # upper codes remove the useless words 'hypothetical and putative'


    add_capsid = [elements for elements in removed_protein_final if 'capsid' in elements]
    #print(len(add_capsid)) # checking for original number of capsid proteins
    add_capsid = list(dict.fromkeys(add_capsid))  # remove the duplicate protein in the list
    #print(len(add_capsid)) # checking for removed duplicate proteins 587->56
    #print(add_capsid)
    add_nucleic_acid = [elements for elements in removed_protein_final if 'nucle' in elements]
    #print(len(add_nucleic_acid)) # checking for original number of nucleic_acid proteins
    add_nucleic_acid = list(dict.fromkeys(add_nucleic_acid)) # remove the duplicate protein in the list
    #print(len(add_nucleic_acid)) # checking for removed duplicate proteins 2041->239
    #print(add_nucleic_acid)
    add_collar = [elements for elements in removed_protein_final if 'collar' in elements]
    #print(len(add_collar)) # checking for original number of collar proteins
    add_collar = list(dict.fromkeys(add_collar)) # remove the duplicate protein in the list
    #print(len(add_collar)) # checking for removed duplicate proteins 7->5
    #print(add_collar)
    add_whiskers = [elements for elements in removed_protein_final if 'whiskers' in elements]
    #print(len(add_whiskers)) # checking for original number of whiskers proteins
    add_whiskers = list(dict.fromkeys(add_whiskers)) # remove the duplicate protein in the list
    #print(len(add_whiskers)) # checking for removed duplicate proteins 19->6
    #print(add_whiskers)
    add_sheath = [elements for elements in removed_protein_final if 'sheath' in elements]
    #print(len(add_sheath)) # checkin for original number of sheath proteins
    add_sheath = list(dict.fromkeys(add_sheath)) # remove the duplicate protein in the list
    #print(len(add_sheath)) # checking for removed duplicate proteins 298->19
    #print(add_sheath)
    add_baseplate = [elements for elements in removed_protein_final if 'baseplate' in elements]
    #print(len(add_baseplate)) # checking for original number of basplate proteins
    add_baseplate = list(dict.fromkeys(add_baseplate)) # remove the duplicate protein in the list
    #print((len(add_baseplate))) # checking for removed duplicate proteins 1172->83
    #print(add_baseplate)
    add_tail_fiber = [elements for elements in removed_protein_final if 'tail fiber' in elements]
    #print(len(add_tail_fiber)) # checking for original number of tail fiber proteins
    add_tail_fiber = list(dict.fromkeys(add_tail_fiber)) # remove the duplicate protein in the list
    #print(len(add_tail_fiber)) # checking for removed duplicate proteins 915->109
    #print(add_tail_fiber)
    add_spike = [elements for elements in removed_protein_final if 'spike' in elements]
    #print(len(add_spike))  # checing for original number of spike proteins
    add_spike = list(dict.fromkeys(add_spike)) # remove the duplicate protein in the list
    #print(len(add_spike)) # checking for removed duplicate proteins    31->15
    #print(add_spike)

    classified_proteins = itertools.chain(add_capsid,add_nucleic_acid,add_collar,add_whiskers,add_sheath,add_baseplate,add_tail_fiber,add_spike)
    classified_proteins=list(classified_proteins) # merge the proteins which are only related to phage structures


    origin_seq=[]      # This code will retrun the lower case alphabet protein to origin alphabet
    for i in range(len(lower_protein_revised)):
        if lower_protein_revised[i] in classified_proteins:
            origin_seq.append(i)   # This code extract the sequence of proteins which were converted into lower case

    origin_protein=[]
    for number in origin_seq:
        origin_protein.append(protein_revised[int(number)]) # extract origin name of proteins


    #origin_protein = list(dict.fromkeys(origin_protein))  # remove the duplicate protein in the list
                                                           # The value was increased, This can be guessed because it did not merge lower cases alphabet and upper cases alphabet

    protein_dic={}     # make the dictionary for format like protein_name: sequence IDs (Some sequences can have the same protein)
    for key in origin_protein:   # Protein name
        for value in Number_of_ID:   # Sequence ID
            if key in protein_dic:   # If there are not protein name in the dictionary, add the new protein name
                protein_dic[key]+=[value]
            else:
                protein_dic[key]=[value]  # If there are protein name in the dictionary, add the new Sequence ID in the value part
            Number_of_ID.remove(value)

            break    # stop the loop if there are no Sequence ID in the Protein name list
    return protein_dic
    '''
    This Part is used to distinguish unique words, but it is not yet complete code
    
    #print(sorted(protein_dic.keys()))
    #words= sorted(protein_dic.keys())
    #print(words)
    #words_list=[]
    #for word in words:
    #    words_repeat=word.replace('_',' ')
    #    words_repeat= words_repeat.replace('-',' ').split(' ')
    #    words_repeat = [x for x in words_repeat if x != '' and  x !='protein']
     #   words_list.append(words_repeat)

    #print(words_list)
    #count = 0
    #counting_list = []
    #checking_list=[]
    #for compare in words:
    #    for repeated_example in words_list:
     #       for compare2 in repeated_example:
     #           if compare2 in compare:
      #              count+=1
      #  counting_list.append(count)
      #  if count<300:
       #     checking_list.append(compare)
       # count=0

    #print(checking_list)



    #unique =[]
    #count=0
    #counting_list=[]
    #print(len(words))
    #for i in range(len(words)):
     #   for words_counting in words:
      #      if str(words[i]) in words_counting:
       #         count+=1
       # unique.append(count)
       # if count >1:
        #    counting_list.append(i)
        #count=0
    #print(len(words))
    #print(counting_list)
    #words_remove= sorted(protein_dic.keys())
    #for i in range(len(words)):
     #   for j in range(len(counting_list)):
      #      if i== int(counting_list[j]):
       #         words_remove.remove(words[i])
    #print(len(words))
    #print(len(words_remove))
    #print(words_remove)




    # print(sorted(protein_dic.keys()))

    #print(word)
    #print(len(unique))
    #return protein_dic
    #filt_keys=['riia lysis inhibitor']  # This part is example for validating the dictionary
    #return protein_dic
    #search=[protein_dic[key] for key in filt_keys]     # It shows the number of sequence ID for protein name example
    #for key, value in sorted(protein_dic.items()):    # It shows the number of values in the keys
     #   print(key, len(value))
    #print(sorted(protein_dic.keys()))

#def bracket_remove(file_name):
 #   protein_name_candidate = []
  #  for seq_record in SeqIO.parse(file_name, "fasta"):  # open the fasta file
   #     all_species.append(seq_record.description.split('|'))
    #for protein in all_species:
     #   protein_name_candidate.append(protein[2])
'''

