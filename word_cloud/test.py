from Bio import SeqIO
import re
import itertools
import matplotlib.pyplot as plt
from wordcloud import WordCloud
#from collections import Counter
# file_name = "Ecoli_phage_all.fasta"    #fasta file name
# whole data number = 52261
# removed data number = 20959
# length of dic = 2941
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

# make the list of whole proteins.
########################################################################################################################
# extract all kind of proteins

    list_extend_all=[] # make the list for all kind of proteins
    for words in lower_protein_revised:
        words_repeat= words.split(' ')
        words_repeat = [x for x in words_repeat if x != ''] # many proteins contains ''(not meaningful), and remove them
        list_extend_all.extend(words_repeat) # merge the multiple lists as one list



    protein_dic_all_words = {}  # To count the duplicate words, make the dictionary
    for key in list_extend_all:
       if key in protein_dic_all_words:
           protein_dic_all_words[key]+=1  # If duplicated proteins are, add +1 for counting
       else:
           protein_dic_all_words[key]=1   # If there are not duplicated proteins, remain as 1

    sorted_proteins_all = sorted(protein_dic_all_words.items(), key=lambda x: x[1], reverse=True) #To check how many duplicated proteins are, It was sorted with arphabetically order

    all_candidate= ['tail', 'dna', 'baseplate', 'fiber', 'endonuclease', 'polymerase', 'reductase', 'capsid', 'rna', 'wedge', 'terminase', 'assembly', 'head', 'atpase', 'hub', 'prohead', 'recombination', 'aaa', 'conserved', 'helicase', 'transcription', 'core', 'ligase', 'lysis', 'inhibitor', 'regulator', 'hydrolase', 'protease', 'membrane', 'completion', 'kinase', 'neck', 'transcriptional', 'holin', 'sheath', 'portal', 'thioredoxin', 'connector', 'synthase', 'anaerobic', 'binding', 'thymidylate', 'outer', 'factor', 'hnh', 'nuclease', 'internal', 'host', 'ribonucleoside-triphosphate']
    # all_candidate consists of proteins which are over at least 200

    dic_list_all = dict(
        (name, protein_dic_all_words[name]) for name in all_candidate if name in protein_dic_all_words)
    # This code extracts proteins which are over 200

    # This code is used to make wordcloud figure
    #wordcloud = WordCloud(width=1000, height=500).generate_from_frequencies(dic_list_all)

    #plt.figure(figsize=(15, 8))
    #plt.imshow(wordcloud)

    #names = list(dic_list_all.keys())
    #values = list(dic_list_all.values())
    #plt.barh(names, values)
    #plt.ylabel('All', fontweight='bold', color='black', fontsize='15', horizontalalignment='center')
    #plt.yticks(fontsize=10)
    #plt.show()
########################################################################################################################

    removed_protein = [elements for elements in lower_protein_revised if 'hypothetical' not in elements]
    removed_protein_final = [elements for elements in removed_protein if 'putative' not in elements]
    # upper codes remove the useless words 'hypothetical and putative'

    #non_duplicate = list(dict.fromkeys(removed_protein_final))
    removed_protein_final = [elements for elements in removed_protein_final if 'putatove' not in elements]

    list_extend_function=[]
    for words in removed_protein_final:
        words_repeat= words.split(' ')
        words_repeat = [x for x in words_repeat if x != '']
        list_extend_function.extend(words_repeat)

    protein_dic_function_words = {}  # make the dictionary for format like protein_name: sequence IDs (Some sequences can have the same protein)
    for key in list_extend_function:  # Protein name
        if key in protein_dic_function_words:
            protein_dic_function_words[key] += 1
        else:
            protein_dic_function_words[key] = 1

    sorted_proteins_function = sorted(protein_dic_function_words.items(), key=lambda x: x[1], reverse=True)
    #print(sorted_proteins_function)
    function_candidate= ['tail', 'dna', 'baseplate', 'fiber', 'endonuclease', 'polymerase', 'reductase', 'capsid', 'wedge', 'rna', 'terminase', 'assembly', 'atpase', 'head', 'hub', 'prohead', 'aaa', 'recombination', 'lysis', 'core', 'inhibitor', 'transcription', 'ligase', 'hydrolase', 'regulator', 'helicase', 'completion', 'protease',  'kinase', 'sheath', 'thioredoxin', 'connector', 'transcriptional', 'portal', 'synthase', 'anaerobic', 'neck', 'binding']
    print(len(function_candidate))
    dic_list_function = dict(
        (name, protein_dic_function_words[name]) for name in function_candidate if name in protein_dic_function_words)

    #print(dic_list_function)
    #wordcloud = WordCloud(width=1000, height=500).generate_from_frequencies(dic_list_function)

    #plt.figure(figsize=(15, 8))
    #plt.imshow(wordcloud)
    #names = list(dic_list_function.keys())
    #values = list(dic_list_function.values())
    #plt.barh(names, values)
    #plt.ylabel('Function', fontweight='bold', color='black', fontsize='15', horizontalalignment='center')
    #plt.yticks(fontsize=10)
    #plt.show()

########################################################################################################################
# hypothetical classification

    removed_protein_hypothetical = [elements for elements in lower_protein_revised if 'hypothetical' in elements]

    list_extend_hypothetical = []
    for words in removed_protein_hypothetical:
        words_repeat = words.split(' ')
        words_repeat = [x for x in words_repeat if x != '']
        list_extend_hypothetical.extend(words_repeat)

    protein_dic_hypothetical_words = {}  # make the dictionary for format like protein_name: sequence IDs (Some sequences can have the same protein)
    for key in list_extend_hypothetical:  # Protein name
        if key in protein_dic_hypothetical_words:
            protein_dic_hypothetical_words[key] += 1
        else:
            protein_dic_hypothetical_words[key] = 1

    sorted_proteins_hypothetical = sorted(protein_dic_hypothetical_words.items(), key=lambda x: x[1], reverse=True)
    #print(sorted_proteins_hypothetical)
########################################################################################################################
# putative classification

    removed_protein_putative = [elements for elements in lower_protein_revised if 'putative' in elements]
    removed_protein_plus = [elements for elements in lower_protein_revised if 'putatove' in elements]
    removed_protein_putative.extend(removed_protein_plus)


    list_extend_putative = []
    for words in removed_protein_putative:
        words_repeat = words.split(' ')
        words_repeat = [x for x in words_repeat if x != '']
        list_extend_putative.extend(words_repeat)

    protein_dic_putative_words = {}  # make the dictionary for format like protein_name: sequence IDs (Some sequences can have the same protein)
    for key in list_extend_putative:  # Protein name
        if key in protein_dic_putative_words:
            protein_dic_putative_words[key] += 1
        else:
            protein_dic_putative_words[key] = 1

    sorted_proteins_putative = sorted(protein_dic_putative_words.items(), key=lambda x: x[1], reverse=True)
    #print(sorted_proteins_putative)
    putative_candidate=['tail', 'dna', 'endonuclease', 'fiber', 'head', 'helicase', 'recombination', 'holin', 'baseplate', 'membrane', 'neck', 'terminase', 'reductase', 'internal', 'assembly', 'structural', 'outer', 'hnh', 45, 'exonuclease', 'atp-dependent', 'polymerase', 'regulator', 'phosphatase', 'homing', 'anti-sigma', 'protease', 'transcription', 'factor', 'hub', 'transcriptional', 'ligase', 'rna', 'serine/threonine', 'portal',  'prohead', 'capsid', 'methylase', 'tape', 'measure', 'primase', 'ssdna', 'exodeoxyribonuclease', 'srd', 'methyltransferase', 'binding', 'atpase', 'binding,', 'repair,', 'pre-synthesis', 'thymidylate', 'synthase',  'core', 'sir2-like', 'activator', 'middle', 'period']
    print(len(putative_candidate))
    dic_list_putative = dict(
        (name, protein_dic_putative_words[name]) for name in putative_candidate if name in protein_dic_putative_words)


    #wordcloud = WordCloud(width=1000, height=500).generate_from_frequencies(dic_list_putative)

    #plt.figure(figsize=(15, 8))
    #plt.imshow(wordcloud)
    #names = list(dic_list_putative.keys())
    #values = list(dic_list_putative.values())
    #plt.barh(names, values)
    #plt.ylabel('Putative', fontweight='bold', color='black', fontsize='15', horizontalalignment='center')
    #plt.yticks(fontsize=10)
    #plt.show()
########################################################################################################################


