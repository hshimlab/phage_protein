from Bio import SeqIO
import re
import itertools
from operator import itemgetter
import matplotlib.pyplot as plt
from wordcloud import WordCloud
# file name = "A_baumannii_phage_protein.fasta"
# Whole proteins number : 22297
# number of hypothetical :
# number of all proteins : 8571 ->1189
# number of function proteins : 793
# number of putative proteins :396
def fasta_seq(file_name):
    all_species=[]    # remove the '|' in the description part
    protein_name=[]   # make the list for classifying protein Name
    Number_of_ID=[]   # make the list for classifying sequence ID

    for seq_record in SeqIO.parse(file_name, "fasta"):   # open the fasta file
        all_species.append(seq_record.description.split('|')) # remove '|' in the description part
        Number_of_ID.append(seq_record.id) # add the sequence IDs to the list

    for protein in all_species:
        protein_name.append(protein[1])     # extract the protein name from the description part

    protein_revised = []
    for candidate in protein_name:  # This code is used to remove symbols which disturb interpreting
        re_candidate = re.sub("[\(\[].*?[\)\]]", "",
                              candidate)  # This code remove bracket [] which contains the kinds of bacteria
        protein_revised.append(re_candidate)

    lower_protein_revised = []  # Python distinguishes alphabet with lower case and upper case
    for lower_alphabet in protein_revised:  # This code is used to change the alphabet form with lower cases
        lower_alphabet = lower_alphabet.lower()
        lower_protein_revised.append(lower_alphabet)

    removed_protein = [elements for elements in lower_protein_revised if 'hypothetical' not in elements]

########################################################################################################################
# All proteins example

    protein_dic_all_words = {}  # To count the duplicate words, make the dictionary
    for key in removed_protein:
        if key in protein_dic_all_words:
            protein_dic_all_words[key]+=1  # If duplicated proteins are, add +1 for counting
        else:
            protein_dic_all_words[key]=1   # If there are not duplicated proteins, remain as 1

    sorted_proteins_all = sorted(protein_dic_all_words.items(), key=lambda x: x[1],
                                 reverse=True)  # To check how many duplicated proteins are, It was sorted with arphabetically order

    For_bar_chart_all = dict(sorted(protein_dic_all_words.items(), key=itemgetter(1), reverse=True)[:30]) # print 30 elements as the form of bar_chart in ascending order

    #names = list(For_bar_chart_all.keys())
    #values = list(For_bar_chart_all.values())
    #plt.barh(names, values)
    #plt.ylabel('All', fontweight='bold', color='black', fontsize='15', horizontalalignment='center')
    #plt.yticks(fontsize=6)
    #plt.show()


    list_split_all = []  # make the list for all kind of proteins
    for words in removed_protein:
        words_repeat = words.split(' ')
        words_repeat = [x for x in words_repeat if x != '']  # many proteins contains ''(not meaningful), and remove them
        list_split_all.extend(words_repeat)  # merge the multiple lists as one list

    protein_dic_all_words_split = {}  # To count the duplicate words, make the dictionary
    for key in list_split_all:
        if key in protein_dic_all_words_split:
            protein_dic_all_words_split[key] += 1  # If duplicated proteins are, add +1 for counting
        else:
            protein_dic_all_words_split[key] = 1  # If there are not duplicated proteins, remain as 1

    sorted_proteins_all_split = dict(sorted(protein_dic_all_words_split.items(), key=itemgetter(1), reverse=True)[:30])  # To check how many duplicated proteins are, It was sorted with arphabetically order

    # This code is used to make wordcloud figure
    #wordcloud = WordCloud(width=1000, height=500).generate_from_frequencies(sorted_proteins_all_split)

    #plt.figure(figsize=(15, 8))
    #plt.imshow(wordcloud)


#########################################################################################################################
# Function proteins example

    removed_protein_function = [elements for elements in removed_protein if 'putative' not in elements]

    protein_dic_function_words = {}  # To count the duplicate words, make the dictionary
    for key in removed_protein_function:
        if key in protein_dic_function_words:
            protein_dic_function_words[key] += 1  # If duplicated proteins are, add +1 for counting
        else:
            protein_dic_function_words[key] = 1  # If there are not duplicated proteins, remain as 1

    sorted_proteins_function = sorted(protein_dic_function_words.items(), key=lambda x: x[1],
                                 reverse=True)  # To check how many duplicated proteins are, It was sorted with alphabetically order

    For_bar_chart_function = dict(sorted(protein_dic_function_words.items(), key=itemgetter(1), reverse=True)[:30]) # print 30 elements as the form of bar_chart in ascending order

    # print(For_bar_chart_All)
    #names = list(For_bar_chart_function.keys())
    #values = list(For_bar_chart_function.values())
    #plt.barh(names, values)
    #plt.ylabel('Function', fontweight='bold', color='black', fontsize='15', horizontalalignment='center')
    #plt.yticks(fontsize=6)
    #plt.show()


    list_split_function = []  # make the list for all kind of proteins
    for words in removed_protein_function:
        words_repeat = words.split(' ')
        words_repeat = [x for x in words_repeat if
                        x != '']  # many proteins contains ''(not meaningful), and remove them
        list_split_function.extend(words_repeat)  # merge the multiple lists as one list

    protein_dic_function_words_split = {}  # To count the duplicate words, make the dictionary
    for key in list_split_function:
        if key in protein_dic_function_words_split:
            protein_dic_function_words_split[key] += 1  # If duplicated proteins are, add +1 for counting
        else:
            protein_dic_function_words_split[key] = 1  # If there are not duplicated proteins, remain as 1

    sorted_proteins_function_split = dict(sorted(protein_dic_function_words_split.items(), key=itemgetter(1), reverse=True)[:30])  # To check how many duplicated proteins are, It was sorted with alphabetically order

    # This code is used to make wordcloud figure
    #wordcloud = WordCloud(width=1000, height=500).generate_from_frequencies(sorted_proteins_function_split)

    #plt.figure(figsize=(15, 8))
    #plt.imshow(wordcloud)


########################################################################################################################
# Putative proteins example

    removed_protein_putative = [elements for elements in removed_protein if 'putative' in elements]

    protein_dic_putative_words = {}  # To count the duplicate words, make the dictionary
    for key in removed_protein_putative:
        if key in protein_dic_putative_words:
            protein_dic_putative_words[key] += 1  # If duplicated proteins are, add +1 for counting
        else:
            protein_dic_putative_words[key] = 1  # If there are not duplicated proteins, remain as 1

    sorted_proteins_putative = sorted(protein_dic_putative_words.items(), key=lambda x: x[1],
                                      reverse=True)  # To check how many duplicated proteins are, It was sorted with arphabetically order

    For_bar_chart_putative = dict(sorted(protein_dic_putative_words.items(), key=itemgetter(1), reverse=True)[:30])

    #names = list(For_bar_chart_putative.keys())
    #values = list(For_bar_chart_putative.values())
    #plt.barh(names, values)
    #plt.ylabel('Putative', fontweight='bold', color='black', fontsize='15', horizontalalignment='center')
    #plt.yticks(fontsize=4)
    #plt.show()

    list_split_putative = []  # make the list for all kind of proteins
    for words in removed_protein_putative:
        words_repeat = words.split(' ')
        words_repeat = [x for x in words_repeat if
                        x != '']  # many proteins contains ''(not meaningful), and remove them
        list_split_putative.extend(words_repeat)  # merge the multiple lists as one list

    protein_dic_putative_words_split = {}  # To count the duplicate words, make the dictionary
    for key in list_split_putative:
        if key in protein_dic_putative_words_split:
            protein_dic_putative_words_split[key] += 1  # If duplicated proteins are, add +1 for counting
        else:
            protein_dic_putative_words_split[key] = 1  # If there are not duplicated proteins, remain as 1

    sorted_proteins_putative_split = dict(
        sorted(protein_dic_putative_words_split.items(), key=itemgetter(1), reverse=True)[
        :30])  # To check how many duplicated proteins are, It was sorted with alphabetically order

    # This code is used to make wordcloud figure
    #wordcloud = WordCloud(width=1000, height=500).generate_from_frequencies(sorted_proteins_putative_split)

    #plt.figure(figsize=(15, 8))
    #plt.imshow(wordcloud)


########################################################################################################################

















