import random
from collections import namedtuple
import glob
import os
import matplotlib
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
os.chdir(r'/Users/file location') # this code is used to run all files in the specific directory
input_list= glob.glob('*.fasta') # select all fasta files to run graphic viewer

from dna_features_viewer import (GraphicFeature, GraphicRecord, CircularGraphicRecord) # graphic viewer package


for input_file in input_list: # run multiple files by using for loop
    Color=list(matplotlib.colors.cnames.values()) # to give the graphic color, it generates all random color codes

    seq = SeqIO.parse(input_file,"fasta") # start to the read fasta files

    name_separation=[] # to save the figure with the original fasta file name
    ID_information=[] # this list is used to check if the start and end points were well separated
    start_end=[] # this part is needed to separate the format of '>KY268295.1:941..1099'
    description_information=[] # this list is used to check description part of the specific sequence and will give the name to further process
    start_point=[] # get the start point which starts translation
    end_point=[] # get the end point which finishes translation

    name_separation.append(str(input_file.split('.fasta')[0])) # for example, the file name is 'A_baumannii_protein_A_baumannii_KY082667.1.fasta'.
                                                               # So, the code splits name with '.fasta' and picks first element in the list


    for info in seq: # In this loop, sequences which contain 'hypothetical' will be removed
        if 'ypothetical' not in info.description:  # the word 'ypothetical' was used because python can not distinguish Hypothetical and hypothetical
            add_point = info.id.split(':')  # the format of '>KY268295.1:941..1099' is separated by symbol ':'
            start_end.append(add_point)  # save the content which form is '941..1099'
            description_information.append(info.description)  # save the data
            ID_information.append(info.id)  # save the ID of the specific sequence


    complement_list = []  # name the strand as '-1' if their names contain complement
    for i in range(len(ID_information)):
        if 'complement' in ID_information[i]:
            complement_list.append(-1)  # if sequence ID contain 'complement' add '-1' for further process
        else:
            complement_list.append(1)  # else, list will add '+1' for further process

    for decision in start_end:  # this loop is used to separate start_end list's elements
        split_point = decision[1].split('..')  # '941..1099' This form will be separated into start and end point
        start_point.append(split_point[0].replace("(", "").replace("<", ""))
        end_point_remove_overlapped = split_point[1].replace(")", "").replace(">", "").split(',')[0]  # remove overlapped proteins -> some fasta files contain overlapped proteins, these proteins will be treated
        end_point.append(end_point_remove_overlapped)  # -> some fasta files contain overlapped proteins, these proteins will be treated at the further study


    name_only = []  # This list is used to give the input data, when graphic viewer gives the name on the box
    for name in description_information:  # from the description_information list, it can extract name for the graphic viewer box
        name_bracket = name.split('|')  # list has the form of the '>KY268295.1:2411..2524 |putative membrane protein [Acinetobacter phage vB_AbaP_AS12]'
        name_add = re.sub("\(.*?\)|\[.*?\]", "", name_bracket[1])  # get the form of 'putative membrane protein [Acinetobacter phage vB_AbaP_AS12]' and remove the bracket []
        name_add_extra = re.sub(", partial", "", name_add)  # some fasta files contain 'partial', so remove the word 'partial'
        name_only.append(name_add_extra)

    if end_point == [] or () or None:  # This code was used to ignore zero-byte file, if zero-byte files were removed, you don't need to used this code
        continue

    features_list = []
    # for example, below code is the original example of graphic viewer
    # [GraphicFeature(start=0, end=20, strand=+1, color="#ffd700",
    #                    label="Small feature"),
    #     GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
    #                    label="Gene 1 with a very long name"),
    #     GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
    #                    label="Gene 2"),
    #     GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
    #                    label="Gene 3")]
    # If above code is directly input, it will be very hard, So, features_list was used to make format of above code for every fasta files
    for i in range(len(start_point)):
        features_list.append(GraphicFeature(start=int(start_point[i]), end=int(end_point[i]), strand=int(complement_list[i]),color=random.choice(Color), label=str(name_only[i])))
        # By using for loop, the format of function (features) will be made
        # ->GraphicFeature(start=0, end=20, strand=+1, color="#ffd700", label="Small feature"),

    features = features_list  # this code can be recognized as original form, in addition, it can recognize more range than example

    name = name_separation[0]  # This list only has one element (fasta file name) -> [0] is used to extract name as the type 'string'


    record = GraphicRecord(sequence_length=max([int(i) for i in end_point]), features=features)  # this code is used to make linear graph
    # the code 'sequence_length=max([int(i) for i in end_point]' was used, because there were some endpoints which were not ordinary sorted
    # To get the maximum value of the end points, the process of sorting was used.
    ax, _ = record.plot(figure_width=150)  # this size is the maximum value for making the graph without assertion error

    ax.figure.savefig(f'{name}.png')  # this code save the file as the name of fasta file
    record.plot(figure_width=150)

