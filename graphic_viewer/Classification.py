import random
from collections import namedtuple

import matplotlib
from Bio import SeqIO
import re

Color=list(matplotlib.colors.cnames.values())
input_file="A_baumannii_protein_A_baumannii_KT588073.1.fasta"

seq = SeqIO.parse(input_file,"fasta")
name_separation=[]

name_separation.append(str(input_file.split('_')[3]+'_'+str(input_file.split('_')[4]+'_'+str(input_file.split('_')[5].split('.fasta')[0]))))
#print(name_separation)

ID_information=[]
start_end=[]
description_information=[]
start_point=[]
end_point=[]
for info in seq:
    add_point=info.id.split(':')
    start_end.append(add_point)
    description_information.append(info.description)
    ID_information.append(info.id)
#print(ID_information)

complement_list=[]   # name the strand as '-1' if their names are complement
for i in range(len(ID_information)):
    if 'complement' in ID_information[i]:
        complement_list.append(i)
#print(complement_list)
#print(start_end)
#print(description_information)
for decision in start_end:
    split_point=decision[1].split('..')
    start_point.append(split_point[0].replace("(",""))
    end_point.append(split_point[1].replace(")",""))
#print(start_point)
#print(end_point)

name_only=[]
for name in description_information:
    name_bracket=name.split('|')
    name_add = re.sub("\(.*?\)|\[.*?\]","",name_bracket[1])
    name_only.append(name_add)
#print(name_only)

final_scripts=[]


for i in range(len(start_end)):
    form1='GraphicFeature(start={}, '.format(start_point[i])
    form2='end={}, '.format(end_point[i])
    if i in complement_list:
        form3='strand={}, '.format('-1')
    else:
        form3='strand={}, '.format('+1')
    form4='color="{}", '.format(random.choice(Color))
    if 'ypothetical' in name_only[i]:
        form5='label="{}")'.format("HP")
    else:
        form5 = 'label="{}")'.format(name_only[i])
        if i < len(start_end)-1:
            final_scripts.append(form1+form2+form3+form4+form5+",\n")
        if i==len(start_end)-1:
            final_scripts.append(form1 + form2 + form3 + form4 + form5)


script_file=open('script_of_{}'.format(name_separation[0]),'w')
script_file.writelines(final_scripts)
script_file.close()





