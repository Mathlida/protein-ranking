import os
import random


drigenes = [line.strip() for line in open("KDPs.txt").readlines()]
genelist_hint = [line.strip() for line in open("genelist.txt").readlines()]
difgene = [line.strip() for line in open("DEPs.txt").readlines()]

dict_all = {}
for gene_hint in genelist_hint:
    dict_all[gene_hint] = 0
with open('hint_entrz.txt') as infile:
    for line in infile:
        first, second, third = line.strip().split('\t')
        if first == second and first in drigenes:
            next
        else:
            if first in drigenes:
                dict_all[second] += 1
            if second in drigenes:
                dict_all[first] += 1
with open('rank_LR_KDPs.txt','w') as of:
    for key, value in dict_all.iteritems():
        of.write("{}\t{}\n".format(key,value))