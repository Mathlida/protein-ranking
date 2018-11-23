disnerdif = []
with open('PPI_hint.txt') as inf:
    for line in inf:
        item = line.strip().split('\t')
        if item[0] in drigenes and item[1] in difgene:
            disnerdif.append(item[1])
        elif item[0] in difgene and item[1] in drigenes:
            disnerdif.append(item[0])
for gene_hint in genelist_hint:
        dict_all[gene_hint] = 0
with open('PPI_hint.txt') as infile:
    for line in infile:
        first, second, third = line.strip().split('\t')
        if first == second and first in drigenes:
             next
        else:
            if first in disnerdif:
                dict_all[second] += 1
            if second in disnerdif:
                dict_all[first] += 1
with open('rank_LR_eKDPs.txt','w') as of:
    for key, value in dict_all.iteritems():
        of.write("{}\t{}\n".format(key,value))