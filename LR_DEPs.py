dict_all = {}
for gene_hint in genelist_hint:
    dict_all[gene_hint] = 0
with open('PPI_hint.txt') as infile:
    for line in infile:
        first, second,third = line.strip().split('\t')
        if first == second and first in drigenes:
            next
        else:
            if first in difgene:
                dict_all[second] += 1
            if second in difgene:
                dict_all[first] += 1
with open('rank_LR_DEPs.txt','w') as of:
    for key, value in dict_all.iteritems():
        of.write("{}\t{}\n".format(key,value))