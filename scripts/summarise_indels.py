import pickle
import gap_checker
import ete3
import csv
from math import log
from scipy.special import rel_entr


print (snakemake.input)
print (snakemake.output)
print (snakemake.wildcards)

print ('taxon was ')
print (snakemake.wildcards['taxon'])

tree_path = snakemake.input.tree
aln_path = snakemake.input.aln


indelible_tree_path = snakemake.input.indelible_tree
indelible_aln_path = snakemake.input.indelible_aln

indelible_leaves_path = snakemake.input.leaves
indelible_gaps_path = snakemake.input.gaps






taxa = snakemake.wildcards['taxon']
print (snakemake.params)
if 'method' in snakemake.wildcards.keys():
    method = snakemake.wildcards['method']

elif str(snakemake.params):
    print ('yeah it was')
    method = snakemake.params
else:
    method = 'indelible'

print ('method is ')
print (method)

rep = snakemake.wildcards['rep']

tree = ete3.PhyloTree(tree_path, aln_path, format=1)





for node in tree.traverse():

    if not node.is_leaf():  


    if not node.is_root():

        parent = node.get_ancestors()[0]
        indels = gap_checker.collect_indels(parent.sequence, node.sequence)
        gap_dict[parent.name + "_" + node.name] = indels



print ('gap dict is ')

print (gap_dict)


# Write the details to a summary file
with open(snakemake.output.gaps_path, 'wb') as handle:
    pickle.dump(gap_dict, handle, protocol=3)

