import ete3
import os


taxon_num = snakemake.wildcards['taxon']
rep = snakemake.wildcards['rep']


tree = ete3.PhyloTree(snakemake.input.tree, snakemake.input.aln, format=1)


for node in tree.traverse("preorder"):

	if not node.is_leaf():
	    curr = int(node.name.split("N")[1]) 
	    node.name = "N" + str(curr - 1)

renamed = {x.name : x.sequence for x in tree.traverse()}

# # Write the cleaned sequence
# if not os.path.isdir(f'fastml_results_ml/{taxon_num}/concatenated/'):
# 	os.mkdir(f'fastml_results_ml/{taxon_num}/concatenated/')

# if not os.path.isdir(f'fastml_results_ml/{taxon_num}/concatenated/'):
# 	os.mkdir(f'fastml_results_ml/{taxon_num}/concatenated/')

# Write the cleaned sequence
if not os.path.isdir(f'fastml_results_parsimony/{taxon_num}/concatenated/{rep}'):
	os.makedirs(f'fastml_results_parsimony/{taxon_num}/concatenated/{rep}')

if not os.path.isdir(f'fastml_results_ml/{taxon_num}/concatenated/{rep}'):
	os.makedirs(f'fastml_results_ml/{taxon_num}/concatenated/{rep}')


# with open(f'fastml_results/{taxon_num}/concatenated/{rep}/FastML_ancestors.fasta', 'w') as cleaned_aln:
with open(str(snakemake.output.aln[0]), 'w') as cleaned_aln:

	for name, seq in renamed.items():
		cleaned_aln.write(">" + name + "\n" + seq+ "\n")

# if not os.path.isdir(f'fastml_results_parsimony/{taxon_num}/cleaned_trees/'):
# 	print ('lets make it')
# 	os.mkdir(f'fastml_results_parsimony/{taxon_num}/cleaned_trees/')


# if not os.path.isdir(f'fastml_results_ml/{taxon_num}/cleaned_trees/'):
# 	print ('lets make it')
# 	os.mkdir(f'fastml_results_ml/{taxon_num}/cleaned_trees/')



if not os.path.isdir(f'fastml_results_parsimony/{taxon_num}/cleaned_trees/{rep}/'):
	print ('lets make it')
	os.makedirs(f'fastml_results_parsimony/{taxon_num}/cleaned_trees/{rep}/')

if not os.path.isdir(f'fastml_results_ml/{taxon_num}/cleaned_trees/{rep}/'):
	print ('lets make it')
	os.makedirs(f'fastml_results_ml/{taxon_num}/cleaned_trees/{rep}/')


# if not os.path.isdir(f'fastml_results/{taxon_num}/cleaned_N0_tree/{rep}'):
# 	os.mkdir(f'fastml_results/{taxon_num}/cleaned_N0_tree/{rep}')

# tree.write(outfile=f'fastml_results/{taxon_num}/cleaned_trees/{rep}/tree.nwk', format=3)

print ('write tree to ')
print (snakemake.output.cleaned_tree)
tree.write(outfile=snakemake.output.cleaned_tree, format=3)


os.system(f"sed 's/;/N0;/g' {snakemake.output.cleaned_tree} >| {snakemake.output.cleaned_N0_tree}")


# Write the cleaned sequences

# if not os.path.isdir(f'fastml_results_parsimony/{taxon_num}/cleaned_trees/{rep}/'):
# 	print ('lets make it')
# 	os.makedirs(f'fastml_results_parsimony/{taxon_num}/cleaned_trees/{rep}/')

# if not os.path.isdir(f'fastml_results_ml/{taxon_num}/cleaned_trees/{rep}/'):
# 	print ('lets make it')
# 	os.makedirs(f'fastml_results_ml/{taxon_num}/cleaned_trees/{rep}/')





# # Write the cleaned sequence
# if not os.path.isdir(f'indelible_output/{taxon_num}/cleaned_aln'):
# 	os.mkdir(f'indelible_output/{taxon_num}/cleaned_aln')

with open(snakemake.output.aln, 'w') as cleaned_aln:
	for name, seq in renamed.items():
		cleaned_aln.write(">" + name + "\n" + seq.replace("*", "-") + "\n")

