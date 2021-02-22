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

# Write the cleaned sequence
if not os.path.isdir(f'fastml_results/{taxon_num}/concatenated/{rep}'):
	os.mkdir(f'fastml_results/{taxon_num}/concatenated/{rep}')

with open(f'fastml_results/{taxon_num}/concatenated/{rep}/FastML_ancestors.fasta', 'w') as cleaned_aln:
	for name, seq in renamed.items():
		cleaned_aln.write(">" + name + "\n" + seq+ "\n")

if not os.path.isdir(f'fastml_results/{taxon_num}/cleaned_trees/'):
	print ('lets make it')
	os.mkdir(f'fastml_results/{taxon_num}/cleaned_trees/')



if not os.path.isdir(f'fastml_results/{taxon_num}/cleaned_trees/{rep}/'):
	print ('lets make it')
	os.mkdir(f'fastml_results/{taxon_num}/cleaned_trees/{rep}/')

# if not os.path.isdir(f'fastml_results/{taxon_num}/cleaned_N0_tree/{rep}'):
# 	os.mkdir(f'fastml_results/{taxon_num}/cleaned_N0_tree/{rep}')

tree.write(outfile=f'fastml_results/{taxon_num}/cleaned_trees/{rep}/tree.nwk', format=3)

os.system(f"sed 's/;/N0;/g' fastml_results/{taxon_num}/cleaned_trees/{rep}/tree.nwk >| fastml_results/{taxon_num}/cleaned_N0_trees/{rep}/tree.nwk")


# Write the cleaned sequence
if not os.path.isdir(f'indelible_output/{taxon_num}/cleaned_aln'):
	os.mkdir(f'indelible_output/{taxon_num}/cleaned_aln')

# with open(f'indelible_output/{taxon_num}/cleaned_aln/{rep}.fasta', 'w') as cleaned_aln:
# 	for name, seq in leaves.items():
# 		cleaned_aln.write(">" + name + "\n" + seq.replace("*", "-") + "\n")

