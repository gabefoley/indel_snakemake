print ('run it output')
print (snakemake.input)
print (snakemake.output)
print ('choggers')
print (snakemake.input.tree)
print (snakemake.input.aln)
print (snakemake.wildcards['taxon'])
import ete3
import os


taxon_num = snakemake.wildcards['taxon']
rep = snakemake.wildcards['rep']



print (snakemake.input.aln)
print (os.getcwd())

print (os.path.exists(snakemake.input.aln))


# os.system(f"sed 's/>/>sequence/g; s/*/-/g' {snakemake.input.aln} >| {snakemake.input.aln}")
# os.system (f"trimal -in {snakemake.output} -out {snakemake.output} -noallgaps")


os.system(f"sed 's/\*/-/g;' ./{str(snakemake.input.aln)}")

print ('done')

deva = os.system (f"trimal -in {snakemake.input.aln} -out {snakemake.input.aln} -noallgaps")

print (deva)

print ('cleaning apparently')


# Add seq to the names of the tree (to help GRASP and FastML process the output)
# And change internal numbering to better match GRASP / FastML output
count = 0

tree = ete3.PhyloTree(snakemake.input.tree, snakemake.input.aln, format=1)


for node in tree.traverse("preorder"):

	if not node.is_leaf():
	    node.name = "N" + str(count)
	    count +=1
	else:
		node.name = "sequence" + node.name

tree.get_tree_root().name = "N0"

print ('***')
for g in tree.traverse():
	print (g.name)
	print (g.sequence)

# Get the sequences of the newly renamed nodes

leaves = {x.name : x.sequence for x in tree.traverse() if x.is_leaf()}
renamed = {x.name : x.sequence for x in tree.traverse()}

# Write the cleaned sequence
if not os.path.isdir(f'indelible_output/{taxon_num}/cleaned_aln'):
	os.mkdir(f'indelible_output/{taxon_num}/cleaned_aln')

with open(f'indelible_output/{taxon_num}/cleaned_aln/{rep}.fasta', 'w') as cleaned_aln:
	for name, seq in leaves.items():
		cleaned_aln.write(">" + name + "\n" + seq.replace("*", "-") + "\n")

# Write the cleaned concatenated files
if not os.path.isdir(f'indelible_output/{taxon_num}/concatenated'):
	os.mkdir(f'indelible_output/{taxon_num}/concatenated')

with open(f'indelible_output/{taxon_num}/concatenated/{rep}.fasta', 'w') as concatenated:
	for name, seq in renamed.items():
		concatenated.write(">" + name + "\n" + seq.replace("*", "-") + "\n")


# Write the cleaned tree to path
if not os.path.exists(f'indelible_output/{taxon_num}/cleaned_trees'):
	os.mkdir(f"indelible_output/{taxon_num}/cleaned_trees/")

tree.write(outfile=f"indelible_output/{taxon_num}/cleaned_trees/{rep}.nwk", format=3)


# Write the tree with no internal nodes to path (for FastML)
if not os.path.exists(f'indelible_output/{taxon_num}/cleaned_no_internal_trees'):
	os.mkdir(f"indelible_output/{taxon_num}/cleaned_no_internal_trees/")

tree.write(outfile=f"indelible_output/{taxon_num}/cleaned_no_internal_trees/{rep}.nwk", format=5)


# Write the cleaned tree with added N0 to path
if not os.path.exists(f'indelible_output/{taxon_num}/cleaned_N0_trees'):
	os.mkdir(f"indelible_output/{taxon_num}/cleaned_N0_trees/")


os.system(f"sed 's/;/N0;/g' indelible_output/{taxon_num}/cleaned_trees/{rep}.nwk >| indelible_output/{taxon_num}/cleaned_N0_trees/{rep}.nwk")

