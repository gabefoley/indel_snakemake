import pickle
import gap_checker
import ete3

taxa = {snakemake.wildcards['taxon']}
rep = {snakemake.wildcards['rep']}

# indelible_aln_path = list({snakemake.input.indelible_aln_path})[0]
# indelible_tree_path = list({snakemake.input.indelible_tree_path})[0]

# leaves_path = list({snakemake.output.leaves})[0]
# gaps_path = list({snakemake.output.gaps})[0]


indelible_aln_path = snakemake.input.indelible_aln_path
indelible_tree_path = snakemake.input.indelible_tree_path

leaves_path = snakemake.output.leaves
gaps_path = snakemake.output.gaps
print (indelible_aln_path)
print (indelible_tree_path)

indelible_tree = ete3.PhyloTree(indelible_tree_path, indelible_aln_path, format=1)

leaves_to_name = {}
gap_dict = {}

indelible_tree = ete3.PhyloTree(indelible_tree_path, indelible_aln_path, format=1)

for node in indelible_tree.traverse():


            if not node.is_leaf():
                leaves = "".join([x.name for x in node.get_descendants() if x.is_leaf()])
                leaves_to_name[node.name] = leaves

            if not node.is_root():
                parent = node.get_ancestors()[0]
                indels = gap_checker.collect_indels(parent.sequence, node.sequence)
                gap_dict[parent.name + "_" + node.name] = indels


with open(leaves_path, 'wb') as handle:
    pickle.dump(leaves_to_name, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open(gaps_path, 'wb') as handle:
    pickle.dump(gap_dict, handle, protocol=3)
