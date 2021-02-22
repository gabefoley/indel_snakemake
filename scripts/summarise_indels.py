import pickle
import gap_checker
import ete3
import csv
from math import log
from scipy.special import rel_entr

# If we include 15 in true and not in recon
# and we include 16 in true and not in recon
# With just a pseudo count 

# Use INDELible max

def get_kl_divergence(recon, true, pseudocount=1):

    # Get the highest length we need to consider. We need to consider cases where one or both of the reconstructions
    # don't actually contain content
    if not true and not recon:
        return 0
    elif not true:
        max_val = max(recon)
    elif not recon:
        max_val = max(true)
    else:
        max_val = max(max(true), max(recon))

    # Get a list of the possible event lengths
    total = [x for x in range(1, max_val + 1)]

    print (f'Total possible insertion lengths {total}\n')

    true_adjust = [true.count(x) + pseudocount for x in total]
    recon_adjust = [recon.count(x) + pseudocount for x in total]

    print(f'True with pseudocount {true_adjust}')
    print(f'Reconstructed with pseudocount {recon_adjust}\n')

    true_probs = [x / sum(true_adjust) for x in true_adjust]
    recon_probs = [x / sum(recon_adjust) for x in recon_adjust]

    print (f'True probs - {true_probs}\n')
    print (f'Reconstructed probs - {recon_probs}\n')

    return sum(recon_probs[i] * log(recon_probs[i]/true_probs[i]) for i in range(len(recon_probs)))



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


outpath = snakemake.output[0]



taxa = snakemake.wildcards['taxon']
print ('chocolate')
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


# indelible_aln_path = f"../Indel_Evaluation_Snakemake/indelible_output/{taxa}/concatenated/{rep}.fasta"
# indelible_tree_path = f"../Indel_Evaluation_Snakemake/indelible_output/{taxa}/cleaned_N0_trees/{rep}.nwk"

# indelible_leaves_path = 'leaves.p'

# indelible_gaps_path = 'gaps.p'

# outpath = "./face"

with open(indelible_leaves_path, 'rb') as handle:
    indelible_leaves_dict = pickle.load(handle)

with open(indelible_gaps_path, 'rb') as handle:
    indelible_gaps_dict = pickle.load(handle)


tree = ete3.PhyloTree(tree_path, aln_path, format=1)

# Get root length
root_len = len(tree.get_tree_root().sequence.replace("-", ""))

tree = ete3.PhyloTree(tree_path, aln_path, format=1)

indelible_tree = ete3.PhyloTree(indelible_tree_path, indelible_aln_path, format=1)

leaves_to_name = {}
gap_dict = {}
ins_accs = []
del_accs = []

for node in tree.traverse():

            if not node.is_leaf():

                leaves = "".join([x.name for x in node.get_descendants() if x.is_leaf()])

                if indelible_leaves_dict[node.name] != leaves:
                    raise RuntimeError(f'{node.name} is different')

                indelible_node = indelible_tree&node.name
                ins_acc, del_acc = gap_checker.get_indel_errors(indelible_node.sequence, node.sequence)
                ins_accs.append(ins_acc)
                del_accs.append(del_acc)

            if not node.is_root():

                parent = node.get_ancestors()[0]
                indels = gap_checker.collect_indels(parent.sequence, node.sequence)
                gap_dict[parent.name + "_" + node.name] = indels




# Get insertion accuracy
avg_ins_acc = sum(ins_accs) / len(ins_accs)
# Get deletion accuracy
avg_del_acc = sum (del_accs) / len(del_accs)

# Get insertion distribution
ins_len = []
del_len = []

for indel_list in gap_dict.values():
    for indel in indel_list:

        if indel.category == 'insertion':
            ins_len.append(int(indel.width))
        elif indel.category == 'deletion':
            del_len.append(int(indel.width))


# Write the details to a summary file


# Calculate KL divergnce
print ('hoggins')
indelible_ins_len = []
indelible_del_len = []
for indel_list in indelible_gaps_dict.values():
    for indel in indel_list:

        if indel.category == 'insertion':
            indelible_ins_len.append(int(indel.width))
        elif indel.category == 'deletion':
            indelible_del_len.append(int(indel.width))

print ('Here are the raw ins and del lengths')
print(ins_len)
print (del_len)
print (indelible_ins_len)
print (indelible_del_len)


ins_kl = get_kl_divergence(ins_len, indelible_ins_len)
del_kl = get_kl_divergence(del_len, indelible_del_len)


with open(outpath, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow(['name', 'taxon', 'method', 'rep', 'root_len', 'ins_acc', 'del_acc', 'ins_len', 'del_len', 'ins_kl', 'del_kl'])
    csvwriter.writerow([outpath, taxa, method, rep, root_len, avg_ins_acc, avg_del_acc, [x for x in ins_len], [x for x in del_len], ins_kl, del_kl])
