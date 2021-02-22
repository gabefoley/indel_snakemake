from tabulate import tabulate
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from ast import literal_eval
import pickle
import os

# taxa = [6, 8]

# summaries = ['../Indel_Evaluation_Snakemake/concatenated_summaries/6_BEP.csv',
#              '../Indel_Evaluation_Snakemake/concatenated_summaries/6_SICP.csv',
#              '../Indel_Evaluation_Snakemake/concatenated_summaries/6_SICML.csv',
#              '../Indel_Evaluation_Snakemake/concatenated_summaries/8_BEP.csv',
#               '../Indel_Evaluation_Snakemake/concatenated_summaries/8_SICP.csv',
#               '../Indel_Evaluation_Snakemake/concatenated_summaries/8_SICML.csv'
#              ]

print (snakemake.params)
print (snakemake.input)
print (snakemake.input.indelible_summary)
taxa = [int(x) for x in snakemake.params[0]]
summaries = snakemake.input 
print (summaries)
print (list(summaries))

print ('FIRE###############')
print (snakemake.input.indelible_summary)

output_path = snakemake.output[0]
output_dict = defaultdict(dict)

def generate_indel_distribution_plot(ins_lens, del_lens, method_name, taxa_num, outpath):
    # ins_len = [x.width for x in ins_lens]
    # del_len = [x.width for x in del_lens]

    print (ins_lens)

    ins_count = len(ins_lens)
    del_count = len(del_lens)

    fig = plt.figure(figsize=(20, 12))

    sns.distplot(ins_lens, color='purple');
    sns.distplot(del_lens, color='blue');

    # fig.legend(labels=[f'{method_name} insertions ({ins_count} total)', f'{method_name} deletions ({del_count} total)'])
    # fig.suptitle(f'{method_name} with {taxa_num} taxa')

    plt.savefig(outpath)





for taxa_num in taxa:
    methods = []
    indel_dist_paths = []
    root_lens = []
    ins_accs = []
    del_accs = []
    ins_kls = []
    del_kls = []


    for summary_path in reversed(summaries):

        summary = pd.read_csv(summary_path, converters={'ins_len' : eval,
                                                        'del_len' : eval})

        print (summary)

        if summary['taxon'][0] == taxa_num:

            method = summary['method'][0]

            print ('method is ')
            print (method)



            methods.append(method)



            root_lens.append(summary['root_len'].mean())
            ins_accs.append(summary['ins_acc'].mean())
            del_accs.append(summary['del_acc'].mean())

            # The insertion and deletion distributions
            ins_lens = summary['ins_len'].sum()
            del_lens = summary['del_len'].sum()

            # The KL divergences
            ins_kls.append(summary['ins_kl'].mean())
            del_kls.append(summary['del_kl'].mean())

            if not (os.path.exists('images')):
                os.mkdir('images')

            # The path to the distribution plot that we will store
            indel_dist_path = f'images/{method}_{taxa_num}.png'


            generate_indel_distribution_plot(ins_lens,
                                             del_lens,
                                             method, taxa_num, indel_dist_path)

            indel_dist_paths.append(indel_dist_path)


    print ('and now methods is ')
    print (methods)
    output_dict[taxa_num]['methods'] = methods
    output_dict[taxa_num]['root_lens'] = root_lens
    output_dict[taxa_num]['ins_accs'] = ins_accs
    output_dict[taxa_num]['del_accs'] = del_accs
    output_dict[taxa_num]['ins_kls'] = ins_kls
    output_dict[taxa_num]['del_kls'] = del_kls
    output_dict[taxa_num]['indel_dist_paths'] = indel_dist_paths

print ('output dict')

print (output_dict)

with open(output_path, 'wb') as handle:
    pickle.dump(output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


