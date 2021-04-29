import pickle
import numpy as np
import os
import time

output_dict_path = snakemake.input[0]

with open(output_dict_path, 'rb') as handle:
    output_dict = pickle.load(handle)

first_entry = next(iter(output_dict))

entry_num = len(output_dict[first_entry]['methods'])

print ('entry num is ' + str(entry_num))

print (output_dict)

with open(snakemake.output[0], "w") as latex_file:


    # Define the ending rules
    midrule = ' \\\\ \midrule\n'
    cmidrule = ' \\\\ \cmidrule(l){2-' + str(entry_num + 2) + '}\n'
    bottomrule = ' \\\\ \bottomrule\n'

    latex_file.write('\documentclass[11pt]{article}'
                     '\n\\usepackage[landscape, paperwidth=80cm, paperheight=50cm, left=0mm, top=0mm, bottom=0mm, right=0mm, margin=5mm]{geometry}'
                     '\n\\usepackage{graphicx}'
                     '\n\\usepackage{booktabs}'
                     '\n\\usepackage{multirow}'
                     '\n\\begin{document}'
                     '\n\\begin{table}[]'
                     '\n\\begin{tabular}{@{}' + 'l' * (entry_num + 2) + '@{}}'
                     '\n\\toprule\n')

    methods_latex = '& '

    for entry in output_dict[first_entry]['methods']:
        methods_latex += " & " + entry

    latex_file.write(methods_latex + midrule)

    for taxa_num in output_dict:

        # Work out which methods have the most similar root lengths to INDELible

        # Get the true root length
        true_root_len = output_dict[taxa_num]['root_lens'][0]

        best_distance = 10000
        best_root_idxs = []
        # Check which others are closest to this
        for idx, entry in enumerate(output_dict[taxa_num]['root_lens'][1:]):
        	print ('entry')
        	print (entry)
        	print ('best distance')
        	print (best_distance)

        	if abs(float(entry) - float(true_root_len)) == best_distance:
        		best_root_idxs.append(idx + 1)


        	elif abs(float(entry) - float(true_root_len)) < best_distance:
        		best_root_idxs = [idx + 1]
        		best_distance = abs(float(entry) - float(true_root_len))


        # Also add the INDELible index to be bolded as well
        best_root_idxs.append(0)

        # Work out which methods have the lowest indel inaccuracies

        print (output_dict[taxa_num]['ins_accs'])
        print (output_dict[taxa_num]['ins_kls'])



        min_ins_acc_val = min(output_dict[taxa_num]['ins_accs'][1:])
        min_ins_acc_idxs = [i for i, x in enumerate(output_dict[taxa_num]['ins_accs']) if x == min_ins_acc_val]

        min_dels_acc_val = min(output_dict[taxa_num]['del_accs'][1:])
        min_dels_acc_idxs = [i for i, x in enumerate(output_dict[taxa_num]['del_accs']) if x == min_dels_acc_val]

        # Work out which methods have the lowest average KL divergence
        min_ins_kls_val = min(output_dict[taxa_num]['ins_kls'][1:])
        min_ins_kls_idxs = [i for i, x in enumerate(output_dict[taxa_num]['ins_kls']) if x == min_ins_kls_val]

        min_dels_kls_val = min(output_dict[taxa_num]['del_kls'][1:])
        min_dels_kls_idxs = [i for i, x in enumerate(output_dict[taxa_num]['del_kls']) if x == min_dels_kls_val]

        # Also add the INDELible index to be bolded as well
        min_ins_acc_idxs.append(0)
        min_dels_acc_idxs.append(0)
        min_ins_kls_idxs.append(0)
        min_dels_kls_idxs.append(0)


        # Write out the LaTeX code
        indel_dist_latex = '\multicolumn{1}{|l|}{\multirow{6}{*}{' + str(taxa_num) + '}} & \multicolumn{1}{l|}{Indel ' \
                                                                                     'distribution}       '

        for entry in output_dict[taxa_num]['indel_dist_paths']:
            indel_dist_latex += '& \multicolumn{1}{l|}{\\includegraphics[width = 80pt, height = 80pt]{' + str(entry) + \
                                '}}'

        # KL Divergence (insertions)

        ins_kls_latex = '\multicolumn{1}{|l|}{} & \multicolumn{1}{l|}{Mean KL-divergence (insertions)}'



        for idx, entry in enumerate(output_dict[taxa_num]['ins_kls']):
            bold = False
            if idx in min_ins_kls_idxs:
                bold = True

            ins_kls_latex += '& \multicolumn{1}{l|}{'
            ins_kls_latex += '\\textbf{' if bold else ''
            ins_kls_latex += str(np.round(entry, decimals=5)) + '}'
            ins_kls_latex += '}' if bold else ''


        # KL Divergence (deletions)

        del_kls_latex = '\multicolumn{1}{|l|}{} & \multicolumn{1}{l|}{Mean KL-divergence (deletions)}'



        for idx, entry in enumerate(output_dict[taxa_num]['del_kls']):
            bold = False
            if idx in min_dels_kls_idxs:
                bold = True

            del_kls_latex += '& \multicolumn{1}{l|}{'
            del_kls_latex += '\\textbf{' if bold else ''
            del_kls_latex += str(np.round(entry, decimals=5)) + '}'
            del_kls_latex += '}' if bold else ''



        root_len_latex = '\multicolumn{1}{|l|}{} & \multicolumn{1}{l|}{Mean root length}'
        
        for idx, entry in enumerate(output_dict[taxa_num]['root_lens']):
            bold = False
            if idx in best_root_idxs:
                bold = True

            root_len_latex += '& \multicolumn{1}{l|}{'
            root_len_latex += '\\textbf{' if bold else ''
            root_len_latex += str(entry) + '}'
            root_len_latex += '}' if bold else ''

        ins_acc_latex = '\multicolumn{1}{|l|}{} & \multicolumn{1}{l|}{Mean insertion inaccuracy}'
        for idx, entry in enumerate(output_dict[taxa_num]['ins_accs']):
            bold = False
            if idx in min_ins_acc_idxs:
                bold = True

            ins_acc_latex += '& \multicolumn{1}{l|}{'
            ins_acc_latex += '\\textbf{' if bold else ''
            ins_acc_latex += str(np.round((entry * 100), decimals=5)) + '}'
            ins_acc_latex += '}' if bold else ''

        print(ins_acc_latex)

        del_acc_latex = '\multicolumn{1}{|l|}{} & \multicolumn{1}{l|}{Mean deletion inaccuracy}'
        for idx, entry in enumerate(output_dict[taxa_num]['del_accs']):
            bold = False
            if idx in min_dels_acc_idxs:
                bold = True

            del_acc_latex += '& \multicolumn{1}{l|}{'
            del_acc_latex += '\\textbf{' if bold else ''
            del_acc_latex += str(np.round((entry * 100), decimals=5)) + '}'
            del_acc_latex += '}' if bold else ''

        # Write out all the parts to the latex file

        latex_file.write(indel_dist_latex + cmidrule)
        latex_file.write(root_len_latex + cmidrule)
        latex_file.write(ins_kls_latex + cmidrule)
        latex_file.write(del_kls_latex + cmidrule)

        latex_file.write(ins_acc_latex + cmidrule)
        latex_file.write(del_acc_latex + midrule)

    # End the table and document
    latex_file.write('\n\\end{tabular}'
                     '\n\\end{table}'
                     '\n\end{document}')


