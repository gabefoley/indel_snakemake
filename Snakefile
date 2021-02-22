
rule all:
        input:
            expand("summaries/{taxon}/{method}/{rep}.csv",
                taxon=config['TAXA'], 
                method=config['GRASP_METHODS'], 
                rep = [x for x in range(1, config["REPS"] + 1)])

rule create_control_files:
    output:
        "indelible/{taxon}_control.txt"
    script:
        "scripts/create_control_files.py"

rule run_indelible:
    input:
        "indelible/{taxon}_control.txt"
    output:
        "indelible_output/{taxon}/trees/{rep}.nwk",
        "indelible_output/{taxon}/{taxon}_TRUE_{rep}.fasta",
        "indelible_output/{taxon}/{taxon}_ANCESTRAL_{rep}.fasta",
        "indelible_output/{taxon}/N0_trees/{rep}.nwk",


    script:
        "scripts/run_indelible.py"

rule concat_indelible:
    input:
        extants = "indelible_output/{taxon}/{taxon}_TRUE_{rep}.fasta",
        ancestors = "indelible_output/{taxon}/{taxon}_ANCESTRAL_{rep}.fasta"
    output:
        "indelible_output/{taxon}/concatenated/{taxon}_{rep}.fasta"
    shell:
        "cat {input.extants} {input.ancestors} > {output}"



rule clean_indelible:
    input:
        tree="indelible_output/{taxon}/N0_trees/{rep}.nwk", 
        aln="indelible_output/{taxon}/concatenated/{taxon}_{rep}.fasta"
    output:
        aln="indelible_output/{taxon}/cleaned_aln/{rep}.fasta",
        concat = "indelible_output/{taxon}/concatenated/{rep}.fasta",
        tree="indelible_output/{taxon}/cleaned_trees/{rep}.nwk",
        tree_N0="indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk",
        tree_no_internals="indelible_output/{taxon}/cleaned_no_internal_trees/{rep}.nwk" 

    script:
        "scripts/clean_indelible.py"

rule run_fastml:
    input:
        aln="indelible_output/{taxon}/cleaned_aln/{rep}.fasta",
        tree="indelible_output/{taxon}/cleaned_no_internal_trees/{rep}.nwk"
    output:
        parsimony="fastml_results/{taxon}/{rep}/seq.marginal_Chars_ParsimonyIndels.txt",
        ml="fastml_results/{taxon}/{rep}/seq.marginal_IndelAndChars.txt",
        tree="fastml_results/{taxon}/{rep}/tree.newick.txt",
        dir=directory("fastml_results/{taxon}/{rep}/")
    shell:
        "perl FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File /Users/gabefoley/Dropbox/PhD/20210121_Indel_Evaluation_Project/Indel_Evaluation_Snakemake/{input.aln} --seqType AA --Tree /Users/gabefoley/Dropbox/PhD/20210121_Indel_Evaluation_Project/Indel_Evaluation_Snakemake/{input.tree} --outDir '/Users/gabefoley/Dropbox/PhD/20210121_Indel_Evaluation_Project/Indel_Evaluation_Snakemake/{output.dir}'"

rule clean_fastml_parsimony:
    input:
        aln="fastml_results/{taxon}/{rep}/seq.marginal_Chars_ParsimonyIndels.txt",
        tree="fastml_results/{taxon}/{rep}/tree.newick.txt"
    output:
        aln="fastml_results_parsimony/{taxon}/concatenated/{rep}/FastML_ancestors.fasta",
        cleaned_tree="fastml_results_parsimony/{taxon}/cleaned_trees/{rep}/tree.nwk",
        cleaned_N0_tree="fastml_results_parsimony/{taxon}/cleaned_N0_trees/{rep}/tree.nwk"
    script:
        "scripts/clean_fastml.py"

rule clean_fastml_ml:
    input:
        aln="fastml_results/{taxon}/{rep}/seq.marginal_IndelAndChars.txt",
        tree="fastml_results/{taxon}/{rep}/tree.newick.txt"
    output:
        aln="fastml_results_ml/{taxon}/concatenated/{rep}/FastML_ancestors.fasta",
        cleaned_tree="fastml_results_ml/{taxon}/cleaned_trees/{rep}/tree.nwk",
        cleaned_N0_tree="fastml_results_ml/{taxon}/cleaned_N0_trees/{rep}/tree.nwk"
    script:
        "scripts/clean_fastml.py"


# rule concat_fastml:
#     input:
#         extants="indelible_output/{taxon}/cleaned_aln/{rep}.fasta",
#         ancestors="fastml_results/{taxon}/{rep}/GRASP_ancestors.fasta"
#     output:
#         "grasp_results/concatenated/{taxon}/{method}/{rep}/GRASP_ancestors.fasta"
#     shell:
#         "cat {input.extants} {input.ancestors} > {output}"


rule run_grasp:
    input:
        aln="indelible_output/{taxon}/cleaned_aln/{rep}.fasta",
        tree="indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk",
        # method=expand(config["GRASP_METHODS"])

    output:
        dir=directory("grasp_results/{taxon}/{method}/{rep}"),
        aln="grasp_results/{taxon}/{method}/{rep}/GRASP_ancestors.fasta",
        tree = "grasp_results/{taxon}/{method}/{rep}/GRASP_ancestors.nwk"

    shell:
        "grasp -aln {input.aln} -nwk {input.tree} -model LG -out {output.dir} -savetree  -indel {wildcards.method} -inf joint -gap"

# Add the ancestors and extant sequences to one file
rule concat_fasta:
    input:
        extants="indelible_output/{taxon}/cleaned_aln/{rep}.fasta",
        ancestors="grasp_results/{taxon}/{method}/{rep}/GRASP_ancestors.fasta"
    output:
        "grasp_results/concatenated/{taxon}/{method}/{rep}/GRASP_ancestors.fasta"
    shell:
        "cat {input.extants} {input.ancestors} > {output}"


rule concat_summaries:
    input: 
        expand("summaries/{{taxon}}/{{method}}/{rep}.csv", rep = [x for x in range(1, config['REPS'] + 1)])
    output:
        "concatenated_summaries/{taxon}_{method}.csv"
    shell:
        """awk "FNR==1 && NR!=1{{next;}}{{print}}" {input} > {output}"""

rule get_indelible_dict:
    input:
        indelible_aln_path = "indelible_output/{taxon}/concatenated/{rep}.fasta",
        # indelible_tree_path = "indelible_output/{taxon}/cleaned_trees/{rep}.nwk",
        indelible_tree_path = "indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk"

    output:
        leaves = "indelible_dicts/{taxon}/leaves_{rep}.p",
        gaps = "indelible_dicts/{taxon}/gaps_{rep}.p"
    script:
        "scripts/get_indelible_dict.py"


rule summarise_indels:
    input:
        tree = "grasp_results/{taxon}/{method}/{rep}/GRASP_ancestors.nwk",
        aln = "grasp_results/concatenated/{taxon}/{method}/{rep}/GRASP_ancestors.fasta",
        indelible_tree = "indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk",
        indelible_aln = "indelible_output/{taxon}/concatenated/{rep}.fasta",
        leaves = "indelible_dicts/{taxon}/leaves_{rep}.p",
        gaps = "indelible_dicts/{taxon}/gaps_{rep}.p"

    output:
        "summaries/{taxon}/{method}/{rep}.csv"

    script:
        "scripts/summarise_indels.py"

rule summarise_fastml_parsimony_indels:
    input:
        aln="fastml_results_parsimony/{taxon}/concatenated/{rep}/FastML_ancestors.fasta",
        tree="fastml_results_parsimony/{taxon}/cleaned_N0_trees/{rep}/tree.nwk",
        indelible_tree = "indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk",
        indelible_aln = "indelible_output/{taxon}/concatenated/{rep}.fasta",
        leaves = "indelible_dicts/{taxon}/leaves_{rep}.p",
        gaps = "indelible_dicts/{taxon}/gaps_{rep}.p",
        
    params:
        method='fastMLparsimony'

    output:
        "fastml_parsimony_summaries/{taxon}/fastml/{rep}.csv"

    script:
        "scripts/summarise_indels.py"

rule summarise_fastml_ml_indels:
    input:
        aln="fastml_results_ml/{taxon}/concatenated/{rep}/FastML_ancestors.fasta",
        tree="fastml_results_ml/{taxon}/cleaned_N0_trees/{rep}/tree.nwk",
        indelible_tree = "indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk",
        indelible_aln = "indelible_output/{taxon}/concatenated/{rep}.fasta",
        leaves = "indelible_dicts/{taxon}/leaves_{rep}.p",
        gaps = "indelible_dicts/{taxon}/gaps_{rep}.p",
        
    params:
        method='fastMLml'

    output:
        "fastml_ml_summaries/{taxon}/fastml/{rep}.csv"

    script:
        "scripts/summarise_indels.py"


rule concat_fastml_summaries:
    input: 
        expand("fastml_parsimony_summaries/{{taxon}}/fastml/{rep}.csv", rep = [x for x in range(1, config['REPS'] + 1)])
    output:
        "concatenated_fastml_parsimony_summaries/{taxon}.csv"
    shell:
        """awk "FNR==1 && NR!=1{{next;}}{{print}}" {input} > {output}"""


rule concat_fastml__ml_summaries:
    input: 
        expand("fastml_ml_summaries/{{taxon}}/fastml/{rep}.csv", rep = [x for x in range(1, config['REPS'] + 1)])
    output:
        "concatenated_fastml_ml_summaries/{taxon}.csv"
    shell:
        """awk "FNR==1 && NR!=1{{next;}}{{print}}" {input} > {output}"""




rule summarise_indelible:
    input:
        tree = "indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk",
        aln = "indelible_output/{taxon}/concatenated/{rep}.fasta",
        indelible_tree = "indelible_output/{taxon}/cleaned_N0_trees/{rep}.nwk",
        indelible_aln = "indelible_output/{taxon}/concatenated/{rep}.fasta",
        leaves = "indelible_dicts/{taxon}/leaves_{rep}.p",
        gaps = "indelible_dicts/{taxon}/gaps_{rep}.p"

    output:
        "indelible_summaries/{taxon}/{rep}.csv"

    script:
        "scripts/summarise_indels.py"

rule concat_indelible_summaries:
    input: 
        expand("indelible_summaries/{{taxon}}/{rep}.csv", rep = [x for x in range(1, config['REPS'] + 1)])
    output:
        "concatenated_indelible_summaries/{taxon}.csv"
    shell:
        """awk "FNR==1 && NR!=1{{next;}}{{print}}" {input} > {output}"""



rule generate_output_dict:
    input:
        expand("concatenated_summaries/{taxon}_{method}.csv",
            taxon=config['TAXA'], 
            method=config['GRASP_METHODS']),
        # fastml_parsimony_summary = expand('concatenated_fastml_parsimony_summaries/{taxon}.csv', taxon=config['TAXA']),
        # fastml_ml_summary = expand('concatenated_fastml_ml_summaries/{taxon}.csv', taxon=config['TAXA']),

        indelible_summary = expand('concatenated_indelible_summaries/{taxon}.csv', taxon=config['TAXA'])

    params:
            myparam=lambda wildcards: config["TAXA"]
    output:
        "output_dicts/output.p"
        
    script:
        "scripts/generate_output_dict.py"


rule generate_latex:
    input:
        "output_dicts/output.p"

    output:
        "plots/plot.tex"
    script:
        "scripts/generate_latex.py"

rule compile_latex:
    input:
        "plots/plot.tex"
    output:
        "plots/plot.pdf"

    shell:
        "pdflatex -output-directory plots {input}"


# Constraints for rep and taxon wildcards are that they can only contain digits
wildcard_constraints:
    rep="\d+",
    taxon="\d+"