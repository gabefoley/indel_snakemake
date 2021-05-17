
rule all:
        input:
            expand("output_dicts/{dataset}/gaps_{method}.p",
                dataset=config['DATASETS'], 
                method=config['GRASP_METHODS'])

rule run_grasp:
    input:
        aln="input/{dataset}.aln",
        tree="input/{dataset}.nwk"

    output:
        dir=directory("grasp_results/{dataset}/{method}"),
        aln="grasp_results/{dataset}/{method}/GRASP_ancestors.fasta",
        tree = "grasp_results/{dataset}/{method}/GRASP_ancestors.nwk"

    shell:
        "grasp_indel -aln {input.aln} -nwk {input.tree} -model LG -out {output.dir} -savetree  -indel {wildcards.method} -inf joint -gap -forcelinear -threads 4"

# Add the ancestors and extant sequences to one file
rule concat_fasta:
    input:
        extants="input/{dataset}.aln",
        ancestors="grasp_results/{dataset}/{method}/GRASP_ancestors.fasta",
    output:
        "grasp_results/concatenated/{dataset}/{method}/GRASP_ancestors.fasta",
    shell:
        "cat {input.extants} {input.ancestors} > {output}"


# rule concat_summaries:
#     input: 
#         expand("summaries/{{taxon}}/{{method}}/{rep}.csv", rep = [x for x in range(1, config['REPS'] + 1)])
#     output:
#         "concatenated_summaries/{taxon}_{method}.csv"
#     shell:
#         """awk "FNR==1 && NR!=1{{next;}}{{print}}" {input} > {output}"""


rule summarise_indels:
    input:
        tree = "grasp_results/{dataset}/{method}/GRASP_ancestors.nwk",
        aln="grasp_results/{dataset}/{method}/GRASP_ancestors.fasta",
    output:
        gaps_path = "output_dicts/{dataset}/gaps_{method}.p"
    script:
        "scripts/summarise_indels.py"


# Constraints for rep and taxon wildcards are that they can only contain digits
wildcard_constraints:
    rep="\d+",
    taxon="\d+"