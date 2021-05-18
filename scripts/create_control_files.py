print (snakemake.input)
print (snakemake.wildcards)

import os


for outpath in snakemake.output:

	print (outpath)


	taxa_num = snakemake.wildcards['taxon']
	print (taxa_num)
	print (snakemake.config)
	print (snakemake.config['REPS'])
	reps = snakemake.config['REPS']


	with open (f"./{outpath}", "w") as outfile:
		outfile.write("[TYPE] AMINOACID 1\n")
		outfile.write("[SETTINGS] \n [ancestralprint] NEW \n [output] FASTA \n [fastaextension] fasta \n [randomseed] 5622324 \n [markdeletedinsertions] FALSE \n [insertaslowercase] TRUE \n [fileperrep] TRUE \n")
		outfile.write("[MODEL] model1 \n [submodel] WAG \n [indelmodel] POW 1.7 4 \n [indelrate] 0.05 \n")
		# outfile.write("[MODEL] model1 \n [submodel] WAG \n [indelmodel] POW 1.7 20 \n [insertrate] 0.05 [deleterate] 0.01 \n")

		outfile.write(f"[TREE] tree1 \n [rooted] {taxa_num} [treedepth] 1 \n")
		outfile.write("[PARTITIONS] partition1 \n [tree1 model1 200] \n")
		outfile.write(f"[EVOLVE] partition1 {reps} {taxa_num} \n")



	# [TYPE] AMINOACID 1	//  EVERY control file must begin with a [TYPE] command.
	# 			//  The number after "AMINOACID" can be 1 or 2 and chooses the 
	# 			//  algorithm that INDELible uses (see manuscript). Both give 
	# 			//  identical results but in some cases one is quicker.
	# 			//  Other blocks and commands following this statement
	# 			//  can come in any order you like.

	# [SETTINGS]
	#   [ancestralprint] NEW
	#   [output] FASTA
	#   [randomseed] 2232
	#   [markdeletedinsertions] TRUE
	#   [insertaslowercase] TRUE
	#   [fileperrep] TRUE

	# [MODEL]    modelname         //  Evolutionary models are defined in [MODEL] blocks.
	#   [submodel] WAG             //  Here the substitution model is simply set as WAG.
	#   [indelmodel]  POW  1.7 500 //  Power law insertion/deletion length distribution (a=1.7)
	#   [indelrate]   0.1          //  insertion rate = deletion rate = 0.1
	#                              //  relative to average substitution rate of 1.   
	# [TREE] tree1
	#   [unrooted] 2000 

	#   // [unrooted] 10 2.4 1.1 0.2566 0.34  // ntaxa birth death sample mut

	# [PARTITIONS] partitionname             //  [PARTITIONS] blocks say which models go with
	#   [tree1 modelname 500]            //  which trees and define the length of the
	#                                        //  sequence generated at the root (1000 here).

	# [EVOLVE] partitionname 2 outputname  //  This will generate 100 replicate datasets 
	#                                        //  from the [PARTITIONS] block named above.

	# // The true alignment will be output in a file named outputname_TRUE.phy
	# // The unaligned sequences will be output in a file named outputname.fas
	# // To learn how to implement more complicated simulations (or different 
	# // models) please consult the manual or the other example control files.
