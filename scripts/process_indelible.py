import ete3
import os
print ('process')
print (snakemake.wildcards['taxon'])

print (snakemake.input)

taxon = snakemake.wildcards['taxon']

for line in open(str(snakemake.input)):



	if line.startswith(taxon):

			split = line.split("\t")
			rep = split[3]
			tree = ete3.PhyloTree(split[8], format=1)


			if not os.path.exists("indelible_output/{taxon}/trees"):
				os.mkdir(f"indelible_output/{taxon}/trees/")

			print ('written')


			tree.write(outfile=f"indelible_output/{taxon}/trees/{rep}.nwk", format=3)


			# Add the N0 to the root of the tree
			if not os.path.exists("indelible_output/N0_trees"):
				os.mkdir(f"indelible_output/N0_trees/")


			os.system(f"sed 's/;/ROOT;/g' indelible_output/{taxon}/trees/{rep}.nwk >| indelible_output/{taxon}/N0_trees/{rep}.nwk")