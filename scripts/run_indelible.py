import os
import time
import shutil
import ete3
print ('run it')

print (snakemake.input)

for control_file in snakemake.input:
	taxa_num = control_file.split("indelible/")[1].split("_control.txt")[0]
	taxa_folder = f"./indelible_output/{taxa_num}"
	
	# If taxa folder doesn't exist, create it
	if not os.path.isdir(taxa_folder):
		os.mkdir(taxa_folder)

	print (control_file)

	print (taxa_folder)

	# Copy in the specific INDELible control file and rename it
	shutil.copyfile(control_file, f"{taxa_folder}/control.txt")
	
	# Change into specific INDELible run directory and run INDELible
	os.chdir(taxa_folder)
	os.system('indelible')



while not os.path.exists("trees.txt"):
	time.sleep(1)
for line in open(str("trees.txt")):
	if line.startswith(taxa_num):

		print (line)
		split = line.split("\t")
		rep = split[3]
		tree = ete3.PhyloTree(split[8], format=1)




		if not os.path.exists("trees"):
			os.mkdir(f"./trees/")



		tree.write(outfile=f"./trees/{rep}.nwk", format=3)


		# Add the N0 to the root of the tree
		if not os.path.exists("N0_trees"):
			os.mkdir(f"./N0_trees/")


		os.system(f"sed 's/;/ROOT;/g' ./trees/{rep}.nwk >| ./N0_trees/{rep}.nwk")


