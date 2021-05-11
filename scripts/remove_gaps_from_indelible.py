# import ete3
# import os
# import SeqIO

# from sym import Alphabet

# Protein_wGAP  = Alphabet('ACDEFGHIKLMNPQRSTVWY-acdefghiklmnpqrstvwy')


# print (snakemake.input.extants)
# print (snakemake.input.ancestors)


# extants = sequence.readFastaFile(str(snakemake.input.extants), Protein_wGAP)
# ancestors = sequence.readFastaFile(str(snakemake.input.ancestors), Protein_wGAP)

# for x in extants:
# 	print (x)


from Bio import AlignIO, SeqIO

aa = 'ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy'


extants = AlignIO.read(open(snakemake.input.extants), "fasta")
ancestors = AlignIO.read(open(snakemake.input.ancestors), "fasta")

gap_pos = []
for pos in range(len(extants[0])):


	# print (pos)
	# print (type(extants))
	# print(extants[:, pos])
	# print (type(extants[:, pos]))

	# Check to see if it is all gaps

	all_gaps = True
	for x in extants[:, pos]:
		if x != "-":
			all_gaps = False 
			break
	if all_gaps:
		# print ('YYYYFDFFDFD34343')
		gap_pos.append(pos)


# Remove this column

# print (gap_pos)
offset = 0
for pos in gap_pos:
    pos = pos - offset
    extants = extants[:, :pos ] + extants[:, pos  + 1:]
    ancestors = ancestors[:, :pos ] + ancestors[:, pos  + 1:]
    offset += 1



SeqIO.write(extants, snakemake.output.extants, "fasta")

SeqIO.write(ancestors, snakemake.output.ancestors, "fasta")
