from ete3 import PhyloTree
import os
import indel_placement
from Bio import AlignIO
import argparse


def load_tree_with_alignment(tree, alignment):
    tree = PhyloTree(tree, alignment=alignment, format=1, alg_format="fasta")
    return tree


def get_positions_with_content(tree, alignment, node_name="#N0"):
    """
    Return the alignment positions that must be present in the ancestor at the given node
    :param tree: The phylogenetic tree
    :param node: The node of interest (defaults to root node)
    :return:
    """

    selected_node = None

    # Get the node of interest
    if node_name == "#N0" or node_name == "N0":
        selected_node = tree.get_tree_root()

    else:
        for node in tree.get_descendants():
            if node.name == node_name:
                selected_node = node
                break

        if selected_node == None:
            raise NameError("The node you selected is not in this tree")

    # Check each position in the alignment to see if any leaf in both children have content there

    # If a position appears in each child tree, no need to check it at higher nodes - it must appear at this node
    found_positions = []

    # If a position appears in neither child tree, no need to check it at higher nodes - it won't appear at this node
    skip_positions = []

    found_positions = sorted(get_found_positions(selected_node, len(alignment[0]), found_positions, skip_positions))

    return found_positions


def get_found_positions(selected_node, aln_len, found_positions, skip_positions):
    # Get the leaf nodes under each child subtree
    child1_leaves = selected_node.children[0].get_leaves()
    child2_leaves = selected_node.children[1].get_leaves()


    for pos in range(aln_len):
        if pos not in found_positions and pos not in skip_positions:
            found = False
            found_multiple = False

            for child in selected_node.children:

                # # If a child is a leaf and has content, automatically elevate the content
                # if child.is_leaf() and child.sequence[pos] != "-":
                #     found_multiple = True
                # else:

                for leaf in child.get_leaves():
                    if leaf.sequence[pos] != "-":
                        if found:
                            found_multiple = True
                            break

                        found = True
                        break


            # Must appear (no need to search in higher nodes)
            if found_multiple:
                found_positions.append(pos)

            # Doesn't appear (no need to search in higher nodes)
            elif not found_multiple:
                skip_positions.append(pos)

    if selected_node.is_root() or len(found_positions) + len(skip_positions) == aln_len:
        return found_positions

    else:
        return get_found_positions(selected_node.up, aln_len, found_positions, skip_positions)


def get_parsimony_table(tree):
    print("Making total parsimony\n")
    parsimony = indel_placement.make_total_parsimony(tree)
    print("Sorting parsimony table\n")
    parsimony.sort_index(inplace=True)

    return parsimony


def get_parsimony_score(parsimony, node, pos, unbalanced=True):
    score = parsimony[pos][node]

    # Only return the score if cost for presence is higher than cost for absence
    if unbalanced:
        if score[0] > score[1]:
            return score

    else:
        return score


def get_parsimony_scores(parsimony, positions, node="#N0"):
    parsimony_dict = {}
    for pos in positions:
        score = get_parsimony_score(parsimony, node, pos)
        if score:
            parsimony_dict[pos] = score
    return parsimony_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--aln", help="Path to alignment", required=True)
    parser.add_argument("-t", "--tree", help="Path to phylogenetic tree", required=True)
    parser.add_argument("-n", "--node", help="Node to return cost for (default is root)")
    parser.add_argument("-p", "--positions", help="Just return the positions required to be there", action="store_true")
    parser.add_argument("-f", "--fasta_output", help="Return all ancestors as a FASTA file")
    parser.add_argument("-to", "--tree_output", help="Write out the ancestor tree")

    args = parser.parse_args()

    # node_of_interest = 0 if not args.node else args.node
    node_of_interest = "#N0" if not args.node else args.node
    print()

    print(u"\U0001F4B0" + " Running Ancestral Cost " + u"\U0001F4B0" + "\n")

    aln = AlignIO.read(args.aln, "fasta")

    tree = load_tree_with_alignment(args.tree, args.aln)

    if args.fasta_output:
        print ("Getting the positions required to be there for all ancestors")
        with open(args.fasta_output, "w+") as fasta_output:
            for node in tree.traverse():
                if not node.is_leaf():
                    print (node.name)
                    positions = get_positions_with_content(tree, aln, node.name)
                    print (positions)



                    ancestor_gaps = ["A" if x in positions else "-" for x in range(len(aln[0]))]

                    fasta_output.write(f'>{node.name}\n{"".join(ancestor_gaps)}\n')



        with open (args.tree_output, "w+") as tree_output:
            print (tree.write(format=7))
            tree_output.write(tree.write(format=1))




    else:

        print("Getting the positions that are required to be there\n")
        positions = get_positions_with_content(tree, aln, node_of_interest)

        print("These positions are required to be there - ", [x + 1 for x in positions])

        # If we just want the positions, print them here and don't continue the rest of the progam
        if not args.positions:

            print("\nLabelling internal nodes\n")
            previous_labels = indel_placement.label_internal_nodes(tree)

            node_of_interest = node_of_interest if not args.node else previous_labels[args.node]

            print("Getting the parsimony scores\n")
            parsimony = get_parsimony_table(tree)

            print("Getting the parsimony scores for required positions\n")
            parsimony_scores = get_parsimony_scores(parsimony, positions, node_of_interest)


            print("Retrieving values for " + node_of_interest + "\n")

            if not parsimony_scores:
                print ("No positions had unbalanced parsimony scores")

            else:


                for k, v in parsimony_scores.items():
                    print(
                        f"Position {k + 1} was required to be in the ancestor, but the parsimony cost for it being there "
                        f"was {v[0]} "
                        f"compared to the parsimony cost for it not being there which was {v[1]}")
                    print()


