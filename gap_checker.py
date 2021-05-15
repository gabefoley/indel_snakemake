import tree_code as tc
import functools
from collections import defaultdict
from itertools import chain

@functools.total_ordering
class Indel():
    """
    Class to store information about an indel event
    """
    def __init__(self, category, width, start, end, content, gapped_content, terminal):
        self.category = category
        self.width = width
        self.start = start
        self.end = end
        self.content = content
        self.gapped_content = gapped_content
        self.terminal = terminal

    def __eq__(self, other):
        return ((self.start, self.end) ==
                (other.start, other.end))
    def __lt__(self, other):
        return ((self.start, self.end) <
                (other.start, other.end))

    def __hash__(self):
        return (self.start, self.end).__hash__()

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

def get_gap_dict(tree_path, alignment_paths, names):

    """
    Return a gap dictionary that maps all of the branches in a tree to a dictionary of indels on each branch,
    organised by reconstruction
    :param tree: The tree labelled with ancestral nodes
    :param alignments: The alignments containing ancestors and extant sequences, for as many reconstructions as wanted
    :param names: A name for each reconstruction - must match the number of alignments given
    :return: The gap dictionary
    """

    if len(alignment_paths) != len(names):
        raise AttributeError('Error: alignments should be the same length as names')

    gap_dict = defaultdict(dict)


    for idx, alignment_path in enumerate(alignment_paths):
        tree = tc.load_tree(tree_path, alignment_path)

        # gap_dict[names[idx]] = {}

        for node in tree.traverse("postorder"):
            if not node.is_root():

                parent = node.get_ancestors()[0]
                indels = collect_indels(parent.sequence, node.sequence)
                gap_dict[names[idx]][parent.name + "_" + node.name] = indels

                # gap_dict[parent.name + "_" + node.name] = {}
                gap_dict[parent.name + "_" + node.name][names[idx]] = indels



    return gap_dict

def extract_indel(parent, child, start, end, terminal):
    """
    :param parent: Parent sequence
    :param child: Child sequence
    :param start: Start of indel
    :param end: End of indel
    :return: Indel object representing the extracted indel
    """
    category = 'insertion' if parent[start] == "-" else 'deletion'
    gapped_content = child[start:end] if category == 'insertion' else parent[start:end]
    content = gapped_content.replace("-", "")
    width = len(content)

    indel = Indel(category, width, start, end - 1, content, gapped_content, terminal)

    return indel

def collect_indels(parent_seq, child_seq):
    """
    Collect the set of indels for a given branch, as described by a parent-child relationship
    :param parent_seq: The parent of the branch
    :param child_seq: The child of the branch
    :return: A dictionary mapping the branch name to the set of Indels that are present at that branch
    """

    indels = []
    start = None
    end = None
    category = None
    terminal_pos = True

    for idx, (p, c) in enumerate(zip(parent_seq, child_seq)):


        # If we've found a start and both positions have content then we must have the end
        if start != None and ((p != "-") and (c != "-")):
            end = idx

        # Else, if we've found a start and only one of the positions is a gap we might have the end
        elif start != None and ((p == "-") != (c == "-")):

            # If we're tracking an insertion and the parent has content then we've found an end
            if category == 'insertion' and p != "-":
                end = idx

            # Else, if we're tracking a deletion and the parent has a gap then we've found an end
            elif category == 'deletion' and p == "-":
                end = idx

        # If we don't have a start position and only one of the positions is a gap we've found a start
        if start == None and ((p == "-") != (c == "-")):
            start = idx

            # Work out the category based on parent sequence
            category = 'insertion' if p == "-" else 'deletion'

        # If we have both a start and end we've found an indel
        if start != None and end:

            # If terminal_pos is still true then we must be at the start of one of the sequences
            terminal = 'N' if terminal_pos else 'X'

            # Get the indel
            indel = extract_indel(parent_seq, child_seq, start, end, terminal)
            indels.append(indel)

            #If the position we're currently has only one of the positions as a gap, a new indel starts immediately
            if ((p == "-") != (c == "-")):
                start = idx

                # Work out the category based on the parent sequence
                category = 'insertion' if p == "-" else 'deletion'


            # Otherwise reset the search for the start of an indel
            else:
                start = None

            # Reset the search for the end of an indel
            end = None

        # If both parent and child have content then we can't be terminal anymore
        if ((p != "-") and (c != "-")):
            terminal_pos = False

    # If we didn't find an end position for the indel it must mean we're at the end of one of the sequences
    if start != None and not end:

        # Get the indel
        indel = extract_indel(parent_seq, child_seq, start, idx + 1, 'C')
        indels.append(indel)

    return indels


def get_unique_indels_for_branch(branch, dict):

    unique_indels = defaultdict(list)

    names = [x for x in dict[branch].keys()]


    for indel in dict[branch][names[0]]:
        if indel not in dict[branch][names[1]]:
            unique_indels[names[0]].append(indel)


    for indel in dict[branch][names[1]]:
        if indel not in dict[branch][names[0]]:
            unique_indels[names[1]].append(indel)

    return unique_indels


def get_unique_indels_for_tree(unique_dict):
    unique_set = set()
    for branch in unique_dict.keys():
        print (branch)
        print (unique_dict[branch])
        if unique_dict[branch]:
            for name in unique_dict[branch].keys():
                print (name)
                for indel in unique_dict[branch][name]:
                    unique_set.add(indel)

    return list(unique_set)

def get_unique_indels_for_branch(tree, gap_dict):

    unique_dict = {}

    for node in tree.traverse("postorder"):

        if not node.is_root():


            parent = node.get_ancestors()[0]
            branch = parent.name + "_" + node.name

            unique_dict[branch] = get_unique_indels_for_branch(branch, gap_dict)


    unique_indels_for_tree = get_unique_indels_for_tree(unique_dict)

    for unique_indel in unique_indels_for_tree:
        print (unique_indel.start)
        print (unique_indel.end)


def get_indel_errors(true, reconstructed):
    # Calculate the number of residues present in the reconstructed sequence but absent in the true sequence, divided by the length of the alignment.
    ins_count = 0
    del_count = 0

    for res_pair in zip(true, reconstructed):
        if res_pair[0] != "-" and res_pair[1] == "-":
            ins_count +=1
        elif res_pair[0] == "-" and res_pair[1] != "-":
            del_count +=1

    return (ins_count / len(true), del_count / len(true))




#
#
# gap_dict = get_gap_dict('../tests/files/CYP2U1_just_IDs_10_reconstructed-tree_GRASP.nwk',
#              ['../tests/files/CYP2U1_just_IDs_10_joint-ancestors_GRASP.fasta'], ['grasp_10'])
#
#
#
#
#
#
