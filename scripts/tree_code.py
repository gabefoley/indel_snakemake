from ete3 import PhyloTree, TreeStyle, TextFace, add_face_to_node, SeqMotifFace, NodeStyle, faces, ImgFace, CircleFace, \
    AttrFace
import os

def load_tree(tree_path, aln_path=None):
    """
    Load a tree, associate an alignment with it if given
    """

    tree = PhyloTree(tree_path, alignment=aln_path, format=1, alg_format='fasta')
    return tree