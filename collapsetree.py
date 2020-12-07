
from ete3 import Tree 
from ete3 import NodeStyle
from ete3 import TreeStyle
from Bio import SeqIO
import time

# Method to collapse trees
def collapse(node):
    # print(node.name)
    
    # Base case
    if (node.is_leaf()):
        return
    
    # Access all children of current node
    childlist = node.children
    firstchild = childlist[0]

    # check that all the children have the same

    # recursively call method on lower level nodes
    for child in childlist:
        collapse(child)

    # Make sure children are leaves, and check that none
    # of the other sisters nodes differ in color
    for child in childlist:
        if not(child.is_leaf()) or child.img_style["fgcolor"] != firstchild.img_style["fgcolor"]:
            return

    # for child in childlist:
    #     child.delete()

    color = firstchild.img_style["fgcolor"]

    # Delete all redundant nodes
    for child in childlist:
        child.delete()
    
    # set color of current leaf to the previous color
    node.img_style["fgcolor"] = color
    node.img_style["shape"] = "sphere"
    node.img_style["size"] = 10