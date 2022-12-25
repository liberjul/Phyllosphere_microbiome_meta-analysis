from Bio import Phylo
import io, copy
with open("../data/phylo/SupFigs1-3-FigTree_abr.txt", "r") as tree_text:
    tree_str = tree_text.readline().strip()

# with open("../data/phylo/SupFigs1-3-FigTree.txt", "r") as tree_text:
#     tree_str = tree_text.readline().strip()
with open("../data/phylo/Order_tip_names.txt", "r") as ifile:
    tip_names = set([x.strip() for x in ifile.readlines()])
tree = Phylo.read(io.StringIO(tree_str), "newick")
leaves = set([x.name for x in tree.get_terminals()])

tips_to_prune = leaves - tip_names

newtree = copy.deepcopy(tree)
for i in tips_to_prune:
    newtree.prune(i)

print(newtree)
Phylo.write(newtree, "../data/phylo/pruned_tree.tre", "newick")
# Phylo.write(newtree, "../data/phylo/pruned_tree_w_annot.tre", "newick")
