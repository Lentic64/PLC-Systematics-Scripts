library(ape)

# filenames
target_tree     <- "MyTree1.tre"
bootstrap_trees <- "Combined_Boot_Trees.tre"
output_tree     <- "Mytree2_Bootstraps.tre"
#consensus_tree <- "ConsensusTree.tre"

# read the files
tree     <- read.nexus(target_tree)
bs_trees <- read.nexus(bootstrap_trees)

# compute the bootstraps
bs_scores <- 100 * prop.clades(tree, bs_trees) / length(bs_trees)
round_scores <- round(bs_scores, digits = 0)
tree$node.label <- round_scores

# write the tree
write.tree(tree, file = output_tree)


#create consensus tree
#cons_tree <- consensus(bs_trees, p=0.70, check.labels = TRUE, rooted = FALSE)
#write.tree(cons_tree, file = consensus_tree)


