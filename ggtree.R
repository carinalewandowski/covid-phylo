phy<-read.newick("reduced_tree_NY.nh")

collapse_identical_tips <- function(phy,tip_label){
  matching_tips <- which(phy$tip.label==tip_label)
  nt <- length(phy$tip.label) # number of tips in tree
  nm <- length(matching_tips) # Number of tips matching the label
  keep <- numeric(nm)
  
  cur_tip <- 1
  while(cur_tip<=nm){
    if(cur_tip == nm){
      keep[cur_tip] <- 1
      break
    }
    next_tip <- cur_tip + 1
    mrca_ <- getMRCA(phy,c(matching_tips[cur_tip],matching_tips[next_tip]))
    descendants <- getDescendants(phy, mrca_)
    descendant_tips <- descendants[descendants<=nt]
    if(all(descendant_tips %in% matching_tips)){
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + length(descendant_tips)
    }else{
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + 1
    }
  }
  to_drop <- matching_tips[!keep]
  new_phy <- drop.tip(phy,to_drop)
  return(new_phy)
}

plot(phy)
phy <- collapse_identical_tips(phy,"GGG")
phy <- collapse_identical_tips(phy,"AAC")

ggphy<-read.tree("final_ny_aligned.nh")

ggtree(ggphy, mrsd="2020-10-18") + 
  theme_tree2() + 
  xlim(2020, 2021)
