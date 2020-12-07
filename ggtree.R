# phy<-read.newick("reduced_tree_NY.nh")

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

# file <- system.file("extdata/BEAST", "final_ny_aligned.tree", package="ggtree")
ggphy <- read.beast("final_ny_aligned.tree") 

tipcategories = read.csv("sites_28432_28434.csv", 
                         header = TRUE, 
                         stringsAsFactors = FALSE)
dd = as.data.frame(tipcategories)

p = ggtree(ggphy, mrsd="2020-12-05") + theme_tree2()

my_colours <- c("black", "orange", "red", "blue")


p %<+% dd + 
  geom_tiplab(color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
              theme(legend.position = c(0.5,0.2),) # no keys

library(ape)
ntips = length(p$tip.label)

#plot(p, tip.color = my_colours[as.factor(tipcategories$genotype)])

# There is an attribute called show.tip.label, so let's set that to false for now
plot(p,show.tip.label=FALSE)
