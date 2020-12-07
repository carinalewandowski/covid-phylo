covseqs<-read.fasta("final_ny_aligned_name_fixed.txt")

ggphy2 <- read.beast("final_ny_aligned.tree")

p2 = ggtree(ggphy2, mrsd="2020-12-05") + theme_tree2()

plot(ggphy2)
tip <- get.tree(ggphy2)$tip.label
tipcategories = read.csv("sites_28432_28434.csv", 
                         header = TRUE, 
                         stringsAsFactors = FALSE)
newlist <- tipcategories$genotype[order(match(tip, tipcategories$id))]
#tipcategories.rownames<-tipcategories$id

type_list<-c()

for (i in 1:length(tip))
{
 curr_id<-tip[i]
 type_list<-c(type_list, tipcategories$genotype[which(tipcategories$id == curr_id)])
}

beast_tree <- groupOTU(ggphy2, tip[as.factor(type_list)], 
                       group_name = "type")

p2 <- ggtree(beast_tree, aes(color=type), mrsd="2020-12-05") + 
  theme_tree2() + theme(legend.position='none') +
  scale_color_manual(values=c("blue", "red", "green", "yellow"), 
                     labels=c("human", "swine", "a","b"))
p2

new<-match(tipcategories$id, tip)
