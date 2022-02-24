#Projection-Based-PCA
#Input: Your merged evec file from the remote pipeline
#This file contains the PC coordinates for each sample (modern and ancient)
#As well as each sample with chromosome-level jackknifing. The jackknifed samples
#will form the error bars in this plot
#Create 

nchroms = 16 #how many chromosomes did your reference have?

#get your data input
a<-read.table("~/Box Sync/Projects/aDNA_local/review_response/merged.missing_mill.evec.txt",sep="\t",header=T, stringsAsFactors = F)
#get sample IDs
id_splitr <- function(s){
  out <- ifelse(grepl('\\.', s) == T, unlist(strsplit(s, "\\."))[2], "full")
  return(out)
}

#Complete the jackknife computation
jackknife <- function(x,print=T)
{
  c <- mean(x)
  d <- 0
  for (j in 1:length(x))
  {
    d<-d + (x[j]-c)*(x[j]-c)
  }
  sqrt(d*(nchroms-1)/nchroms)
}

#give each sample a replicate number
a$rep <- sapply(a$Sample_ID, id_splitr)

#which samples have all chromosomes?
all_chrs <- a[a$rep == "full",]

#which samples were jackknifed?
jk_chrs <- a[a$rep != "full",]
jk_chrs$ID <- sapply(jk_chrs$Sample_ID, function(x) unlist(strsplit(x, "\\."))[1])


library(tidyverse)
library(ggrepel)
#this is the error about their true locations.... which we get from the orig 
adjusted1 <- jk_chrs %>% group_by(ID) %>% summarise(jkpc1 = jackknife(pc1), jkpc2 = jackknife(pc2), jkpc3 = jackknife(pc3)) %>% 
  merge(., all_chrs[,1:4], by.x = "ID", by.y="Sample_ID") %>% mutate(h_adj_pc1 = pc1 + jkpc1, l_adj_pc1 = pc1 - jkpc1, 
                                                                     h_adj_pc2 = pc2 + jkpc2, l_adj_pc2 = pc2 - jkpc2, 
                                                                     h_adj_pc3 = pc3 + jkpc3, l_adj_pc3 = pc3 - jkpc3)


#grab our meta data, merge it with our PCA locations
meta <- read.csv("~/Box Sync/Projects/aDNA_final_docs/Kitchen_meta_edit.csv", stringsAsFactors = F, header = T) 
meta <- meta %>% separate(LibraryName, into=c("Species", "Site", "ID"), sep = "_")%>% select(Run, Species, Site)
meta$ANC <- "N"
anc<- data.frame(Run=c("S17463", "S17464", "S17465", "S17466"), Species=c("Apalm", "Aspp", "Apalm", "Apalm"), Site=rep("FL", 4), ANC=rep("Y", 4))
all_meta <- rbind(meta, anc)

new_meta <- merge(all_meta, adjusted1, by.x="Run", by.y = "ID")

new_meta$labels <- new_meta$Run
new_meta$labels[new_meta$ANC == "N"] <- ''

#plot these values, with error bars.
#Here I have them colored by site and species. You can change this to whatever you would like

site_pca <- new_meta %>% ggplot(aes(x = pc1, y = pc2, xmin = l_adj_pc1, xmax = h_adj_pc1, 
                        ymin = l_adj_pc2, ymax = h_adj_pc2, col = Site)) + geom_point(size=6, alpha = .4) +
  geom_errorbarh()+ geom_errorbar() + theme_minimal()  + xlab("PC1") + ylab("PC2") +
  geom_label_repel(aes(label=labels), min.segment.length = .25, nudge_x = -.01) + coord_equal() + scale_color_manual(values=c("#99EDCC", "#E2C044", "#034C3C", "#4E5340", "#76B041"))
ggsave("~/Box Sync/Projects/aDNA_final_docs/Reich-based_pca-site.png", site_pca, width = 4, height =  4.5)


species_pca <- new_meta %>% ggplot(aes(x = pc1, y = pc2, xmin = l_adj_pc1, xmax = h_adj_pc1, 
                                    ymin = l_adj_pc2, ymax = h_adj_pc2, col = Species)) + geom_point(size=6, alpha = .4) +
  geom_errorbarh()+ geom_errorbar() + theme_minimal()  + xlab("PC1") + ylab("PC2") +
  geom_label_repel(aes(label=labels), min.segment.length = .25, nudge_x = -.01) + coord_equal() + scale_color_manual(values=c("#99EDCC", "#E2C044", "#034C3C", "#4E5340", "#76B041"))
