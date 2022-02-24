###PCANGSD based pipeline for constructing PCAs and dealing with ADMIXTURE results 
library(tidyverse)
library(ggrepel)
#get covariance matrix output by pcangsd in remote pipeline
C <- as.matrix(read.table("~/Box Sync/Projects/aDNA_local/review_response/pcangsd.cov"))
e <- eigen(C)


#let's get my bamlist- these are your sample names ordered as input into pcangsd
sample_list <- read.table("~/Box Sync/Projects/aDNA_local/review_response/bamlist", stringsAsFactors = F)

names <- sample_list$V1 #names are the names of each sample in the bam list. I did this through bash remotely. 

e_vectors <- data.frame(e$vectors)
e_vectors$names <- names
meta <- read.csv("~/Box Sync/Projects/aDNA_final_docs/Kitchen_meta_edit.csv", stringsAsFactors = F, header = T) 

meta <- meta %>% separate(LibraryName, into=c("Species", "Site", "ID"), sep = "_")%>% select(Run, Species, Site)
meta$ANC <- "N"

#ancient data 
anc<- data.frame(Run=c("S17463", "S17464", "S17465", "S17466"), Species=c("Apalm", "Aspp", "Apalm", "Apalm"), Site=rep("FL", 4), ANC=rep("Y", 4))
all_meta <- rbind(meta, anc)

e_meta <- merge(e_vectors, all_meta, by.x = "names", by.y = "Run")

modnames <- meta$Run
ancnames <- anc$Run
e_meta$X1 <- -e_meta$X1 #flip PC for visual 
sub_anc <- e_meta[e_meta$names %in% ancnames,]
apca <- ggplot() + geom_point(data = e_meta[e_meta$names %in% modnames, ], aes(x=X1, y=X2, col = Species), size = 6, alpha = .4) + theme_bw() +
  scale_color_manual(values=c( "royalblue2","firebrick3", "purple4")) + 
  xlab("PC1") + ylab("PC2") +
  #theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18)) + geom_point(data = sub_anc, aes(x=X1, y =X2), size = 6.5, col = "black", alpha = .8, pch = 24, fill = "green4", stroke = 1.1) + 
  geom_label_repel(data = sub_anc, aes(x=X1, y =X2, label=names), min.segment.length = .25, nudge_x = -.01) + coord_equal()
ggsave("~/Box Sync/Projects/aDNA_final_docs/ANGSD_PCA_amil_bysite.png", apca, width = 4, height =  4.5)

###IF YOU ARE INTERESTED IN PCs BEYOND PC1 & PC2:
apca3 <- ggplot() + geom_point(data = e_meta[e_meta$names %in% modnames, ], aes(x=X1, y=X3, col = Species), size = 6, alpha = .4) + theme_bw() +
  scale_color_manual(values=c( "royalblue2","firebrick3", "purple4")) + 
  xlab("PC1") + ylab("PC3") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18)) + geom_point(data = sub_anc, aes(x=X1, y =X3), size = 6.5, col = "black", alpha = .8, pch = 24, fill = "green4", stroke = 1.1) + 
  geom_label_repel(data = sub_anc, aes(x=X1, y =X3, label=names), min.segment.length = .25, nudge_x = -.01) + coord_equal()

apca23 <- ggplot() + geom_point(data = e_meta[e_meta$names %in% modnames, ], aes(x=X2, y=X3, col = interaction(Species,Site)), size = 6, alpha = .4) + theme_bw() +
  #scale_color_manual(values=c( "royalblue2","firebrick3", "purple4")) + 
  xlab("PC2") + ylab("PC3") +
  #theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18)) + geom_point(data = sub_anc, aes(x=X2, y =X3), size = 6.5, col = "black", alpha = .8, pch = 24, fill = "green4", stroke = 1.1) + 
  geom_label_repel(data = sub_anc, aes(x=X2, y =X3, label=names), min.segment.length = .25, nudge_x = -.01) + coord_equal()


###IF YOU DID HAPLOTYPE CALLS IN ANGSD
#This pipeline is not given in the remote runs - it yielded the same result as 
#the Projection-Based-Pipeline. If you are interested in how to set this up
#please reach out to me: cbscott@utexas.edu. 
#I recommend that you use the given Projection-Based-Pipeline, as setting this up 
#through ANGSD has several dependencies which are a bit of a pain to get working
#together

C <- as.matrix(read.table("haplosub.cov"))
e <- eigen(C)
anc <- c(rep("red", 4), rep("blue", 60))
plot(e$vectors[,c(1,3)],col=anc,xlab="PC1",ylab="PC2", main="individual allele frequency") #check PCA







       