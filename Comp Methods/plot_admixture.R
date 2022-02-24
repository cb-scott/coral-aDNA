###ANALYZE ADMIXTURE
#Import your NGSadmixture results from the remote pipeline
###what about NGSadmixture

input <- as.matrix(read.table("/Users/carlyscott/Box Sync/Projects/aDNA_local/review_response/amil_admixture.4.qopt"))
#input <- q
rownames(input) <- names
sub_meta <- all_meta[all_meta$Run != c("SRR7235995","SRR7235983"),] #these runs had problems
sub_meta <- sub_meta[sub_meta$Run != "SRR7235983",] #this run had problems - it was mislabeled in the metadata
ord <- sub_meta[order(sub_meta$ANC, sub_meta$Species, sub_meta$Site),1]
tool <- sub_meta[order(sub_meta$ANC, sub_meta$Species, sub_meta$Site),]
tool$num <- 1:nrow(tool)
toplot <- t(input)
toplot <- toplot[,colnames(toplot) != c("SRR7236037","SRR7236038")] #these runs did not map well 
toplot <- toplot[,ord]
ancient <- toplot[,(ncol(toplot)-3):ncol(toplot)]
toplot <- toplot[,1:58] #hard coded 
reorder_anc <- cbind(ancient[,"S17463"], ancient[,"S17464"], ancient[,"S17465"], ancient[,"S17466"])
colnames(reorder_anc) <- c("S17463", "S17464", "S17465", "S17466")
toplot <- cbind(toplot, reorder_anc)

png("~/Box Sync/Projects/aDNA_final_docs/NGSAdmix_k3.png", width = 8, height = 3, units = "in", res = 500)
barplot(toplot, space =0, col = c("royalblue2", "firebrick3", "yellow2", "green", "orange"), xaxt = "n", border = NA)
abline(v = seq(1, 61,  1), col = "grey", lwd = .4 )
abline(v = c(20, 45, 58), lwd =2)
#abline(v = c(6, 10, 15, 27, 34, 39, 45, 51, 52, 58 ), lwd = 2, lty = 'dashed')

dev.off()

#label bars here, if wanted
text(seq(0.5, 61.5,  1),-0.1,ord,xpd=T, srt = 90, cex = .8)

