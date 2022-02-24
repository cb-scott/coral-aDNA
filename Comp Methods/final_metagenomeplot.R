#final metagenome plot 
library(stringr)
library(tidyverse)
library(tidyr)

setwd("~/Box Sync/Projects/aDNA_local/review_response/")

fl_meta <- read.csv("PRJNA546259.fullMeta.csv", header = T, stringsAsFactors = F)
westrich_meta <- read.csv("PRJNA299413.fullMeta.csv", header = T, stringsAsFactors = F)
westrich_coral <- westrich_meta[grepl("pox", westrich_meta$SampleName) == F & grepl("H2O", westrich_meta$SampleName) == F & westrich_meta$SampleName != "SampleName",]
westrich_coral <- westrich_coral[grepl("MR", westrich_coral$SampleName) == T,]
#want to divide by reef again and check overlap
create_table <- function(filepath){
all_fl <- list.files(filepath)
read_list <- list()
for(file in 1:length(all_fl)){
  read_list[[file]] <- read.table(paste(filepath, all_fl[[file]], sep = ''))
}

nam <- vector()
for(i in 1:length(all_fl)){
nam <- c(nam, strsplit(all_fl[i], "[.]")[[1]][1])
}

for(i in 1:length(read_list)){
  read_list[[i]]$name <- nam[i]
}
long_df <- bind_rows(read_list, .id = "column_label")
return(long_df)
}

fl_list <- create_table("scaffold_counts_microbe/")
west_list <- create_table("scaffold_stats_westrich/")


#all acroporid metadata 
fl_count_by_scaff <- fl_list %>% filter(name %in% fl_meta$Run) %>% group_by(V1) %>% summarise(count = sum(V3)) %>% filter(count != 0)
west_cbs <- west_list %>% filter(name %in% westrich_coral$Run) %>% group_by(V1) %>% summarise(count = sum(V3)) %>% filter(count != 0)

#ancient metagenomics 
anc_MAGs <- read.csv("~/Box Sync/Projects/aDNA_final_docs/Microbes_by_ancient.csv", header = F, stringsAsFactors = F,
                     na.strings = "")


fl_count_by_scaff$anc <- ifelse(fl_count_by_scaff$V1 %in% anc_MAGs$V1, "Y", "N")
west_cbs$anc <- ifelse(west_cbs$V1 %in% anc_MAGs$V1, "Y", "N")
fl_anc <- fl_count_by_scaff %>% filter(anc == "Y") %>% select(V1)
west_anc <- west_cbs %>% filter(anc == "Y") %>% select(V1)
fl_anc$fl <- "FL"
west_anc$west<- "West"




all_mags <- merge(anc_MAGs, fl_anc, all.x = T)
all_mags <- merge(all_mags, west_anc, all.x = T)
all_mags$anc <- "ANC"
all_mags$order <- rowSums(is.na(all_mags))


all_mags <- all_mags[order(all_mags$order),]
all_mags$xcoor <- 1:nrow(all_mags)
long_mags <- pivot_longer(all_mags,cols = V2:anc)

y_coord_key <- data.frame(cat = c("ANC", "S17463", "S17464", "S17465", "S17466", "FL", "West"), ycor = c(4.75, 4, 3.5, 3, 2.5, 1.75, 1))

testplot <- merge(long_mags, y_coord_key, by.x = "value", by.y = "cat")
dotplot <- ggplot() + geom_point(data = testplot[testplot$value %in% c("ANC", "West", "FL"),], aes(x = xcoor, y = ycor, col = value), size = 7) + ylim(c(0, 5)) + theme_void() +
  geom_point(data = testplot[testplot$value %in% c("S17463", "S17464", "S17465", "S17466"),], aes(x = xcoor, y = ycor, col = value), size = 4) + scale_color_manual(values=c("#832E21", "#254E70", "azure4", "azure4", "azure4", "azure4", "#6E9075")) + 
  theme(legend.position = "none")
ggsave("~/Box Sync/Projects/aDNA_final_docs/dotplot_mags.png", dotplot, width = 4, height = 2, units = "in", dpi = 300)
##GET ANCIENT DATA BY SAMP
#add a breakdown by geography from the westrich meta 
#it does hold, just alter plot... perhaps manually in Powerpoint for time

stat_test <- matrix(c(7, 3, 4, 6), byrow = T, nrow = 2, ncol = 2)
rownames(stat_test) <- c("LK", "MR")
colnames(stat_test) <- c("Present", "Absent")
result <- fisher.test(stat_test)

