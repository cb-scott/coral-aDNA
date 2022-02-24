#!/usr/bin/env Rscript
#needed for remote pipeline, will run on command line 
args = commandArgs(trailingOnly=TRUE)

#first arg = our translator table
#second arg = our modern snp file (or a subset of it)
#third arg = output file 

translate <- read.table(args[1], sep = '\t', header = F)
translate$V5 <- as.numeric(translate$V5)
translate$V3 <- as.numeric(translate$V3)
translate$V4 <- as.numeric(translate$V4)

a <- read.table(args[2], sep="\t", header=F)

searcher <- function(df){
  chr <- as.numeric(df[2])
  #worked
  if(chr %in% 4:14){
    x <- cbind("ScaffCat1", NA)
  }
  if(chr %in% 15:27){
    x <- cbind("ScaffCat2", NA)
  }
  #print(sub$V3 <= val)
  #print(sub$V4 >= val)
  if(chr %in% 1:3){
    val <- as.numeric(df[4])
    sub <- translate[translate$V3 == chr,]  #would be faster to reference existing subsets rather than create new ones?
    x <- sub[sub$V4 <= val & sub$V4 >= val,]
    if(is.na(x[1,1])){
      print(val)
    }
  }

  #did not work 
  return(x[1,1])
} #this takes a long time on the whole file 
a$og_chr <- apply(a, 1, searcher)
write.table(a, args[3], sep = '\t', row.names = F, col.names = F, quote = F)
