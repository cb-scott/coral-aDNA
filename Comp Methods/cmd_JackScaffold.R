#!/usr/bin/env Rscript
#applied Jack Scaffold - you need a RENAMED SNP FILE FROM CHR_NAMER!!!
#Needed for remote pipeline

args = commandArgs(trailingOnly=TRUE)

modernsnp = args[1]
moderngeno = args[2]
modernindivids = args[3]
#args[4] = out prefix 

read_inds <- read.table(modernindivids, sep="\t", header=F)
nind = nrow(read_inds)

a <- read.table(modernsnp, sep = '\t', header = F)
nchroms <- length(unique(a$V7))

#b <- read.fwf(moderngeno,widths=c(1))
b <- readLines(moderngeno)
c <- cbind(a, b, stringsAsFactors=F)
#we're creating the same numbe of columns as chromosomes + 1, plus an additional columns of SNP data 
for (i in 1:nchroms){ #create as many columns as there are chromosomes 
  c <- cbind(c, b, stringsAsFactors=F)
}

replace_string <- paste(rep("9", nind), sep="", collapse="") #works
un_scaff <- unique(a$V7)
for (j in 1:nchroms)
{
  c[which(c[,7]==un_scaff[j]),7+(j+1)] <- replace_string #alter indexing 
}

c$x <- apply(c[ ,c(8:ncol(c))],1,paste, collapse = "")
write.table(c$x, paste(args[4], ".geno", sep=""), sep="\t", row.names=F, col.names=F, quote=F)

orig <- read.table(modernindivids, sep="\t", header=F)
orig$V1 <- as.character(orig$V1)
orig$V2 <- as.character(orig$V2)
orig$V3 <- as.character(orig$V3)
a <- orig
for (j in 1:nchroms){
  tobind <- cbind(paste(orig$V1,".",j,sep=""), "U", paste(orig$V3,".",j,sep=""))
  a <- rbind(a, tobind, stringsAsFactors=F)
}
write.table(a, paste(args[4], ".ind", sep=""), sep="\t", row.names=F, col.names=F, quote=F)


