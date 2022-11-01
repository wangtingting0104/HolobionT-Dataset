otu <- read.csv("2VINE-t1rare.txt", head=TRUE,sep="\t",row.names = 1)
design <- read.table("3GROUP.txt",sep = "\t",header = T,colClasses = c("character"))





library(vegan)
library(ape)

#groups <- as.list(groups)
data <- t(otu)
data[is.na(data)] <- 0

data=decostand(data,method= "hellinger")

bray.dist = vegdist(data, method="bray")
write.table(as.matrix(bray.dist), file = "4bray.txt", sep="\t", col.names = NA)
#
spe=read.table("4bray.txt",header=T,row.names=1)
library(vegan)
adonis(spe~root*scion,data =design,permutations = 999,method="bray")


