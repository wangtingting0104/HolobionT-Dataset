####### Alpha diversity ######
alpha = read.table("3vine-t1-alpha.txt", header=T, row.names=1, sep="\t", comment.char="")
design =read.table("5vine1-group.txt", header=T, row.names=1, sep="\t", comment.char="")
#design=design[,-1]
design$groupID=design$group

# Add design to alpha
index = cbind(alpha, design[match(rownames(alpha), rownames(design)), ])

m = "shannon"
m= "richness"
m= "chao1"

shapiro.test(index[[m]])#N


bartlett.test(index[[m]] ~ groupID, data = index) #



model = aov(index[[m]] ~ groupID, data=index)##
summary(model)



source("boxplerk.R")##Kruskal-Wallis and wilcox test

a=boxplerk(index[[m]],index$groupID, 
           ylab = m,  xlab = "Compartment",
           bcol = "bisque",p.adj = "fdr",las = 1) ##


stat = a$comparison
index$stat=stat[as.character(index$groupID),]$group






library(agricolae)
a = LSD.test(model,"groupID", p.adj="fdr") # alternative fdr

  stat = a$groups
  index$stat=stat[as.character(index$groupID),]$groups



max=max(index[,c(m)])
min=min(index[,c(m)])
x = index[,c("groupID",m)]
library(dplyr)
y = x %>% group_by(groupID) %>% summarise_(Max=paste('max(',m,')',sep=""))
y=as.data.frame(y)
rownames(y)=y$groupID
index$y=y[as.character(index$groupID),]$Max + (max-min)*0.05

library(ggplot2)
p = ggplot(index, aes(x=groupID, y=index[[m]], color=groupID)) +
  scale_colour_manual(values=rev(c("C"="#CDC733","WC"="#E93B81","W"="#B6C9F0","CW"="#3EDBF0")))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste(m, "index")) + theme_classic() +
  geom_text(data=index, aes(x=groupID, y=y, color=groupID, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=3, alpha=1)+
  scale_x_discrete(limits = c("C", "WC", "W","CW" ))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line.x=element_line(size=.5, colour="black"),
        axis.line.y=element_line(size=.5, colour="black"),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=10),
        legend.position="right",
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text= element_text(size=10),
        text=element_text(size=12))



p



ggsave(paste("6vine1-CHAO1-PLOT-ORDER.pdf", sep=""), p, width = 4, height = 3)
####################################################################

