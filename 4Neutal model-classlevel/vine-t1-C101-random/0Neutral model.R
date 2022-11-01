

library(dplyr)
source("Neutral model.r")

metacommunity = read.csv(file ="0vine-t1-c101-pool.csv" ,header=TRUE,sep=",",row.names = 1)
metacommunity[1:10,1:10]

metacommunity = t(metacommunity)


otu = read.csv(file ="0vine-t1-c101-otutable.csv" ,header=TRUE,sep=",",row.names = 1)
otu=t(otu)
design = read.csv(file ="2metadata.csv" ,header=TRUE,row.names=1,sep=",")
head(design)
###CFT1
aa=design[which(design$rotation=="X101CT1"),]

bb=aa$sample

comun=otu[bb,]


Neutral.fit(spp=comun, pool=metacommunity, stats=TRUE)
neu <- Neutral.fit(spp=comun, pool=metacommunity, stats=FALSE)
neu <- transform(neu, plog = log10(p))
neu$OTU<-rownames(neu)
up<-filter(neu,freq>pred.upr)
up$class<-"Abo"
lw<-filter(neu,freq<pred.lwr)
lw$class<-"Bel"
mid<-filter(neu,pred.lwr<=freq&freq<=pred.upr)
mid$class<-"Neu"
plot.data=merge(up,lw,all=T)
plot.data=merge(plot.data,mid,all=T)

write.table(plot.data,file="3X101CT1.csv",sep = ",",row.names=TRUE)

p1=ggplot(plot.data)+geom_point(aes(x=plog,y=freq,colour=class),alpha=0.8,size=2)+
  geom_line(aes(x=plog,y=freq.pred),colour="#1E90FF",size=1)+
  geom_line(aes(x=plog,y=pred.upr),colour="#1E90FF",size=1,linetype="dashed")+
  geom_line(aes(x=plog,y=pred.lwr),colour="#1E90FF",size=1,linetype="dashed")+
  labs(x="Mean relative abundance (log10)",y="Occurrence frequency")+
  scale_colour_manual(values=c("#FFD166","#06D6A0","#8D99AE"))+
  labs(title = "X101CT1")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+guides(colour=FALSE)
p1


###FCT1
aa=design[which(design$rotation=="C101T1"),]

bb=aa$sample

comun=otu[bb,]


Neutral.fit(spp=comun, pool=metacommunity, stats=TRUE)
neu <- Neutral.fit(spp=comun, pool=metacommunity, stats=FALSE)
neu <- transform(neu, plog = log10(p))
neu$OTU<-rownames(neu)
up<-filter(neu,freq>pred.upr)
up$class<-"Abo"
lw<-filter(neu,freq<pred.lwr)
lw$class<-"Bel"
mid<-filter(neu,pred.lwr<=freq&freq<=pred.upr)
mid$class<-"Neu"
plot.data=merge(up,lw,all=T)
plot.data=merge(plot.data,mid,all=T)

write.table(plot.data,file="3C101T1.csv",sep = ",",row.names=TRUE)

p2=ggplot(plot.data)+geom_point(aes(x=plog,y=freq,colour=class),alpha=0.8,size=2)+
  geom_line(aes(x=plog,y=freq.pred),colour="#1E90FF",size=1)+
  geom_line(aes(x=plog,y=pred.upr),colour="#1E90FF",size=1,linetype="dashed")+
  geom_line(aes(x=plog,y=pred.lwr),colour="#1E90FF",size=1,linetype="dashed")+
  labs(x="Mean relative abundance (log10)",y="Occurrence frequency")+
  scale_colour_manual(values=c("#FFD166","#06D6A0","#8D99AE"))+
  labs(title = "C101T1")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+guides(colour=FALSE)
p2



library(ggpubr)
pdf("4random.pdf",width=10,height=10)
ggarrange(p1,p2, ncol=2,nrow=2)
dev.off()

library(reshape2)
library(ggplot2)
library(vegan)
OTU = read.csv(file ="0vine-t1-c101-otutable.csv" ,header=TRUE,sep=",",row.names = 1)
OTU_t=t(OTU)
OTU_t[1:10,1:10]
design = read.csv(file ="2metadata.csv" ,header=TRUE,sep=",")

#3X101CT1
model = read.csv(file ="3X101CT1.csv",header=TRUE,row.names=1,sep=",")
head(model)


aa=design[which(design$rotation=="X101CT1"),]

bb=aa$sample
OTU1=OTU_t[bb,]
OTU1 =  as.data.frame(colMeans(OTU1))
OTU1=as.matrix(OTU1)
head(OTU1)

library("questionr")
a=cprop(OTU1)###
OTU2=as.data.frame(a[!rownames(a) %in% c("Total") , -which(colnames(a) %in% c("All"))])##删除求相对丰度生成的总行和列
head(OTU2)
OTU2$OTU=rownames(OTU2)
OTU3=merge(model[,c("OTU","class")],OTU2,by.x="OTU",all= T)
head(OTU3)
OTU4=aggregate(OTU3[,(-1:-2)],list(OTU3$class),sum)
head(OTU4)


label=c(rep("X101CT1", 3))
OTU5=data.frame(OTU4, label)

p1 = ggplot(OTU5, aes(x=label, y = x, fill = Group.1 )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("")+
  ylab("X101CT1-Cumulative relative abundance")+
  theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Abo"="#FFD166","Bel"="#06D6A0","Neu"="#8D99AE"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p1
write.csv(OTU5,"5X101-com.csv")
#3FCT1
model = read.csv(file ="3C101T1.csv",header=TRUE,row.names=1,sep=",")
head(model)


aa=design[which(design$rotation=="C101T1"),]

bb=aa$sample
OTU1=OTU_t[bb,]
OTU1 =  as.data.frame(colMeans(OTU1))
OTU1=as.matrix(OTU1)
head(OTU1)


library("questionr")
a=cprop(OTU1)###Column percentages##
OTU2=as.data.frame(a[!rownames(a) %in% c("Total") , -which(colnames(a) %in% c("All"))])#\
head(OTU2)
OTU2$OTU=rownames(OTU2)
OTU3=merge(model[,c("OTU","class")],OTU2,by.x="OTU",all= T)
head(OTU3)
OTU4=aggregate(OTU3[,(-1:-2)],list(OTU3$class),sum)
head(OTU4)


label=c(rep("C101T1", 3))
OTU5=data.frame(OTU4, label)

p2 = ggplot(OTU5, aes(x=label, y = x, fill = Group.1 )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("")+
  ylab("3C101T1-Cumulative relative abundance")+
  theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Abo"="#FFD166","Bel"="#06D6A0","Neu"="#8D99AE"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p2

write.csv(OTU5,"5c101-com.csv")

library(ggpubr)
pdf("5RANDOM-composition.pdf",width=6,height=10)
ggarrange(p1,p2, ncol=2,nrow=2)
dev.off()




OTU = read.csv(file ="0vine-t1-c101-otutable.csv" ,header=TRUE,sep=",",row.names = 1)
OTU[1:10,1:10]
OTU_t=t(OTU)
OTU_t[1:10,1:10]
design = read.csv(file ="2metadata.csv" ,header=TRUE,sep=",")



#X101CT1
aa=design[which(design$rotation=="X101CT1"),]

bb=aa$sample
OTU1=as.data.frame(OTU_t[bb,])
OTU1[1:10,1:10]
OTU2=as.data.frame(t(OTU1))
OTU2$OTU=rownames(OTU2)

model = read.csv(file ="3X101CT1.csv",header=TRUE,row.names=1,sep=",")

head(model)

OTU3=merge(model[,c("OTU","class")],OTU2,by.x="OTU",all= F)
OTU3=rename(OTU3,"#OTU"=OTU)
OTU3[1:10,1:10]
Abo=OTU3[which(OTU3$class=="Abo"),]
Bel=OTU3[which(OTU3$class=="Bel"),]
Neu=OTU3[which(OTU3$class=="Neu"),]

Abo=Abo[which(rowSums(Abo[,-1:-2]) > 0),]
Bel=Bel[which(rowSums(Bel[,-1:-2]) > 0),]
Neu=Neu[which(rowSums(Neu[,-1:-2]) > 0),]

Abo=Abo[,-2]
Bel=Bel[,-2]
Neu=Neu[,-2]


write.table(Abo,"6X101CT1_Abo.txt",sep ='\t',quote = FALSE,row.names = F,col.names = T)
write.table(Bel,"6X101CT1_Bel.txt",sep ='\t',quote = FALSE,row.names = F,col.names = T)
write.table(Neu,"6X101CT1_Neu.txt",sep ='\t',quote = FALSE,row.names = F,col.names = T)


#FCT1
aa=design[which(design$rotation=="C101T1"),]

bb=aa$sample
OTU1=as.data.frame(OTU_t[bb,])
OTU1[1:10,1:10]
OTU2=as.data.frame(t(OTU1))
OTU2$OTU=rownames(OTU2)

model = read.csv(file ="3C101T1.csv",header=TRUE,row.names=1,sep=",")

head(model)

OTU3=merge(model[,c("OTU","class")],OTU2,by.x="OTU",all= F)
OTU3=rename(OTU3,"#OTU"=OTU)
OTU3[1:10,1:10]
Abo=OTU3[which(OTU3$class=="Abo"),]
Bel=OTU3[which(OTU3$class=="Bel"),]
Neu=OTU3[which(OTU3$class=="Neu"),]

Abo=Abo[which(rowSums(Abo[,-1:-2]) > 0),]
Bel=Bel[which(rowSums(Bel[,-1:-2]) > 0),]
Neu=Neu[which(rowSums(Neu[,-1:-2]) > 0),]

Abo=Abo[,-2]
Bel=Bel[,-2]
Neu=Neu[,-2]

write.table(Abo,"6C101T1_Abo.txt",sep ='\t',quote = FALSE,row.names = F,col.names = T)
write.table(Bel,"6C101T1_Bel.txt",sep ='\t',quote = FALSE,row.names = F,col.names = T)
write.table(Neu,"6C101T1_Neu.txt",sep ='\t',quote = FALSE,row.names = F,col.names = T)





library(reshape2)
library(ggplot2)
library(vegan)
#CFT1-ABO
OTU <- read.table('6X101CT1_Abo.txt', sep='\t', header=T, comment.char='', check.names=F)


{annotation = read.csv(file ="1vine-annotation.csv" ,header=TRUE,sep=",",row.names = 1)
annotation$"#OTU"=rownames(annotation)

spe1=merge(annotation,OTU,by.x="#OTU",all=F)


Phylum=aggregate(spe1[,(-1:-7)],list(spe1$Class),sum)



#Phylum=Phylum[-which(Phylum$Group.1=="unknow"),]

rownames(Phylum)=Phylum$Group.1
Phylum=as.matrix(Phylum[,-1])



library("questionr")
a=cprop(Phylum)###Column percentages##

OTU2=as.data.frame.array(a[!rownames(a) %in% c("Total") , -which(colnames(a) %in% c("All"))])##

head(OTU2)

OTU3 <- OTU2[which(rownames(OTU2)!="unknow"),]##

# Decreased sort by abundance
OTU4 = as.data.frame(OTU3[(order(-rowSums(OTU3))), ])
# Calculate average relative abundance for each group
OTU5 =  as.data.frame(rowMeans(OTU4))
OTU5$Phylum=rownames(OTU5)
head(OTU5)
# Filter Top 9 , and other group into Low abundance
other = as.data.frame(sum(OTU5[10:dim(OTU5)[1],1 ]))
other$Phylum=rownames(other)
OTU6 = OTU5[1:(10 - 1), ]
library(dplyr)
OTU6=rename(OTU6,"Mean"='rowMeans(OTU4)')
other=rename(other,"Mean"=`sum(OTU5[10:dim(OTU5)[1], 1])`)

mean_sort = rbind(OTU6,other)
rownames(mean_sort)[10] = c("Low abundance")


bb=as.data.frame(rowMeans(OTU2[which(rownames(OTU2)=="unknow"),]))
bb$Phylum=rownames(bb)
bb=rename(bb,"Mean"=`rowMeans(OTU2[which(rownames(OTU2) == "unknow"), ])`)
mean_sort=rbind(mean_sort,bb)##


# data melt for ggplot2
mean_sort$Phylum = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("Phylum")))
# Set taxonomy order by abundance, default by alphabet##

data_all$Phylum  = factor(data_all$Phylum, levels=rownames(mean_sort))

tax=data_all[,1]
tax=as.data.frame(tax)
rownames(tax)<-rownames(data_all)
tax$value<-as.numeric(data_all$value)
tax$variable<-as.character(data_all$variable)
data_all=tax}

###
p1 = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.8)+ 
  scale_y_continuous(labels = scales::percent) + 

  xlab("6X101CT1_Abo")+
  
  ylab("Relative Abundance(%)")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Acetothermia"="#3C3A8D","Acidobacteria"="#BEBADA","Acidobacteria_Gp1"="#BEBADA","Actinobacteria"="#ff6eb4",
                                  "Armatimonadetes"="#000000","Bacteroidetes"="#44A8DB","BRC1"="#57B78C",
                                  "candidate_division_WPS-1"="#EA711B","candidate_division_WPS-2"="#CECE05",
                                  "Candidatus_Saccharibacteria"="#D72226","Chlamydiia"="#AA9953","Chloroflexi"="#093A3A",
                                  "Cyanobacteria"="#02759E","Deinococcus-Thermus"="#3A398C","Firmicutes"="#74B662",
                                  "Gemmatimonadetes"="#FFBB78FF","Latescibacteria"="#B53C0D", "Microgenomates"="#C49C94FF",
                                  "Nitrospirae"="#03F1F7","Planctomycetes"="#F962C3","Proteobacteria"="#F95A5A",
                                  "Tenericutes"="#0808B7","Verrucomicrobia"="#144404","Low abundance"="#cccccc",
                                  "unknow"="#9E7C61","Rokubacteria"="#D3C13C","Alphaproteobacteria"="#ffec8b","Bacilli"="#89d0f5","Betaproteobacteria"="#ff9900","Cytophagia"="#66cc00","Deltaproteobacteria"="#666600","Gammaproteobacteria"="#fccde5","Flavobacteriia"="#ff0000","Sphingobacteriia"="#9f79ee","Unsigned"="#bebebe"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())

p1



####6C101T1_Abo
OTU <- read.table('6C101T1_Abo.txt', sep='\t', header=T, comment.char='', check.names=F)
{annotation = read.csv(file ="1vine-annotation.csv" ,header=TRUE,sep=",",row.names = 1)
  annotation$"#OTU"=rownames(annotation)
  
  spe1=merge(annotation,OTU,by.x="#OTU",all=F)
  
  
  Phylum=aggregate(spe1[,(-1:-7)],list(spe1$Class),sum)
  
  

    #Phylum=Phylum[-which(Phylum$Group.1=="unknow"),]
  
  rownames(Phylum)=Phylum$Group.1
  Phylum=as.matrix(Phylum[,-1])
  
  

    library("questionr")
  a=cprop(Phylum)###Column percentages##
  
  OTU2=as.data.frame.array(a[!rownames(a) %in% c("Total") , -which(colnames(a) %in% c("All"))])##
  
  
  OTU3 <- OTU2[which(rownames(OTU2)!="unknow"),]##
  
  # Decreased sort by abundance
  OTU4 = as.data.frame(OTU3[(order(-rowSums(OTU3))), ])
  # Calculate average relative abundance for each group
  OTU5 =  as.data.frame(rowMeans(OTU4))
  OTU5$Phylum=rownames(OTU5)
  head(OTU5)
  # Filter Top 9 , and other group into Low abundance
  other = as.data.frame(sum(OTU5[10:dim(OTU5)[1],1 ]))
  other$Phylum=rownames(other)
  OTU6 = OTU5[1:(10 - 1), ]
  library(dplyr)
  OTU6=rename(OTU6,"Mean"='rowMeans(OTU4)')
  other=rename(other,"Mean"=`sum(OTU5[10:dim(OTU5)[1], 1])`)
  
  mean_sort = rbind(OTU6,other)
  rownames(mean_sort)[10] = c("Low abundance")
  

    bb=as.data.frame(rowMeans(OTU2[which(rownames(OTU2)=="unknow"),]))
  bb$Phylum=rownames(bb)
  bb=rename(bb,"Mean"=`rowMeans(OTU2[which(rownames(OTU2) == "unknow"), ])`)
  mean_sort=rbind(mean_sort,bb)##
  
  
  # data melt for ggplot2
  mean_sort$Phylum = rownames(mean_sort)
  data_all = as.data.frame(melt(mean_sort, id.vars=c("Phylum")))
  # Set taxonomy order by abundance, default by alphabet##y
  
  data_all$Phylum  = factor(data_all$Phylum, levels=rownames(mean_sort))
  
  tax=data_all[,1]
  tax=as.data.frame(tax)
  rownames(tax)<-rownames(data_all)
  tax$value<-as.numeric(data_all$value)
  tax$variable<-as.character(data_all$variable)
  data_all=tax}
p4 = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.8)+ 
  scale_y_continuous(labels = scales::percent) + 
  
  xlab("6C101T1_Abo")+
  
  ylab("Relative Abundance(%)")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Acetothermia"="#3C3A8D","Acidobacteria"="#BEBADA","Acidobacteria_Gp1"="#BEBADA","Actinobacteria"="#ff6eb4",
                                  "Armatimonadetes"="#000000","Bacteroidetes"="#44A8DB","BRC1"="#57B78C",
                                  "candidate_division_WPS-1"="#EA711B","candidate_division_WPS-2"="#CECE05",
                                  "Candidatus_Saccharibacteria"="#D72226","Chlamydiae"="#AA9953","Chloroflexi"="#093A3A",
                                  "Cyanobacteria"="#02759E","Deinococcus-Thermus"="#3A398C","Firmicutes"="#74B662",
                                  "Gemmatimonadetes"="#FFBB78FF","Latescibacteria"="#B53C0D", "Microgenomates"="#C49C94FF",
                                  "Nitrospirae"="#03F1F7","Planctomycetes"="#F962C3","Proteobacteria"="#F95A5A",
                                  "Tenericutes"="#0808B7","Verrucomicrobia"="#144404","Low abundance"="#cccccc",
                                  "unknow"="#9E7C61","Rokubacteria"="#D3C13C","Alphaproteobacteria"="#ffec8b","Bacilli"="#89d0f5","Betaproteobacteria"="#ff9900","Cytophagia"="#66cc00","Deltaproteobacteria"="#666600","Gammaproteobacteria"="#fccde5","Flavobacteriia"="#ff0000","Sphingobacteriia"="#9f79ee","Unsigned"="#bebebe"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())

p4




library(ggpubr)
pdf("7TIunknow.pdf",width=12,height=12)
ggarrange(p1,p4, ncol=2,nrow=2)
#,p4,p5,p6,p7,p8,p9,
dev.off()
##


