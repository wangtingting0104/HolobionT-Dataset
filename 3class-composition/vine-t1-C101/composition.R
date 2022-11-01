data_test <- read.table('6vine C01-T1-CLASS.txt',sep='\t',h=T,row.names=1,check=F,comment='')
statistics = apply(data_test, 1, sum)   # 得到每个样本的观测值总和

plot(statistics)
data_percent = data.frame()   # 建立空数据框

# 每个值除以前面得到的总和获得占比

for (n in 1:4) {
  
  data_percent = rbind( data_percent, data_test[n,] / statistics[n] )  
  
} 



statistics = apply(data_percent, 1, sum)   

plot(statistics)




data_percent$names = c(LETTERS[seq( from = 1, to = 4)]#,
                       
                       #letters[seq( from = 1, to = 4 )]
                       )

library(reshape2)

data_plot = melt(data_percent)

colnames(data_plot) = c('names','variable','value')



#group = c( rep('Upper',15), rep('Lower',15))

#data_plot$group = rep(group,7)





library(ggplot2)



p = ggplot( data_plot, aes( x = names, weight = value, fill = variable))+
  
  geom_bar( position = "stack")  # 


p


order_x =  apply( data_percent[,1:17], 2, sum)   # 


order_x = order_x[order(order_x, decreasing = T)]  # decreasing = T


order_x  # 






data_plot$variable = factor(data_plot$variable,
                            
                            levels = names(order_x) ,  
                            
                            ordered = T )





p2 = ggplot( data_plot,aes(x = names, weight = value, fill = variable))+
  
  geom_bar(position = "stack")  

p2
p3<-p2 +   scale_fill_manual(values =rev(c("Acidobacteria_Gp1"="#BEBADA","Acidobacteria_Gp10"="#BEBADA","Actinobacteria"="#C490E4",
                                      "Alphaproteobacteria"="#FFFFB3","Bacilli"="#FB8072","Betaproteobacteria"="#80B1D3", 
                                      "Chlamdiia"="#B2B8A3","Clostridia"="#FF0000", "Cytophagia"="#B3DE69","Deltaproteobacteria"="#DBE6FD",
                                      "Flavobacteriia"="#FFBB78FF", "Gammaproteobacteria"="#FCCDE5", "Holophagae"="#0A1D37","Sacchar"="#B4846C","Sphingobacteriia"="#CCEBC5","Spirochaetia"="#2F5D62","Unsigned"="#80734C")))

p3

ggsave(paste0("8c101-t1-class.pdf"), p3, width=100, height=100, units="mm")
write.csv(data_percent,"8c101-t1-percent.csv")

