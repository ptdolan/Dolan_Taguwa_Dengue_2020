library(reshape2)
library(ggplot2)
library(limma)
library(ggpubr)
library(dplyr)
# #/Users/ptdolan/Research/CirSeq/Dengue/Shuhei_ExperimentalData
# ###############
# # Focus size
# ###############
# 
setwd("~/Documents/GitHub/Dolan_Taguwa_2020/  ")
focus=read.csv("PhenotypicData_2017/FocusSize_SmallN.txt",sep = "\t")
fm<-melt(focus)
Namevecs<-strsplit(as.character(fm$variable),split = "")

fm$Host<-unlist(lapply(Namevecs,function(X) paste(X[1:2],collapse = "")))
fm$Host<-factor(ifelse(fm$Host=="Hu",yes = "Human","Mosquito"))

fm$Set<-factor(unlist(lapply(Namevecs,function(X) X[3])))
fm$Passage<-factor(unlist(lapply(Namevecs,function(X) X[4])))

plaqueSize<-ggplot(fm)+theme_pubr(border = T)+geom_boxplot(aes(Passage,value,color=Host))+
  geom_jitter(aes(Passage,value,color=Host),width = .1)+
  geom_smooth(aes(Passage,value,group=Host,color=Host))+scale_y_log10()+
  facet_grid(facets = Host~Set,scales = "free_y")+scale_color_brewer(palette = "Set1",direction = 1)

ggsave(filename = "plaqueSize.pdf",plaqueSize, width=5, height=5)

# 
# ################
# # RNA production
# ################
# 
# 

RNA=read.csv("PhenotypicData_2017/RNAmeasurements.txt",sep = "\t")
mRNA<-melt(RNA,measure.vars = 7:10,variable.name = "TechRep")
mRNA$AdaptedHost<-factor(ifelse(mRNA$AdaptedHost=="C6/36","Mosquito","Human"))

ggsave(useDingbats=F,"RNAcontent.pdf",width = 4.8,height = 4,
  ggplot(na.omit(mRNA),mapping = aes(as.numeric(passage),value))+theme_pubr()+
  stat_summary(geom = "polygon",fun.data = mean_cl_boot,aes(fill=AdaptedHost))+
  geom_point(aes(color=AdaptedHost))+
  facet_grid(facets = TestedHost~set)+
  scale_colour_brewer(palette = "Set1",direction = 1)+
  scale_x_continuous(breaks = 1:9)+
  ylab("RNA content (AUs)")
)

meanRNA<-dcast(mRNA,formula = passage+AdaptedHost+set+bioRep+TechRep~TestedHost,function(X) mean(X,na.rm=T))
sdRNA<-dcast(mRNA,formula = passage+AdaptedHost+set~TestedHost,fun.aggregate = plotrix::std.error)

# 
# 
# ################
# # Viability
# ################

viability=read.csv("PhenotypicData_2017/VirusTiter_Viability.txt",sep = "\t")

mviability<-na.omit(melt(viability,measure.vars = 2:10,variable.name = "passage",value.name = "titer"))

castV<-dcast(mviability,Population+Rep+passage~CellLine)
castV<-data.frame(castV,relC636=castV[,4:8]/castV[,4])
castV<-data.frame(castV,relHuh7=castV[,4:8]/castV[,7])

mviability<-melt(castV,variable.name = "CellLine",value.name = "titer",id.vars = 1:3)
mviability$set<-limma::strsplit2(mviability$Population,split = " ")[,2]
mviability$Host<-limma::strsplit2(mviability$Population,split = " ")[,1]

rel<-mviability[limma::strsplit2(mviability$CellLine,"\\.")[,1]%in%c("relC636", "relHuh7"),]
mviability<-mviability[limma::strsplit2(mviability$CellLine,"\\.")[,1]%in%c("relC636", "relHuh7")==F,]

rel$comparison<-factor(limma::strsplit2(rel$CellLine,"\\.")[,1])
rel$testCell<-factor(limma::strsplit2(rel$CellLine,"rel.*\\.")[,2])

ggsave(useDingbats=F,"CellLineTiters_relativeToAdapted_Lines.pdf",width=5,height=5,ggplot(rel,mapping = aes(passage,titer))+
  geom_hline(yintercept = 1)+
  theme_pubr(legend = c("right"))+
  stat_summary(fun.data = mean_se,geom = "smooth",se=T,aes(group=comparison,color=comparison))+
  #geom_point(aes(color=comparison))+
  facet_grid(facets = Host+set~testCell)+
  scale_y_log10(breaks=c(0.1,1,10),labels=c("1:10","1:1","10:1"))+
  scale_colour_brewer(palette = "Set1",direction = 1))

meanrel<-dcast(rel,Population+Host+testCell+set+passage+comparison~.,value.var = "titer",fun.aggregate = function(C){mean(C,na.rm=T)})

ggsave(useDingbats=F,"CellLineTiters_relativeToAdapted_Heatmap.pdf", width=5,height=3.5,
ggplot(meanrel,mapping = aes(passage,testCell,fill=.))+
  theme_pubr(legend = c("right"))+
  geom_tile()+
  facet_grid(facets = Host+set~comparison)+
#  scale_y_log10(breaks=c(0.1,1,10),labels=c("1:10","1:1","10:1"))+
  scale_fill_gradient2(trans="log10",low = "cyan",mid = "black",high="magenta",limits=c(.02,50),midpoint = 0)
)

ggplot(na.omit(mviability))+
  theme_pubr(legend = c("right"))+
  stat_summary(fun.data = mean_se,geom = "pointrange",aes(passage,titer,group=CellLine,fill=CellLine,color=CellLine))+
  stat_summary(fun.data = mean_se,geom = "smooth",aes(passage,titer,group=CellLine))+
  #geom_point(aes(color=CellLine))+
  facet_grid(facets = ~Population)+
  scale_y_log10()+
  scale_colour_brewer(palette = "Set1",direction = -1)+
  scale_fill_brewer(palette = "Set1",direction = -1)

# #RELATIVE VIABILITY
relTiters<-dcast(rel,passage+Host+set+testCell+Rep~comparison,value.var = "titer")

ggplot(rel)+theme_pubr()+
  geom_hline(yintercept = 1)+
  stat_summary(fun.data = mean_se,geom='smooth',aes(passage,titer,group=paste(set,Host,testCell),color=testCell))+
  facet_grid(Host~comparison)+
  scale_y_log10()

mviability$Host<-factor(mviability$Host,levels = c("Huh","C6/36"))
castViability<-dcast(mviability,passage+Host+set~CellLine,fun.aggregate = mean,fun.args=c(na.rm=T),value.var = "titer")
castSEMViability<-dcast(mviability,passage+Host+set~CellLine,fun.aggregate = function(X){plotrix::std.error(X,na.rm=T)},value.var = "titer")

# #####################
# # Titers 2016 (N=?1)
# #####################

titers=read.csv("PhenotypicData_2017/titers_2016.txt",sep = "\t")
mtiters<-melt(titers,id.vars = 1,variable.name="sample",value.name = "titer")
mtiters$testHost<-strsplit2(mtiters$sample,"_")[,2]
mtiters$adaptedHost<-paste(strsplit2(strsplit2(mtiters$sample,"_")[,1],split = "")[,1],strsplit2(strsplit2(mtiters$sample,"_")[,1],split = "")[,2],sep="")
mtiters$adaptedHost<-ifelse(mtiters$adaptedHost=="Hu","Human","Mosquito")
mtiters$rep<-strsplit2(mtiters$sample,"")[,3]

ggplot(mtiters,aes(X,titer))+
  stat_summary(geom = 'smooth',aes(color=adaptedHost,group=sample),fun.data = mean_se)+
  geom_point(aes(pch = rep))+
  facet_grid(~testHost)+
  scale_colour_brewer(palette = "Set1",direction = 1)+
  scale_y_log10()+theme_pubr()

mt<-melt(mtiters,measure.vars = "titer")
mts<-dcast(mt,formula = X+adaptedHost+rep~variable+testHost,value.var = 'value')

# # #
# Make allPheno datatable with all data from plotting
# # #

allPheno<-data.table(castViability,mts[,-c(1:3)],RNA=meanRNA[-c(1:3)])
save(file = "allPheno.Rdata",allPheno)
allPheno[,PFU_RNA_C636:=C6.36/`RNA.C6/36`]
allPheno[,PFU_RNA_Huh7:=Huh7/`RNA.Huh7`]

ggplot(allPheno)+
  #stat_summary(na.rm = F,geom="smooth",fun.data=mean_ci,aes(passage,PFU_RNA_C636,group=set,col=Host,fill=Host))+
  geom_smooth(aes(passage,PFU_RNA_C636,group=set,col=Host))+
  geom_point(aes(passage,PFU_RNA_C636,group=set,col=Host))+facet_grid(~Host)+scale_y_log10()

ggplot(allPheno)+
  geom_smooth(aes(passage,Huh7/RNA.Huh7,group=set,col=Host))+
  geom_point(aes(passage,Huh7/RNA.Huh7,group=set,col=Host),alpha=0.5)+facet_grid(~Host)+scale_y_log10()
