library(ggplot2)
library(Hmisc)
library(data.table)
library(ggpubr)
library(Biostrings)
library(RColorBrewer)

Host.colors<- brewer.pal(name = "Set1",4)[c(1,2)]
status.colors<- c(
  brewer.pal(name = "Greys",9)[5],
  brewer.pal(name = "Purples",9)[7],
  "#000000",
  brewer.pal(name = "YlOrRd",9)[4]
)

readFEs<-function(fname){
  fitnessEstimates=fread(fname,header = T,stringsAsFactors = T,col.names = c("R","W","name"))
  fitnessEstimates[,host:=as.factor(limma::strsplit2(name,"_")[,1])]
  fitnessEstimates[,set:=as.factor(limma::strsplit2(name,"_")[,2])]
  fitnessEstimates[,passage:=as.factor(limma::strsplit2(name,"_")[,3])]
  #print(fitnessEstimates)
  return(fitnessEstimates)
}

plotFEs<-function(fE,name="meanFitness_MaskedGenomes.pdf"){
  fE<-fitnessEstimates[R<1000]
  WP<-ggplot(fE)+theme_bw()
  WPsmooth<-WP+
    stat_summary(aes(passage,y=W,color=.id),geom = "line",fun.y= median)+
    facet_grid(.id~host+set,scales = "free_y")+
    scale_color_manual(values = c("blue",status.colors[c(4,2,3,1)]))+
    scale_x_continuous(breaks = 1:9)
  
  WPsmooth<-WP+
    stat_summary(aes(passage,y=W,color=.id),geom = "line",fun.data=mean_cl_boot)+
    facet_grid(.id~host+set,scales = "free_y")+
    scale_color_manual(values = c("blue",status.colors[c(4,2,3,1)]))+
    scale_x_continuous(breaks = 1:9)
  
  WPribbon<-WP+
    stat_summary(aes(passage,ymin=1,ymax=W,y=W,fill=.id),geom = "ribbon",fun.ymax= mean)+
    facet_grid(.id~host+set,scales = "free_y")+
    scale_fill_manual(values = c("blue",status.colors[c(4,2,3,1)]))+
    scale_x_continuous(breaks = 1:9)
  
  ggsave(WPsmooth,filename = name,width=6,height=4)
  ggsave(WPribbon,filename = paste("Ribbon_",name,sep=""),width=6,height=4)
}


#Analysis by Class counts
analyzeGenomes<-function(sampledGenomes,outfile="MutationsPerGenome.pdf"){
  sampledGenomes[is.na(sampledGenomes)]=0
  
  sampledGenomes[,host:=limma::strsplit2(name,"_")[,1]]
  sampledGenomes[,set:=limma::strsplit2(name,"_")[,2]]
  sampledGenomes[,passage:=limma::strsplit2(name,"_")[,3]]
  
  mSG<-melt(sampledGenomes,id.vars = c("V1","name","set","host","passage"),variable.name = "Class",value.name = "Count")
  
  SG<-ggplot(mSG[Class!="WT",])+theme_pubr()+
    scale_color_manual(values=status.colors[c(4,2,3,1)])+
    scale_fill_manual(values=status.colors[c(4,2,3,1)])
  
  SGhist<-SG+geom_histogram(binwidth = 1,aes(Count, fill=Class))+facet_grid(passage~Class)
  
  ClassSmooth<-SG+
    stat_summary(fun.data= mean_se,geom="area",aes(x=passage,y = Count,group=reorder(,Class),color=Class,fill=Class))+
    facet_grid(host~set)
  
  ggsave(ClassSmooth,filename = paste("ClassCount_",outfile,sep=""))
  
  LethalSmooth<-SG+
    stat_summary(data=mSG[Class=="L"],fun.data= mean_cl_boot,geom="smooth",aes(x=passage,y = Count,group=Class,color=Class,fill=Class))+
    facet_grid(host~set)
  
  ggsave(LethalSmooth,filename = paste("lethal_",outfile,sep=""))
}

##########
sampledGenomes=    fread("/Users/ptdolan/RSandbox/sampledGenomes.csv")
sampledGenomesOPP= fread("/Users/ptdolan/RSandbox/OppsampledGenomes.csv")
sampledGenomesSets=fread("/Users/ptdolan/RSandbox/SetssampledGenomes.csv")
# 
# analyzeGenomes(sampledGenomes,"MutationsPerGenome.pdf")
# analyzeGenomes(sampledGenomesOPP,"MutationsPerGenome_Opp.pdf")
# analyzeGenomes(sampledGenomesSets,"MutationsPerGenome_SetCtrl.pdf")

#Single Class fitness
fitnessEstimatesAll<-readFEs("/Users/ptdolan/Google_Drive/DengueDataAndAnalysis/DenguePanelsAndOutput_submissionDraft/NoneMaskedSamples.csv")
fitnessEstimatesB<-readFEs("/Users/ptdolan/Google_Drive/DengueDataAndAnalysis/DenguePanelsAndOutput_submissionDraft/BMaskedSamples.csv")
fitnessEstimatesL<-readFEs("/Users/ptdolan/Google_Drive/DengueDataAndAnalysis/DenguePanelsAndOutput_submissionDraft/LMaskedSamples.csv")
fitnessEstimatesN<-readFEs("/Users/ptdolan/Google_Drive/DengueDataAndAnalysis/DenguePanelsAndOutput_submissionDraft/NMaskedSamples.csv")
fitnessEstimatesD<-readFEs("/Users/ptdolan/Google_Drive/DengueDataAndAnalysis/DenguePanelsAndOutput_submissionDraft/DMaskedSamples.csv")

fitnessEstimates<-rbindlist(use.names = T,list(All=fitnessEstimatesAll,B=fitnessEstimatesB,N=fitnessEstimatesN,L=fitnessEstimatesL,D=fitnessEstimatesD),idcol = T)

fitnessEstimates[,.id:=as.factor(.id)]
fitnessEstimates[,passage:=as.integer(passage)]
plotFEs(fitnessEstimates,name = "meanFitness_DropOthers_SimulatedGenomes.pdf")

computes<-fitnessEstimates[,list(mn=mean(W),med=median(W)),by=list(passage,host,set,.id)]

WP<-ggplot(computes[!.id%in%c("All","N")])+geom_line(lwd=2,aes(passage,mn,color=.id),position = "identity")+
  facet_grid(set~host,scales = "free_y")+
  scale_color_manual(values = c(status.colors[c(4,2,3,1)]))+
  scale_x_continuous(breaks = 1:9)+scale_y_continuous(trans='log2')+ggpubr::theme_pubr()

ggsave(WP,filename = "~/Dropbox/DENV Cirseq MS/PDF_VERSION/MeanContribution.pdf",width=6 ,height=6)

WPsmooth<-
  ggplot(fitnessEstimates[R<2000&.id!="N",])+theme_bw()+facet_grid(set~host,scales = "free_y")+
  stat_summary(aes(passage,y=W,fill=.id),alpha=.5,geom = "ribbon",fun.data=mean_cl_boot)+
  stat_summary(fill="grey",alpha=0.2,data=fitnessEstimates[R<2000&.id=="All"&.id!="N",],aes(passage,y=W,fill=.id),geom = "area",fun.data=mean_cl_boot)+
  stat_summary(aes(passage,y=W,col=.id),alpha=.5,geom = "smooth",fun.data=mean_cl_boot)+
  scale_color_manual("Fitness Class",values = c(status.colors[c(1,4,2,3)]))+
  scale_fill_manual("Fitness Class", values = c(status.colors[c(1,4,2,3)]))+
  scale_x_continuous(breaks = 1:9)+scale_y_continuous(trans='log2')+ggpubr::theme_pubr()
  
ggsave(WPsmooth,filename="/Users/ptdolan/Google_Drive/DenguePanelsAndOutput_submissionDraft/DropOthers_GenomeSamples.pdf",width=6,height=6)


analyzeGenomes(fitnessEstimatesAll)
#DROP ONE
fitnessEstimatesAll<-readFEs("/Users/ptdolan/Google_Drive/DenguePanelsAndOutput_submissionDraft/NoneDropOne_MaskedSamples.csv")
fitnessEstimatesB<-readFEs("/Users/ptdolan/Google_Drive/DenguePanelsAndOutput_submissionDraft/BDropOne_MaskedSamples.csv")
fitnessEstimatesL<-readFEs("/Users/ptdolan/Google_Drive/DenguePanelsAndOutput_submissionDraft/LDropOne_MaskedSamples.csv")
fitnessEstimatesN<-readFEs("/Users/ptdolan/Google_Drive/DenguePanelsAndOutput_submissionDraft/NDropOne_MaskedSamples.csv")
fitnessEstimatesD<-readFEs("/Users/ptdolan/Google_Drive/DenguePanelsAndOutput_submissionDraft/DDropOne_MaskedSamples.csv")

fitnessEstimates<-rbindlist(use.names = T,
                            list(All=fitnessEstimatesAll,
                                 B=fitnessEstimatesB,
                                 N=fitnessEstimatesN,
                                 L=fitnessEstimatesL,
                                 D=fitnessEstimatesD),idcol = T)

fitnessEstimates[,.id:=as.factor(.id)]
fitnessEstimates[,passage:=as.integer(passage)]
plotFEs(fitnessEstimates,name = "meanFitness_DropOne_SimulatedGenomes.pdf")

WPsmooth<-
  ggplot(fitnessEstimates[R<1000&.id!="N",])+theme_bw()+
  stat_summary(aes(passage,y=W,color=.id),lwd=1.4,geom = "line",fun.y= function(X){mean(X)})+
  facet_grid(host~set,scales = "free_y")+
  scale_color_manual(values = c("lightgrey",status.colors[c(4,2,3,1)]))+
  scale_x_continuous(breaks = 1:9)+scale_y_log10(breaks=c(0.25,.5,1,2,4))

ggsave(WPsmooth,filename="/Users/ptdolan/Google_Drive/DenguePanelsAndOutput_submissionDraft/DropOne_GenomeSamples.pdf")

fitnessEstimates=readFEs("~/RSandbox/SetsFitnessSamples.csv")
plotFEs(fitnessEstimates,"meanFitness_SimulatedGenomes_Sets.pdf")

fitnessEstimates=readFEs("~/RSandbox/OppsFitnessSamples.csv")
plotFEs(fitnessEstimates,"meanFitness_SimulatedGenomes_Opps.pdf")
