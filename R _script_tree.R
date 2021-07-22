library(ggtree)
library(tidytree)
library(phangorn)
library(tidyr)
library(dplyr)
library(tidytree)
library(ggplot2)
library(treeio)
library(phangorn)
library(aplot)
library(cowplot)
library(xlsx)

setwd("$Datadir")

#### 1 SNP tree ####

snp <- read.csv("$Datadir/snp.fixrel.tsv",header = T, sep = "\t",stringsAsFactors = F)

snp <- subset(snp, snp$ALT_FREQ>=0.03)

str(snp)

nCoVsnp<- subset(snp, select=c("label","sample","POS","ALT_FREQ"))

colnames(nCoVsnp)<- c("label","sample","pos","freq")

class(nCoVsnp$pos)

nCoVsnp$pos <- as.numeric(nCoVsnp$pos)

class(nCoVsnp$freq)

nCoVsnp$label <- factor(nCoVsnp$label)


#Convert the alter frequency into a categorical variable
attach(nCoVsnp)

nCoVsnp$group[freq>=0.03&freq<0.1]<-"0"
nCoVsnp$group[freq>=0.1&freq<0.5]<-"1"
nCoVsnp$group[freq>=0.5&freq<0.9]<-"2"
nCoVsnp$group[freq>=0.9]<-"3"

detach(nCoVsnp)

class(nCoVsnp$group)

nCoVsnp$group<-factor(nCoVsnp$group,levels=c("0","1","2","3"),labels=c("3-10","10-50","50-90","≥90"),ordered=T)


#####1-1 Tree Visualization####
CoVtree <- read.newick("$Datadir/0706ncov.nwk")

CoVtree[["tip.label"]]

CoVtree2 <- root(CoVtree, outgroup="Illumina_5137_GZ_2021/5/21", edgelabel = TRUE)

nCoVsnp1 <- select(nCoVsnp, label)

p<- ggtree(CoVtree2,ladderize = T,size=0.5)

p

p1 <- p %<+% nCoVsnp1

p1[["data"]]

p1[["data"]][["x"]]<-30000*p1[["data"]][["x"]]

a<-as.data.frame(p1[["data"]])

p2<-p1+

  theme_tree2(legend.position=c(0.8,0.8))+geom_rootedge(0.08)+ylim2(p)

p2

p2[["data"]]

b<-as.data.frame(p2[["data"]])


#####1-2 Panel:SNP#####

nCoVsnp2 <- select(nCoVsnp,label,pos,freq, group)

d<-filter(p,isTip)%>%select(c(label,y))

nCoVsnp3<-left_join(nCoVsnp2,d,by="label")

p3<-ggplot(nCoVsnp3,aes(x=pos,y=y,color=group))+geom_point(size=2)+
  
  scale_color_manual("Frequency(%)",values=c("#FDBE85","#FD8D3C","#E6550D","#A63603"),na.value=NA)+
  
  labs(x="",y ='')+theme_classic()+ylim2(p)+theme(legend.position=c(0.2,0.8))

p3

p3[["data"]]

c<-as.data.frame(p3[["data"]])

p4<-plot_grid(p2,p3,ncol=2,align="h",rel_widths = c(1,2))

p4

ggsave("tree-SNP.pdf",p4,width=13,height=8)
