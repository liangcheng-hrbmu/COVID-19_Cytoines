library(Seurat)
library(ggplot2)
library(stringr)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(MAST)
library(ggpubr)
library(monocle)
library(ggsci)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(SingleR)
library(harmony)
library(tidyverse)
library(ggplot2)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(DOSE)
library(limma)




result_path<-"path"

COVID_CCAresult<-readRDS(paste0(result_path,COVID_CCAresult.rds))

####删除线粒体大于20%的细胞####
COVID_CCAresult<-subset(COVID_CCAresult,percent.mt <20)

####删除低质量样本####

##查看单个样本的细胞数目##
table(COVID_CCAresult$orig.ident)

##mid1,S4,S12,S16,S17样本的细胞数目不足500为低质量样本##
mid1<-grep("*COVID.19.Mild.1",colnames(COVID_CCAresult))

S4<-grep("*COVID.19.Severe.4",colnames(COVID_CCAresult))
S12<-grep("*COVID.19.Severe.12",colnames(COVID_CCAresult))
S16<-grep("*COVID.19.Severe.16",colnames(COVID_CCAresult))
S17<-grep("*COVID.19.Severe.17",colnames(COVID_CCAresult))
##删除低质量样本##
all_de<-c(mid1,S4,S12,S16,S17)
COVID_CCAresult2<-COVID_CCAresult[,-all_de]
locmono<-which(COVID_CCAresult2$celltype%in%"Monocytes")
COVID_CCAresult2$celltype[locmono]<-"myeloid cells"
saveRDS(COVID_CCAresult2,paste0(result_path,"COVID_CCAresult2_del.rds"))

####F2B####
COVID_CCAresult2<-readRDS(paste0(result_path,"COVID_CCAresult2_del.rds"))
COVID_CCAresult2$celltype<-factor(COVID_CCAresult2$celltype,
                                  levels = c("Neutrophils","T cells","myeloid cells",
                                             "Mono-like cells","B cells","Erythrocytes"))


fig2B<-DimPlot(COVID_CCAresult2, reduction = "tsne",group.by = "celltype", label=F,cols = pal_nejm()(6))

####F2C####
COVID_CCAresult2 = SetIdent( COVID_CCAresult2, value = "celltype")
cluster_marker_data_allcells <- FindAllMarkers(COVID_CCAresult2,only.pos = T, slot = "data",test.use = "MAST",
                                      assay = "SCT",logfc.threshold = 0.25,min.pct = 0.25)

marker_gene<-cluster_marker_data_allcells %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
marker<-unique(marker_gene$gene)

marker<-c(marker_gene$gene[1:10],
          marker_gene$gene[31:40],
          marker_gene$gene[51:60],
          marker_gene$gene[41:50],
          marker_gene$gene[11:20],
          marker_gene$gene[21:30])

marker<-unique(marker)

##保存marker基因##
write.table(cluster_marker_data_allcells,paste0(result_path,"result2/Allcells_markergene.txt"),sep="\t",row.names=F,quote=F)
cluster_marker_data_allcells<-read.table(paste0(result_path,"result2/Allcells_markergene.txt"),sep="\t",header = T)

COVID_CCAresult2$celltype<-factor(COVID_CCAresult2$celltype,
                                  levels = c("Neutrophils","Mono-like cells","Erythrocytes",
                                             "B cells","T cells","myeloid cells"))

fig2C<-DotPlot(COVID_CCAresult2, features = marker,group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+  
  scale_color_gradientn(colours = 
                          c("#584B97","#4266A5","#3584B0","#4FA2A7","#6BBD9E",
                            "#8CCA9E","#AFD69E","#CCE199","#E4EBA0","#F0F2B4",
                            "#F7F1AD","#F8E093","#F8E093","#F6CA7B","#F4B167",
                            "#EE9154","#E96F46","#DA5646","#CA3F4B","#B02045",
                            "#951D41"))+theme(legend.key.size = unit(5, "pt"),
                                              legend.title =element_text( size = 8),
                                              axis.text.x = element_text(size = 8,angle = 90))

####F2D####
COVID_CCAresult2$celltype<-factor(COVID_CCAresult2$celltype,
                                  levels = c("Neutrophils","T cells","myeloid cells",
                                             "Mono-like cells","B cells","Erythrocytes"))
cell_numtable<-as.data.frame(table(COVID_CCAresult2$celltype, COVID_CCAresult2$orig.ident))
colnames(cell_numtable)<-c("Celltype","Sample","Number")
cell_numtable$Sample<-factor(cell_numtable$Sample,
                             levels = c("Healthy.1","Healthy.6","Healthy.5",
                                        "Healthy.2","Healthy.4","Healthy.3",
                                        "COVID.19.Mild.6","COVID.19.Mild.7",
                                        "COVID.19.Mild.4","COVID.19.Mild.5",
                                        "COVID.19.Mild.2","COVID.19.Mild.3",
                                        "COVID.19.Severe.5","COVID.19.Severe.20",
                                        "COVID.19.Severe.9","COVID.19.Severe.6",
                                        "COVID.19.Severe.14","COVID.19.Severe.1",
                                        "COVID.19.Severe.10","COVID.19.Severe.11",
                                        "COVID.19.Severe.13","COVID.19.Severe.15",
                                        "COVID.19.Severe.18","COVID.19.Severe.19",
                                        "COVID.19.Severe.2","COVID.19.Severe.3",
                                        "COVID.19.Severe.7","COVID.19.Severe.8" )        )



##单个样本细胞比例##
p5<-ggplot(cell_numtable,aes(x=Sample,y=Number,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+scale_fill_nejm()+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+  
  guides(fill=guide_legend(title=NULL))


cell_numtable2<-as.data.frame(table(COVID_CCAresult2$celltype, COVID_CCAresult2$Status))
colnames(cell_numtable2)<-c("Celltype","Sample","Number")
cell_numtable2$Sample <- factor(cell_numtable2$Sample,levels = c("Healthy", "COVID.19.Mid", "COVID.19.Severe"))

##细胞比例## ##scale_fill_manual(values=c(pal_nejm()(7)))
p6<-ggplot(cell_numtable2,aes(x=Sample,y=Number,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+ scale_fill_nejm()+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank())+  
  guides(fill=guide_legend(title=NULL))


fig2D<-ggarrange(p5,p6,
                 common.legend = T,legend = "right",
                 ncol = 2, nrow = 1)

####F2E####
CCL2_allcells<-FeaturePlot(COVID_CCAresult2, reduction = "tsne",slot ="scale.data",features = c("CCL2"),label = F,label.size = 3)+
  scale_color_gradientn(colours = 
                          c("#F8F6BF","#FFDCA0","#FAB880","#F69B6F","#F47F62",
                            "#F06C63","#E55065","#CD4271","#BF3B77","#AC327B",
                            "#9E2D7E","#7C2980","#4F2977","#2C1D57","#191135"))
CCL3_allcells<-FeaturePlot(COVID_CCAresult2, reduction = "tsne",slot ="scale.data",features = c("CCL3"),label = F,label.size = 3)+
  scale_color_gradientn(colours = 
                          c("#F8F6BF","#FFDCA0","#FAB880","#F69B6F","#F47F62",
                            "#F06C63","#E55065","#CD4271","#BF3B77","#AC327B",
                            "#9E2D7E","#7C2980","#4F2977","#2C1D57","#191135"))
CCL7_allcells<-FeaturePlot(COVID_CCAresult2, reduction = "tsne",slot ="scale.data",features = c("CCL7"),label = F,label.size = 3)+
  scale_color_gradientn(colours = 
                          c("#F8F6BF","#FFDCA0","#FAB880","#F69B6F","#F47F62",
                            "#F06C63","#E55065","#CD4271","#BF3B77","#AC327B",
                            "#9E2D7E","#7C2980","#4F2977","#2C1D57","#191135"))

fig2E<-ggarrange(CCL2_allcells,CCL3_allcells,CCL7_allcells,
                 common.legend = T,legend = "right",
                 ncol = 3, nrow = 1)


####F2F####
fig2F<-DotPlot(COVID_CCAresult2, features = c("CCL2","CCL3","CCL7"),group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+  
  scale_color_gradientn(colours = 
                          c("#584B97","#4266A5","#3584B0","#4FA2A7","#6BBD9E",
                            "#8CCA9E","#AFD69E","#CCE199","#E4EBA0","#F0F2B4",
                            "#F7F1AD","#F8E093","#F8E093","#F6CA7B","#F4B167",
                            "#EE9154","#E96F46","#DA5646","#CA3F4B","#B02045",
                            "#951D41"))+theme(legend.key.size = unit(5, "pt"),
                                              legend.title =element_text( size = 8),
                                              axis.text.x = element_text(size = 8,angle = 90))


####F2G####
##CCL3在所有细胞中的表达情况##
geneloc<-COVID_CCAresult2[which(rownames(COVID_CCAresult2)%in%"CCL3"),]
orig<-unique(geneloc$orig.ident)
result_exp<-c()
for(i in 1:length(orig)){
  cell<-which(geneloc$orig.ident%in%orig[i])
  exp0length<-length(geneloc@assays$RNA@counts[cell])  
  exp<-sum(geneloc@assays$RNA@counts[cell])/exp0length
  exp_orig<-cbind(orig[i],exp)
  result_exp<-rbind(result_exp,exp_orig)
}
colnames(result_exp)<-c("Sample","Exp")
result_exp<-as.data.frame(result_exp)
result_exp[["Status"]]<-"NA"
h<-grep("*Healthy.",result_exp$Sample)
result_exp[h,]$Status<-"Healthy"
S<-grep("*COVID.19.Severe.",result_exp$Sample)
result_exp[S,]$Status<-"COVID.19.Severe"
M<-grep("*COVID.19.Mild.",result_exp$Sample)
result_exp[M,]$Status<-"COVID.19.Mid"

result_exp$Status <- factor(result_exp$Status,levels = c("Healthy", "COVID.19.Mid", "COVID.19.Severe"))
result_exp<-as.data.frame(result_exp)
result_exp$Exp<-as.numeric(result_exp$Exp)

my_comparisons=list(c("Healthy","COVID.19.Mid"),c("COVID.19.Mid","COVID.19.Severe"),c("Healthy","COVID.19.Severe"))

CCL3_allcells_plot<-ggplot(result_exp,aes(x=Status,y=Exp,fill=Status))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,fill=c("#98C3E5","#88B382","#FDCB86"),outlier.color = "white")+
  geom_jitter(aes(fill=Status),width=0.2,shape=21,size=3.5)+
  scale_fill_manual(values = c("#005393","#18641E","#FF9300"))+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in all cells")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))

##将髓系细胞挑出##
myeloid_CCAresult<-subset(COVID_CCAresult2,celltype=="myeloid cells")
##CCL3在髓系细胞中的表达情况##
celltable <- function(object,str,gene){
  ####基因在细胞中的表达情况####
  geneloc<-object[which(rownames(object)%in%gene),]
  orig<-unique(geneloc$orig.ident)
  result_exp<-c()
  for(i in 1:length(orig)){
    cell<-which(geneloc$orig.ident%in%orig[i])
    exp0length<-length(geneloc@assays$RNA@counts[cell])  
    exp<-sum(geneloc@assays$RNA@counts[cell])/exp0length
    exp_orig<-cbind(orig[i],exp)
    result_exp<-rbind(result_exp,exp_orig)
  }
  colnames(result_exp)<-c("Sample","Exp")
  result_exp<-as.data.frame(result_exp)
  result_exp[["Status"]]<-"NA"
  result_exp[["cellstate"]] <- str
  h<-grep("*Healthy.",result_exp$Sample)
  result_exp[h,]$Status<-"Healthy"
  S<-grep("*COVID.19.Severe.",result_exp$Sample)
  result_exp[S,]$Status<-"COVID.19.Severe"
  M<-grep("*COVID.19.Mild.",result_exp$Sample)
  result_exp[M,]$Status<-"COVID.19.Mid"
  
  result_exp$Status <- factor(result_exp$Status,levels = c("Healthy", "COVID.19.Mid", "COVID.19.Severe"))
  result_exp<-as.data.frame(result_exp)
  result_exp$Exp<-as.numeric(result_exp$Exp)
  return(result_exp)
}

CCL3_myeloid_exptable <- celltable(myeloid_CCAresult,str='myeloid cells',gene = "CCL3")

my_comparisons=list(c("Healthy","COVID.19.Mid"),c("COVID.19.Mid","COVID.19.Severe"),c("Healthy","COVID.19.Severe"))

CCL3_myeloid_plot<-ggplot(CCL3_myeloid_exptable,aes(x=Status,y=Exp,fill=Status))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,fill=c("#98C3E5","#88B382","#FDCB86"),outlier.color = "white")+
  geom_jitter(aes(fill=Status),width=0.2,shape=21,size=3.5)+
  scale_fill_manual(values = c("#005393","#18641E","#FF9300"))+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in myeloid cells")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))

##将T细胞细胞挑出##
Tcells_CCAresult<-subset(COVID_CCAresult2,celltype=="T cells")
##CCL3在T细胞细胞中的表达情况##

CCL3_Tcells_exptable <- celltable(Tcells_CCAresult,str='T cells',gene = "CCL3")


CCL3_Tcells_plot<-ggplot(CCL3_Tcells_exptable,aes(x=Status,y=Exp,fill=Status))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,fill=c("#98C3E5","#88B382","#FDCB86"),outlier.color = "white")+
  geom_jitter(aes(fill=Status),width=0.2,shape=21,size=3.5)+
  scale_fill_manual(values = c("#005393","#18641E","#FF9300"))+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in T cells")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))

fig2G<-ggarrange(CCL3_allcells_plot,CCL3_myeloid_plot,CCL3_Tcells_plot,
                        common.legend = T,legend = "right",
                        ncol = 3, nrow = 1)


####F3A####
##亚型划分##
myeloid_CCAresult<-readRDS(paste0(result_path,"myeloid_defined.rds"))
Tcells_CCAresult<-readRDS(paste0(result_path,"Tcells_defined.rds"))

myeloid_CCAresult<-subset(myeloid_CCAresult,percent.mt<20)

myeloid_CCAresult$subpopulation<-factor(myeloid_CCAresult$subpopulation,
                                        levels = c("classical_monocyte",
                                                   "non-classical_monocyte",
                                                   "Macrophages",
                                                   "intermediate_monocyte"))

fig3A<-DimPlot(myeloid_CCAresult,reduction = "tsne",group.by = "subpopulation",
        cols = c("#F4C781","#42B1D3","#B54E93","#5FB46F"))

####F3B####
myeloid_CCAresult = SetIdent( myeloid_CCAresult, value = "subpopulation")
cluster_marker_data_myeloidsubtype <- FindAllMarkers(myeloid_CCAresult,only.pos = T, slot = "data",test.use = "MAST",
                                      assay = "SCT",logfc.threshold = 0.25,min.pct = 0.25)

write.table(cluster_marker_data_myeloidsubtype,paste0(result_path,"result3/myeloidsubtypes_markergene.txt"),quote = F,sep = "\t",row.names = F)

marker<-c("CD14",
          "FCGR3A","FCGR3B",
          "MS4A6A",
          "C1QC","HLA-DQA1",
          "PPBP","RGS18")

vlnmatrix<-GetAssayData(myeloid_CCAresult,assay = "SCT", slot = "data")

markerloc<-which(rownames(vlnmatrix)%in%marker)
metadata<-cbind(rownames(myeloid_CCAresult@meta.data),myeloid_CCAresult@meta.data)
colnames(metadata)[1]<-"barcode"
metadatavln<-cbind(metadata$barcode,as.character(metadata$subpopulation))
colnames(metadatavln)<-c("Barcode","Subpopulation")
metadatavln<-as.data.frame(metadatavln)
metadatavln$Subpopulation<-factor(metadatavln$Subpopulation,
                                  levels = c("classical_monocyte",
                                             "non-classical_monocyte",
                                             "Macrophages",
                                             "intermediate_monocyte"))
allgene_vln<-c()

for(i in 1:length(markerloc)){
  markerlocnew<-which(rownames(vlnmatrix)%in%marker[i])
  exp<-data.frame(marker[i],vlnmatrix[markerlocnew,])
  colnames(exp)<-c("Gene","log(TPM+1)")
  vln<-cbind(exp,metadatavln)
  allgene_vln<-rbind(allgene_vln,vln)
}

allgene_vln$Gene<-factor(allgene_vln$Gene,
                         levels = marker)

fig3B<-ggplot(allgene_vln,aes(x=Subpopulation,
                       y=`log(TPM+1)`,
                       fill=Subpopulation))+geom_violin()+
  scale_fill_manual(values = c("#F4C781","#42B1D3","#B54E93","#5FB46F"))+
  facet_grid(Gene~ .)+theme_bw()+
  theme(axis.text = element_text(size=8),  # 坐标轴上的文字
        axis.title = element_text(size=8, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position ="none")
####F3C####
cell_numtable<-as.data.frame(table(myeloid_CCAresult$celltype, 
                                   myeloid_CCAresult$orig.ident))

colnames(cell_numtable)<-c("Celltype","Sample","Number")

cell_numtable[["Status"]]<-c("NA")
cell_numtable$Status[1:24]<-("COVID.19.Mid")
cell_numtable$Status[25:88]<-("COVID.19.Severe")
cell_numtable$Status[89:112]<-("Healthy")


celltypeall<-c()
for(i in 1:length(unique(cell_numtable$Sample))){
  celltype1<-subset(cell_numtable,Sample==unique(cell_numtable$Sample)[i],select = c("Celltype","Sample","Number","Status"))
  prop<-prop.table(celltype1$Number)
  celltype1<-cbind(celltype1,prop)
  celltypeall<-rbind(celltypeall,celltype1)
}

celltypeall$Celltype<-factor(celltypeall$Celltype,
                             levels = c("classical_monocyte",
                                        "non-classical_monocyte",
                                        "intermediate_monocyte",
                                        "Macrophages"))

celltypeall$Status<-factor(celltypeall$Status,
                           levels = c("Healthy","COVID.19.Mid",
                                      "COVID.19.Severe"))

celltypeall$Celltype<-factor(celltypeall$Celltype,
                             levels = c("classical_monocyte",
                                        "non-classical_monocyte",
                                        "Macrophages",
                                        "intermediate_monocyte"))




fig3C<-ggplot(celltypeall,aes(x=Status,y=prop,fill=Celltype))+
  geom_bar(stat = "summary",
           fun=mean,
           position = position_dodge(),width = 0.9)+
  scale_fill_manual(values = c("#F4C781","#42B1D3","#B54E93","#5FB46F"))+
  stat_summary(fun.data = "mean_se",geom = "errorbar", 
               colour="black",
               position = position_dodge(0.9),width=0.4,size=0.8)+
  stat_compare_means()+
  theme_bw()+
  theme(axis.text = element_text(size=8),  # 坐标轴上的文字
        axis.title = element_text(size=8, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position ="none")




####F3D####
fig3C<-FeaturePlot(myeloid_CCAresult, reduction = "tsne",slot ="scale.data",
                   features = c("CCL3"),label = F,label.size = 3,pt.size = 1)+
  scale_color_gradientn(colours = 
                          c("#F8F6BF","#FFDCA0","#FAB880","#F69B6F","#F47F62",
                            "#F06C63","#E55065","#CD4271","#BF3B77","#AC327B",
                            "#9E2D7E","#7C2980","#4F2977","#2C1D57","#191135"))


####F3E####
fig3E<-DotPlot(myeloid_CCAresult, features = c("CCL3"),group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+  
  scale_color_gradientn(colours = 
                          c("#584B97","#4266A5","#3584B0","#4FA2A7","#6BBD9E",
                            "#8CCA9E","#AFD69E","#CCE199","#E4EBA0","#F0F2B4",
                            "#F7F1AD","#F8E093","#F8E093","#F6CA7B","#F4B167",
                            "#EE9154","#E96F46","#DA5646","#CA3F4B","#B02045",
                            "#951D41"))+theme(legend.key.size = unit(5, "pt"),
                                              legend.title =element_text( size = 8),
                                              axis.text.x = element_text(size = 8,angle = 90))

####F3F####
##将经典单核挑出##
cm_CCAresult<-subset(myeloid_CCAresult,celltype=="classical_monocyte")
##CCL3在经典单核细胞中的表达情况##
CCL3_cm_exptable <- celltable(cm_CCAresult,str='classical monocyte',gene = "CCL3")

my_comparisons=list(c("Healthy","COVID.19.Mid"),c("COVID.19.Mid","COVID.19.Severe"),c("Healthy","COVID.19.Severe"))

CCL3_cm_plot<-ggplot(CCL3_cm_exptable,aes(x=Status,y=Exp,fill=Status))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,fill=c("#98C3E5","#88B382","#FDCB86"),outlier.color = "white")+
  geom_jitter(aes(fill=Status),width=0.2,shape=21,size=3.5)+
  scale_fill_manual(values = c("#005393","#18641E","#FF9300"))+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in classical monocyte")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))

##将非经典单核挑出##
nm_CCAresult<-subset(myeloid_CCAresult,celltype=="non-classical_monocyte")
##CCL3在非经典单核细胞中的表达情况##
CCL3_nm_exptable <- celltable(nm_CCAresult,str='non-classical_monocyte',gene = "CCL3")

CCL3_nm_plot<-ggplot(CCL3_nm_exptable,aes(x=Status,y=Exp,fill=Status))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,fill=c("#98C3E5","#88B382","#FDCB86"),outlier.color = "white")+
  geom_jitter(aes(fill=Status),width=0.2,shape=21,size=3.5)+
  scale_fill_manual(values = c("#005393","#18641E","#FF9300"))+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in non-classical monocyte")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))

##将巨噬单核挑出##
ma_CCAresult<-subset(myeloid_CCAresult,celltype=="Macrophages")
##CCL3在巨噬细胞中的表达情况##
CCL3_ma_exptable <- celltable(ma_CCAresult,str='Macrophages',gene = "CCL3")

my_comparisons=list(c("Healthy","COVID.19.Mid"),c("COVID.19.Mid","COVID.19.Severe"),c("Healthy","COVID.19.Severe"))

CCL3_ma_plot<-ggplot(CCL3_ma_exptable,aes(x=Status,y=Exp,fill=Status))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,fill=c("#98C3E5","#88B382","#FDCB86"),outlier.color = "white")+
  geom_jitter(aes(fill=Status),width=0.2,shape=21,size=3.5)+
  scale_fill_manual(values = c("#005393","#18641E","#FF9300"))+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in Macrophages")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))

##合并CCL3在髓系亚型中的表达##
fig3F<-ggarrange(CCL3_cm_plot,CCL3_nm_plot,CCL3_ma_plot,
                 common.legend = T,legend = "right",
                 ncol = 3, nrow = 1)

####F3G####
classical_CCAresult<-subset(myeloid_CCAresult,subset= subpopulation == "classical_monocyte")

classical_CCAresult = SetIdent(classical_CCAresult, value = "Status")
cluster_marker_data_classical_tumor <- FindMarkers(classical_CCAresult,only.pos = T,ident.1 = "COVID.19.Severe",ident.2 = "Healthy",slot = "data",test.use = "MAST",
                                                   assay = "SCT",logfc.threshold = 0.25,min.pct = 0.25)

cluster_marker_data_classical_tumor<-cbind(cluster_marker_data_classical_tumor,rownames(cluster_marker_data_classical_tumor))
colnames(cluster_marker_data_classical_tumor)[6]<-"Gene"

cluster_marker_data_classical_control <-FindMarkers(classical_CCAresult,only.pos = T,ident.1 = "Healthy",ident.2 = "COVID.19.Severe",slot = "data",test.use = "MAST",
                                                    assay = "SCT",logfc.threshold = 0.25,min.pct = 0.25)

cluster_marker_data_classical_control<-cbind(cluster_marker_data_classical_control,rownames(cluster_marker_data_classical_control))
colnames(cluster_marker_data_classical_control)[6]<-"Gene"

write.table(cluster_marker_data_classical_tumor,paste0(result_path,"result3/cluster_marker_data_classical_tumor.txt"),quote = F,sep = "\t",row.names = F)
write.table(cluster_marker_data_classical_control,paste0(result_path,"result3/cluster_marker_data_classical_control.txt"),quote = F,sep = "\t",row.names = F)

##巨噬##
macrophages_CCAresult<-subset(myeloid_CCAresult,subset= subpopulation == "Macrophages")

macrophages_CCAresult = SetIdent(macrophages_CCAresult, value = "Status")
cluster_marker_data_macrophages_tumor <- FindMarkers(macrophages_CCAresult,only.pos = T,ident.1 = "COVID.19.Severe",ident.2 = "Healthy",slot = "data",test.use = "MAST",
                                                     assay = "SCT",logfc.threshold = 0.25,min.pct = 0.25)
cluster_marker_data_macrophages_tumor<-cbind(cluster_marker_data_macrophages_tumor,rownames(cluster_marker_data_macrophages_tumor))
colnames(cluster_marker_data_macrophages_tumor)[6]<-"Gene"

genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

classical_CCAresult_cs<-subset(classical_CCAresult,subset= Status == "COVID.19.Mid",invert =T)

classical_CCAresult_cs<-NormalizeData(classical_CCAresult_cs)

Idents(classical_CCAresult_cs) <- "orig.ident" 
expr <- AverageExpression(classical_CCAresult_cs, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

gene_con<-cluster_marker_data_classical_control$Gene
gene_tum<-cluster_marker_data_classical_tumor$Gene
geneall<-c(gene_con,gene_tum)

geneloc<-which(rownames(expr)%in%geneall)

expr2<-expr[geneloc,]

gsva.res <- gsva(expr2, genesets, method="gsva") 
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)

pheatmap::pheatmap(gsva.res, show_colnames = T, 
                   scale = "row",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

group_list <- data.frame(sample = colnames(gsva.df)[-1], 
                         group = c(rep("case", 11), rep("con", 6),rep("case",3)))

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva.res)

contrast.matrix <- makeContrasts(con-case, levels = design)

fit <- lmFit(gsva.res, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)

cutoff <- 0
df$group <- cut(df$score, breaks = c(-Inf, cutoff, Inf),labels = c(1,2))

sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)


fig3G<-ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c("#F0CEB3", "#BE7A89"), guide = FALSE) + 
  #画2条虚线
  geom_hline(yintercept = c(-1,1), 
             color="white",
             linetype = 2, #画虚线
             size = 0.5) + #线的粗细
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0.1, label=ID, color = group),
            size = 2, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 2, hjust = "outward") +  
  scale_colour_manual(values = c("black","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴




####F4A####
mar<-subset(celltypeall,subset=Celltype=="Macrophages")

my_comparisons=list(c("Healthy","COVID.19.Mid"),c("COVID.19.Mid","COVID.19.Severe"),c("Healthy","COVID.19.Severe"))

fig4B<-ggplot(mar,aes(x=Status,y=prop,fill=Celltype))+
  geom_bar(stat = "summary",
           fun=mean,
           position = position_dodge(),width = 0.5)+
  scale_fill_manual(values = c("#B54E93"))+
  stat_summary(fun.data = "mean_se",geom = "errorbar", 
               colour="black",
               position = position_dodge(0.9),width=0.2,size=1)+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  theme_bw()+
  theme(axis.text = element_text(size=8),  # 坐标轴上的文字
        axis.title = element_text(size=8, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position ="none")

####F4B####
marker_myeloidtop<-cluster_marker_data_myeloidsubtype %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

gene_ma<-marker_myeloidtop$gene[1:10]

gene_ma_go <- bitr(gene_ma,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
gene_ma_goresult <- enrichGO(gene_ma_go$ENTREZID,OrgDb = "org.Hs.eg.db",
                             keyType = "ENTREZID",ont = "ALL",readable = T)
ma_go.res<-data.frame(gene_ma_goresult)

ma_goBP <- subset(ma_go.res,subset = (ONTOLOGY == "BP"))[1:5,]
ma_goCC <- subset(ma_go.res,subset = (ONTOLOGY == "CC"))[1:5,]
ma_goMF <- subset(ma_go.res,subset = (ONTOLOGY == "MF"))[1:5,]
ma_go.df <- rbind(ma_goBP,ma_goCC,ma_goMF)

ma_go.df$Description <- factor(ma_go.df$Description,levels = rev(ma_go.df$Description))

fig4C<-ggplot(ma_go.df,aes(x=Description,y=Count))+
  geom_point(aes(size=5,color=-log(p.adjust),shape=ONTOLOGY))+
  coord_flip()+theme_bw()+
  scale_color_gradientn(colours = c("#8390A0","#C6B9AF","#E6CBB7",
                                    "#F2D6BB","#F0CEB3","#E8BDAC",
                                    "#C89191","#BE7A89","#9B516A",
                                    "#8F3261","#892859"))

####F4C####
expr_matrix <- as(as.matrix(GetAssayData(myeloid_CCAresult,slot="counts")), 'sparseMatrix')
f_data <- data.frame(gene_short_name = rownames(myeloid_CCAresult),gene_id = row.names(myeloid_CCAresult))
rownames(f_data)<-rownames(myeloid_CCAresult)
fd <- new('AnnotatedDataFrame', data = f_data)


p_data <- myeloid_CCAresult@meta.data 
rownames(p_data)<-colnames(myeloid_CCAresult)
pd <- new('AnnotatedDataFrame', data = p_data) 
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds_allmeyloid <- newCellDataSet(expr_matrix,
                                 phenoData = pd,
                                 featureData = fd,
                                 expressionFamily = negbinomial.size())


cds_allmeyloid <- estimateSizeFactors(cds_allmeyloid)
cds_allmeyloid <- estimateDispersions(cds_allmeyloid)
#筛选基因,这里可以根据自己的需要筛选特定的基因
express_genes <- VariableFeatures(myeloid_CCAresult)
cds_allmeyloid <- setOrderingFilter(cds_allmeyloid, express_genes)

cds_allmeyloid <- reduceDimension(cds_allmeyloid, max_components = 2,
                                  method = 'DDRTree')

cds_allmeyloid <- orderCells(cds_allmeyloid,reverse = F)
saveRDS(cds_allmeyloid,paste0(result_path,"allmeyloid_monocle.rds"))

cds_allmeyloid<-readRDS(paste0(result_path,"allmeyloid_monocle.rds"))

plot_cell_trajectory(cds_allmeyloid,color_by="Status", 
                     size=1,show_backbone=TRUE)

plot_cell_trajectory(cds_allmeyloid,color_by="Pseudotime", 
                     size=1,show_backbone=TRUE,)

fig4A1<-plot_complex_cell_trajectory(cds_allmeyloid,x=1,y=2,color_by = "subpopulation",
                                     show_backbone=TRUE,root_states = 5)+
  scale_colour_manual(values = c("#F4C781","#42B1D3","#B54E93","#5FB46F"))

cds_allmeyloid$Status<-factor(cds_allmeyloid$Status,levels = c("Healthy","COVID.19.Mid","COVID.19.Severe"))

fig4A2<-plot_complex_cell_trajectory(cds_allmeyloid,x=1,y=2,color_by = "Status",
                                     show_backbone=TRUE,root_states = 5)+
  scale_colour_manual(values = c("#98C3E5","#88B382","#FDCB86"))

ggarrange(fig4A1,fig4A2, ncol = 2, nrow = 1)

####F4D####
monocyte_table <- celltable(cm_CCAresult,str='classical_monocyte',gene = "CCL3")
Macrophages_table <- celltable(ma_CCAresult,str='Macrophages',gene = "CCL3")

allcelltable<-rbind(monocyte_table,Macrophages_table)

my_comparisons<-list(c("classical_monocyte","Macrophages"))

fig4D<-ggplot(allcelltable,aes(x=cellstate,y=Exp,fill=cellstate))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,outlier.color = "white")+
  scale_fill_manual(values = c("#F4C781","#B54E93"))+
  geom_jitter(aes(fill=cellstate),width=0.1,shape=21,size=3.5)+
  facet_grid(. ~ Status)+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in each group")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))


####F4F####
fig4f<-FeaturePlot(COVID_CCAresult2, reduction = "tsne",slot ="data",
            features = c("CCR1"),label = F,label.size = 3,pt.size = 0.5)+
  scale_color_gradientn(colours = 
                          c("#F8F6BF","#FFDCA0","#FAB880","#F69B6F","#F47F62",
                            "#F06C63","#E55065","#CD4271","#BF3B77","#AC327B",
                            "#9E2D7E","#7C2980","#4F2977","#2C1D57","#191135"))


####F4G####
FeaturePlot(myeloid_CCAresult, reduction = "tsne",slot ="data",
            features = c("CCR1"),label = F,label.size = 3,pt.size = 1)+
  scale_color_gradientn(colours = 
                          c("#F8F6BF","#FFDCA0","#FAB880","#F69B6F","#F47F62",
                            "#F06C63","#E55065","#CD4271","#BF3B77","#AC327B",
                            "#9E2D7E","#7C2980","#4F2977","#2C1D57","#191135"))

####补充####
cell_numtable<-as.data.frame(table(COVID_CCAresult2$celltype, 
                                   COVID_CCAresult2$orig.ident))

colnames(cell_numtable)<-c("Celltype","Sample","Number")

cell_numtable[["Status"]]<-c("NA")
cell_numtable$Status[1:36]<-("mild")
cell_numtable$Status[37:132]<-("severe")
cell_numtable$Status[133:168]<-("healthy")


celltypeall<-c()
for(i in 1:length(unique(cell_numtable$Sample))){
  celltype1<-subset(cell_numtable,Sample==unique(cell_numtable$Sample)[i],select = c("Celltype","Sample","Number","Status"))
  prop<-prop.table(celltype1$Number)
  celltype1<-cbind(celltype1,prop)
  celltypeall<-rbind(celltypeall,celltype1)
}


col_celltypefine<-c("#98C3E5","#88B382","#FDCB86")
celltype<-unique(COVID_CCAresult2$celltype)
my_comparisons=list(c("healthy","mild"),c("mild","severe"),c("healthy","severe"))



for (i in 1:length(celltype)) {
  celltype2<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
  p1<-ggplot(celltype2,aes(Status,prop,fill=Status))+
    geom_boxplot()+
    ggtitle("")+
    scale_fill_manual(values = col_celltypefine)+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+
    ggtitle(celltype[i])+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(result_path,"cell/cell_prop/",celltype[i],".pdf"),plot = p1,width = 5, height = 4)
  
  healthy_median_data<-subset(celltype2,Status == "healthy")
  healthy_median<-median(healthy_median_data$prop)
  mild_median_data<-subset(celltype2,Status == "mild")
  mild_median<-median(mild_median_data$prop)
  severe_median_data<-subset(celltype2,Status == "severe")
  severe_median<-median(severe_median_data$prop)
  result_median<-cbind(healthy_median,mild_median,severe_median)
  colnames(result_median)<-c("healthy","mild","severe")
  write.table(result_median,paste0(result_path,"cell/cell_prop/",celltype[i],".txt"),quote=F,sep="\t",row.names=F)
  
}


cell_numtable<-as.data.frame(table(myeloid_CCAresult$celltype, 
                                   myeloid_CCAresult$orig.ident))

colnames(cell_numtable)<-c("Celltype","Sample","Number")

cell_numtable[["Status"]]<-c("NA")
cell_numtable$Status[1:24]<-("mild")
cell_numtable$Status[25:88]<-("severe")
cell_numtable$Status[89:112]<-("healthy")


celltypeall<-c()
for(i in 1:length(unique(cell_numtable$Sample))){
  celltype1<-subset(cell_numtable,Sample==unique(cell_numtable$Sample)[i],select = c("Celltype","Sample","Number","Status"))
  prop<-prop.table(celltype1$Number)
  celltype1<-cbind(celltype1,prop)
  celltypeall<-rbind(celltypeall,celltype1)
}


col_celltypefine<-c("#98C3E5","#88B382","#FDCB86")
celltype<-unique(myeloid_CCAresult$celltype)
my_comparisons=list(c("healthy","mild"),c("mild","severe"),c("healthy","severe"))


for (i in 1:length(celltype)) {
  celltype2<-subset(celltypeall,Celltype==celltype[i],select = c("Celltype","Sample","Number","Status","prop"))
  p1<-ggplot(celltype2,aes(Status,prop,fill=Status))+
    geom_boxplot()+
    ggtitle("")+
    scale_fill_manual(values = col_celltypefine)+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),panel.grid = element_blank())+
    guides(fill=guide_legend(title=NULL))+
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+
    ggtitle(celltype[i])+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(result_path,"cell/cell_prop2/",celltype[i],".pdf"),plot = p1,width = 5, height = 4)
  
  healthy_median_data<-subset(celltype2,Status == "healthy")
  healthy_median<-median(healthy_median_data$prop)
  mild_median_data<-subset(celltype2,Status == "mild")
  mild_median<-median(mild_median_data$prop)
  severe_median_data<-subset(celltype2,Status == "severe")
  severe_median<-median(severe_median_data$prop)
  result_median<-cbind(healthy_median,mild_median,severe_median)
  colnames(result_median)<-c("healthy","mild","severe")
  write.table(result_median,paste0(result_path,"cell/cell_prop2/",celltype[i],".txt"),quote=F,sep="\t",row.names=F)
  
}

monocyte_table <- celltable(cm_CCAresult,str='classical_monocyte',gene = "CCL3")
Macrophages_table <- celltable(ma_CCAresult,str='Macrophages',gene = "CCL3")
ncms_CCAresult<-subset(myeloid_CCAresult,celltype=="non-classical_monocyte")
ncms_table<- celltable(ma_CCAresult,str='non-classical_monocyte',gene = "CCL3")
im_CCAresult<-subset(myeloid_CCAresult,celltype=="intermediate_monocyte")
im_table<- celltable(im_CCAresult,str='intermediate_monocyte',gene = "CCL3")


allcelltable<-rbind(monocyte_table,ncms_table,Macrophages_table)

my_comparisons<-list(
  c("classical_monocyte","non-classical_monocyte"),
  c("non-classical_monocyte","Macrophages"),
  c("classical_monocyte","Macrophages"))

allcelltable$cellstate<-factor(allcelltable$cellstate,
                               levels = c("classical_monocyte",
                                          "non-classical_monocyte",
                                          "Macrophages"))

ggplot(allcelltable,aes(x=cellstate,y=Exp,fill=cellstate))+
  stat_boxplot(geom = "errorbar",width= 0.2,size=1)+
  geom_boxplot(size=0.8,outlier.color = "white")+
  scale_fill_manual(values = c("#F4C781","#42B1D3","#B54E93"))+
  geom_jitter(aes(fill=cellstate),width=0.1,shape=21,size=3.5)+
  stat_compare_means(comparisons = my_comparisons,aes(label= ..p.signif..))+
  ggtitle("CCL3 in each group")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=15,hjust = 0.5))+  
  guides(fill=guide_legend(title=NULL))
