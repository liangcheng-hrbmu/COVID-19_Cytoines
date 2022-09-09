library(TwoSampleMR)
library(RadialMR)

result_path<-c("W:/Mypaper/COVID19_SingleCell/COVID19-brief report/COVID_19---IL/")

########MR analysis#########
MRdata<-list.files(paste0(result_path,"2.MRdata/"))

for(i in 1:length(MRdata)){
  MRfilepath<-paste0(result_path,c("2.MRdata/"),MRdata[i])
  MRfile<-read.table(MRfilepath,header=T,sep="\t")
  p0.05loc<-which(MRfile[,8]>0.05)
  MRanalysis_data<-MRfile[p0.05loc,]
  MRanalysis_data2<-cbind(MRanalysis_data[,1:3],MRanalysis_data[,12],MRanalysis_data[,4:11])
  colnames(MRanalysis_data2)[4]<-c("base_pair_location_Grch38")
  MRoutpath<-paste0(result_path,c("3.MRdata_pvalue0.05/MRdata_P_value0.05_"),MRdata[i])
  write.table(MRanalysis_data2,MRoutpath,quote=F,sep="\t",row.names = F)
  print(i)
}

##part2###
analysisdata<-list.files(paste0(result_path,"3.MRdata_pvalue0.05"))

re_result<-c()
heter_result<-c()
pleiotropy_result<-c()



result_MR_12cytokines<-c(3,12,14,15,
                         20,23,26,27,
                         28,31,36,40)
length(result_MR_12cytokines)
as.numeric(result_MR_12cytokines[i])

for(i in 1:length(analysisdata)){
  
  filename1<-paste0(result_path,c("3.MRdata_pvalue0.05/"),
              analysisdata[i])
  
  a<-read_exposure_data(filename =filename1,sep = "\t",snp_col = "rsID", 
                        beta_col ="beta_COVID19", 
                        se_col = "SE_COVID19", 
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        pval_col ="P.value_COVID19")
  
  outcomename<-strsplit(analysisdata[i],".txt")[[1]]
  
  b<-read_outcome_data(filename =filename1,sep = "\t",snp_col = "rsID",
                       beta_col =paste0("beta_",strsplit(outcomename,"MRdata_P_value0.05_")[[1]][2]),
                       se_col = paste0("SE_",strsplit(outcomename,"MRdata_P_value0.05_")[[1]][2]), 
                       effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele",
                       pval_col =paste0("P.value_",strsplit(outcomename,"MRdata_P_value0.05_")[[1]][2]) )
  
  
  
  a2<-clump_data(
    a,
    clump_kb = 1,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1
  )
  
  
  b2<-b[which(b$SNP%in%a2$SNP),]
  
  
  data<-harmonise_data(exposure_dat = a2,outcome_dat = b2)
  data$exposure<-"COVID19"
  data$outcome<-strsplit(outcomename,"MRdata_P_value0.05_")[[1]][2]
  res<-mr(data,method_list = c("mr_ivw","mr_weighted_median"))
  re_result<-rbind(res,re_result)
  
  heter<-mr_heterogeneity(data)
  heter_result<-rbind(heter,heter_result)
  
  pleiotropy<-mr_pleiotropy_test(data)
  pleiotropy_result<-rbind(pleiotropy,pleiotropy_result)
  
  res_single <- mr_singlesnp(data)
  forestplot<-mr_forest_plot(res_single)
  
  leaveoneout_result<-mr_leaveoneout(data)
  leaveoneout_plot<-mr_leaveoneout_plot(leaveoneout_result)
  
  mrplot<-mr_scatter_plot(res,data)
  
  
  
  write.table(res_single,paste0(result_path,c("4.MRresult/singleSNP/singlesnp_file/"),analysisdata[i]),quote=F,sep="\t",row.names = F)
  write.table(leaveoneout_result,paste0(result_path,c("4.MRresult/leave_one_out/leaveoneout_file/"),analysisdata[i]),quote=F,sep="\t",row.names = F)
  
  plotresult_name<-strsplit(MRdata[i],".txt")[[1]]
  pdf(paste0(result_path,c("4.MRresult/singleSNP/forestplot/"),plotresult_name,c("_forest.pdf")))
  forestplot
  dev.off()
  
  pdf(paste0(result_path,c("4.MRresult/leave_one_out/leaveoneout_plot/"),plotresult_name,c("_leaveoneout.pdf")))
  leaveoneout_plot
  dev.off()  

  pdf(paste0(result_path,c("4.MRresult/mrplot/"),plotresult_name,c("_scatterplot.pdf")))
  mrplot
  dev.off()
  }

write.table(re_result,paste0(result_path,"4.MRresult/MRresult2_allMethod(Exposure_COVID19).txt"),quote=F,sep="\t",row.names = F)
write.table(heter_result,paste0(result_path,"4.MRresult/MRresult_heter(Exposure_COVID19).txt"),quote=F,sep="\t",row.names = F)
write.table(pleiotropy_result,paste0(result_path,"4.MRresult/MRresult_pleiotropy(Exposure_COVID19).txt"),quote=F,sep="\t",row.names = F)




########part3#########
filename1<-paste0(result_path,"3.MRdata_pvalue0.05/MRdata_P_value0.05_IL1b.txt")
a<-read_exposure_data(filename =filename1,sep = "\t",snp_col = "rsID", 
                      beta_col ="beta_COVID19", 
                      se_col = "SE_COVID19", 
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      pval_col ="P.value_COVID19")

outcomename<-strsplit("MRdata_P_value0.05_IL1b.txt",".txt")[[1]]

b<-read_outcome_data(filename =filename1,sep = "\t",snp_col = "rsID",
                     beta_col =paste0("beta_",strsplit(outcomename,"MRdata_P_value0.05_")[[1]][2]),
                     se_col = paste0("SE_",strsplit(outcomename,"MRdata_P_value0.05_")[[1]][2]), 
                     effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele",
                     pval_col =paste0("P.value_",strsplit(outcomename,"MRdata_P_value0.05_")[[1]][2]) )




a2<-clump_data(
  a,
  clump_kb = 1,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1
)


b2<-b[which(b$SNP%in%a2$SNP),]


data<-harmonise_data(exposure_dat = a2,outcome_dat = b2)

res_single <- mr_singlesnp(data)
mr_forest_plot(res_single)


##4##
data_MRresult<-read.table("D:/1.Wang/Second_postgraduate/Term_two/2021.4/COVID19-brief report/COVID_19---IL/4.MRresult/IVW_WeightMedian/IVW_WeightMedian_plot1.txt",header=T,sep="\t")

data_MRresult<-as.data.frame(data_MRresult)

data_MRresult$method<-factor(data_MRresult$method,
                             levels =c( "MR_Egger_intercept",
                                        "Weighted median","Inverse variance weighted"))

data_MRresult$sig<-factor(data_MRresult$sig,
                          levels =c( "Promote",
                                     "Inhibition","p>0.05","p<0.05"))

m1<-ggplot(data_MRresult,aes(x = outcome, y = method)) +
  
  geom_point(aes(size=b,colour=sig))+scale_colour_nejm()  +
  theme_bw()+theme(panel.grid=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y = element_blank())+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45))




fig1c<-read.table("W:/Mypaper/COVID19_SingleCell/COVID19-brief report/COVID_19---IL/4.MRresult/fig1_Cggplot.txt",header = T,sep="\t")

cyto<-c("Eotaxin","IL16","IL18","IL1b",
        "IL5","IL8","MCP1","MCP3",
        "MCSF","MIP1a","SCGFb","TRAIL")

loc<-which(fig1c$outcome%in%cyto)

fig1C_12<-fig1c[loc,]

xmin = fig1C_12$b + 1.96* fig1C_12$se
xmax = fig1C_12$b - 1.96* fig1C_12$se

fig1C_12<-cbind(fig1C_12,xmin,xmax)


fig1C_12$outcome<-factor(fig1C_12$outcome,levels = c("MCP3","MCSF","SCGFb","IL1b",
                                                     "MIP1a","MCP1","IL16","IL5",
                                                     "IL8","TRAIL","Eotaxin","IL18"))

ggplot(fig1C_12,aes(x=b,y=outcome,colour=outcome,
                    xmin=xmin,xmax=xmax))+
  geom_vline(xintercept = 0, 
             linetype = "dotted",size=0.5)+
  geom_errorbarh(height=0.3) +
  geom_point(aes(x=b),size=2)+
  facet_wrap(vars(method))+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  theme(panel.grid=element_blank())+
  guides(colour = guide_legend(reverse=TRUE))+
  scale_colour_manual(values = c("#4D76A0","#9CC7E0","#E98A2A","#F4B87A",
                                 "#8AC27A","#B09433","#489490","#D65557",
                                 "#EC9796","#C96D8E","#AA769B","#CFA3C3"))



