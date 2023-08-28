
EUR_sumfilegout_males<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/EUR_maleGRS_summary_table.txt")
EUR_resulttable_males<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/EUR_maleGRS_results_table.txt")
EUR_sumfilegout_females<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/EUR_femaleGRS_summary_table.txt")
EUR_resulttable_females<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/EUR_femaleGRS_results_table.txt")
EUR_roc_males<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/EUR_maleGRS_ROC.txt")
EUR_roc_females<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/EUR_femaleGRS_ROC.txt")



TAMA_sumfilegout_females<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/TAMA_femaleGRS_summary_table.txt")
TAMA_resulttable_females<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/TAMA_femaleGRS_results_table.txt")
TAMA_sumfilegout_males<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/TAMA_maleGRS_summary_table.txt")
TAMA_resulttable_males<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/TAMA_maleGRS_results_table.txt")
TAMA_roc_males<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/TAMA_maleGRS_ROC.txt")
TAMA_roc_females<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/TAMA_femaleGRS_ROC.txt")



EUR_males<-merge(EUR_sumfilegout_males, EUR_resulttable_males, by.x="Var1",by.y="name", all=T)
EUR_males$ancestry<-"European"
EUR_females<-merge(EUR_sumfilegout_females, EUR_resulttable_females, by.x="Var1",by.y="name", all=T)
EUR_females$ancestry<-"European"

TAMA_males<-merge(TAMA_sumfilegout_males, TAMA_resulttable_males, by.x="Var1",by.y="name", all=T)
TAMA_males$ancestry<-"Trans-ancestry "
TAMA_females<-merge(TAMA_sumfilegout_females, TAMA_resulttable_females, by.x="Var1",by.y="name", all=T)
TAMA_females$ancestry<-"Trans-ancestry"

males<-rbind(EUR_males,TAMA_males)
females<-rbind(EUR_females,TAMA_females)
rm(EUR_females,EUR_males, EUR_resulttable_females,EUR_resulttable_males, EUR_sumfilegout_females,EUR_sumfilegout_males, TAMA_females, TAMA_males, TAMA_resulttable_females, TAMA_resulttable_males, TAMA_sumfilegout_females, TAMA_sumfilegout_males)

library(tidyverse)
library(magrittr)
numcols = c("OR","LCI","UCI", "P","N")
library(ggplot2)
library(ggrepel)
library(ggh4x)
library(gridExtra)
library(scales)


males$plotb<-c(1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10)
females$plotb<-c(1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10)

EUR_roc_males$ancestry<-"European"
TAMA_roc_males$ancestry<-"Trans-ancestry"

ggroc_males<-rbind(EUR_roc_males,TAMA_roc_males)
ggroc_males$adjusted_specificity<-1-ggroc_males$specificity
ggroc_males$name<-factor(ggroc_males$name, levels = c("Male PRS", "Demographic", "Combined"))




# set colour palette
#devtools::install_github("G-Thomson/Manu")
library(Manu)

gwas_palette <- c("#7ACCD7", "#115896", "#7C6C65", "#4C4C53", "#BA2F00", "#B865A1", "#4C4C53")
gwas_palette <- sort(gwas_palette)

A<-ggplot(males, aes(as.factor(plotb),Freq))+theme_light(base_size = 12)+
  geom_col(color="black",alpha=1, aes(fill=ancestry))+
  ylab("Count")+
  scale_fill_manual(values = c("#4C4C53","#4C4C53")) +
  scale_x_discrete(expand = c(0.1,0.1))+
  theme(axis.text.x = element_blank(), legend.position = "none", axis.title.x = element_blank())+
  facet_wrap(~ancestry)
A

B<-ggplot(males, aes(plotb,prevalencegout))+
  geom_line(aes(colour=ancestry))+theme_light(base_size = 12)+
  scale_colour_manual(values=c("#4C4C53","#4C4C53"))+
  ylab("Gout Prevalence (%)")+
  theme(axis.text.x = element_blank(), legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 18, b = 0, l = 0)))+    
  scale_x_continuous(labels = NULL, breaks = 1:10, minor_breaks=NULL, expand = c(0.1, 0.1))+ 
  facet_wrap(~ancestry)
B

C<-ggplot(data=males, aes(x=as.factor(Var1), y=OR)) + theme_light(base_size = 12)+
  geom_bar(stat="identity", color="black", alpha=1, na.rm = T, aes(fill=ancestry)) +
  scale_fill_manual(values=c("#4C4C53","#4C4C53"))+
    labs( y= "Odds Ratio", x = "Polygenic Risk Score (bins)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none", axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(t = 0, r = -1, b = 0, l = 0)))+ 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.2, position=position_dodge(.9), na.rm = T) +
  geom_text_repel(aes(x=Var1, y=labelpos, label = sig), na.rm = T, nudge_y = ifelse(males$labelpos <=1, -0.1, +0.1),size =3)+
  scale_y_log10(expand=c(0.1,0.1), labels = comma)+
  scale_x_discrete(expand = c(0.1,0.1))+
  facet_grid(~ancestry, scale="free", space="free_x")
C
D <-ggplot(data = ggroc_males, aes(x=adjusted_specificity, y=sensitivity))+
  geom_line(aes(colour=name))+theme_light(base_size = 12)+
  scale_colour_manual(values=c("#7ACCD7", "#B865A1","#7C6C65"))+
  labs(x = "1-Specificity", y = "Sensitivity", colour="AUROC model")+
  theme(legend.position = "bottom", legend.margin = margin(t = 0, r = 0, b = 0, l = 0), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.x = element_text(margin = margin(t = 10, b=0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+    
  facet_wrap(~ancestry)
D

#strip.background =element_rect(fill="Black"),

png(filename = "/Volumes//biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/male_GRS_all.png", width = 7, height = 12, units = "in", res = 300)
grid.arrange(A + 
               theme(plot.margin = unit(c(0,0.25,0,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "A"),
             B + 
               theme(plot.margin = unit(c(0,0.25,0,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "B"),
             C + 
               theme(plot.margin = unit(c(0,0.25,0,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "C"),
             D +
               theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "D"),
             layout_matrix = matrix( c(1,2,3,4), ncol = 1), heights = c(5,5,7,6))
dev.off()




EUR_roc_females$ancestry<-"European"
TAMA_roc_females$ancestry<-"Trans-ancestry"

ggroc_females<-rbind(EUR_roc_females,TAMA_roc_females)
ggroc_females$adjusted_specificity<-1-ggroc_females$specificity
ggroc_females$name<-factor(ggroc_females$name, levels = c("Female PRS", "Demographic", "Combined"))





A<-ggplot(females, aes(as.factor(plotb),Freq))+theme_light(base_size = 12)+
  geom_col(color="black",alpha=1, aes(fill=ancestry))+
  ylab("Count")+
  scale_fill_manual(values = c("#6AA086","#6AA086")) +
  scale_x_discrete(expand = c(0.1,0.1))+
  theme(axis.text.x = element_blank(), legend.position = "none", axis.title.x = element_blank())+ 
  facet_wrap(~ancestry)
A

B<-ggplot(females, aes(plotb,prevalencegout))+
  geom_line(aes(colour=ancestry))+theme_light(base_size = 12)+
  scale_colour_manual(values=c("#6AA086","#6AA086"))+
  ylab("Gout Prevalence (%)")+
  theme(axis.text.x = element_blank(), legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 26, b = 0, l = 0)))+  
  scale_x_continuous(labels = NULL, breaks = 1:10, minor_breaks=NULL, expand = c(0.1, 0.1))+ 
  facet_wrap(~ancestry)
B

C<-ggplot(data=females, aes(x=as.factor(Var1), y=OR)) + theme_light(base_size = 12)+
  geom_bar(stat="identity", color="black", alpha=1, na.rm = T, aes(fill=ancestry)) +
  scale_fill_manual(values=c("#6AA086","#6AA086"))+
  labs(y="Odds Ratio", x = "Polygenic Risk Score (bins)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none", axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)))+ 
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.2, position=position_dodge(.9), na.rm = T) +
  geom_text_repel(aes(x=Var1, y=labelpos, label = sig), na.rm = T, nudge_y = ifelse(females$labelpos <=1, -0.1, +0.1),size =3)+
  scale_y_log10(expand=c(0.1,0.1), labels = comma)+
  scale_x_discrete(expand = c(0.1,0.1))+
  facet_grid(~ancestry, scale="free", space="free_x")
C
D <-ggplot(data = ggroc_females, aes(x=adjusted_specificity, y=sensitivity))+
  geom_line(aes(colour=name))+theme_light(base_size = 12)+
  scale_colour_manual(values=c("#7ACCD7", "#B865A1","#7C6C65"))+
  labs(x = "1-Specificity", y = "Sensitivity", colour="AUROC model")+
  theme(legend.position = "bottom", legend.margin = margin(t = 0, r = 0, b = 0, l = 0), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.x = element_text(margin = margin(t = 10, b=0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+    
  facet_wrap(~ancestry)
D
#strip.background =element_rect(fill="Black"),

png(filename = "/Volumes//biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/PRS/female_GRS_all.png", width = 7, height = 12, units = "in", res = 300)
grid.arrange(A + 
               theme(plot.margin = unit(c(0,0.25,0,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "A"),
             B + 
               theme(plot.margin = unit(c(0,0.25,0,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "B"),
             C + 
               theme(plot.margin = unit(c(0,0.25,0,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "C"),
             D +
               theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"), plot.tag = element_text()) + 
               labs(tag = "D"),
             layout_matrix = matrix( c(1,2,3,4), ncol = 1), heights = c(5,5,7,6))
dev.off()

