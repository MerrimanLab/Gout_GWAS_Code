#What did I do?
  
  #find regions of significance - ie reducing to loci

### reducing snps into loci (note changes to file names and locations may have occured)

#if (!require("BiocManager"))
install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#BiocManager::install("RIPSeeker")
library(GenomicRanges)
library(tidyverse)

poplist<-c("AFR","EAS","EUR","LAT")
sexlist<-c("female", "male")#full
for (s in sexlist){
for( POP in poplist){
  filename<-paste0("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/",POP,"/", s, "/",POP,"_meta_",s,"1_clean_rsid.nfiltered.biallelic.txt")
  GWAS<-read.delim(filename, stringsAsFactors = FALSE, header = TRUE)
  GWAS_2<-GWAS %>% filter(P<1e-7) %>%  mutate(STARTwide = BP - 50000, STOPwide= BP + 50000)
  padded<-GRanges(seqnames = GWAS_2$CHR, ranges = IRanges(start = GWAS_2$STARTwide, end = GWAS_2$STOPwide, names = GWAS_2$SNP))
  padded_reduced<-GenomicRanges::reduce(padded)
  GWAS_gr <- GRanges(seqnames = GWAS_2$CHR, ranges = IRanges(start = GWAS_2$BP, width = 1), names = GWAS_2$SNP)
  reduced <- padded_reduced %>% as.data.frame()
  reduced$SNP_number<-countOverlaps(padded_reduced, GWAS_gr)
  reduced$locus<-paste("chr", reduced$seqnames,"_", round(reduced$start/1000000, 2),"_",round(reduced$end/1000000, 2),"MB", sep="")
  write.table(reduced, paste0("/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/1_reduced_loci/",s,"_",POP,"GWAS_SNP_reduced_loci.txt"), sep = "\t", na="", quote=F, row.names = F)
}
}


GWAS<-read.delim("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/AFR/female/AFR_meta_full1_clean_rsid.nfiltered.biallelic.txt")

sexlist<-c("full","female", "male")
for (s in sexlist){
TAMA_GWAS_2021<-read.delim("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_meta_analysis/TAMA/", s,"/tama_",s,"_clean.nfiltered.biallelic.txt", stringsAsFactors = FALSE, header = TRUE)
TAMA_GWAS_2021_2<-TAMA_GWAS_2021 %>% filter(logBF>5) %>%  mutate(STARTwide = POS - 50000, STOPwide= POS + 50000)
TAMA_padded<-GRanges(seqnames = TAMA_GWAS_2021_2$CHR, ranges = IRanges(start = TAMA_GWAS_2021_2$STARTwide, end = TAMA_GWAS_2021_2$STOPwide, names = TAMA_GWAS_2021_2$gene_name))
TAMA_padded_reduced<-GenomicRanges::reduce(TAMA_padded)
TAMA_GWAS_gr <- GRanges(seqnames = TAMA_GWAS_2021_2$CHR, ranges = IRanges(start = TAMA_GWAS_2021_2$POS, width = 1), names = TAMA_GWAS_2021_2$SNP)
TAMA_reduced <- TAMA_padded_reduced %>% as.data.frame()
TAMA_reduced$SNP_number<-countOverlaps(TAMA_padded_reduced, TAMA_GWAS_gr)
TAMA_reduced$locus<-paste("chr", TAMA_reduced$seqnames,"_", round(TAMA_reduced$start/1000000, 2),"_",round(TAMA_reduced$end/1000000, 2),"MB", sep="")
write.table(TAMA_reduced, "/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS results/reduced_loci/", s, "_TAMA_GWAS_SNP_reduced_loci.txt", sep = "\t", na="", quote=F, row.names = F)
}

### These results went on to Tanya for clumping on the regions of significance
### from the results of clumping I merged back in the GWAS results to help with the decision making

poplist<-c("AFR","EAS","EUR","LAT")
sexlist<-c("full","female", "male")
for (s in sexlist){
for(POP in poplist){
  results_full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/data/summarystats/",POP,"/",POP,"_meta_", s,"1_clean_rsid.nfiltered.biallelic.txt"), stringsAsFactors = FALSE, header = TRUE)
  full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/results/loci_clumping/",POP,"/",POP,"_", s,"_sigregion_0.01clumping_resultsummary.txt"))
  full<-merge(full, results_full, by = "SNP", all.x = T)
  full$result<-NA
  full$result[is.na(full$signals)]="sub significant"
  full$result[which(full$signals==1 & full$SNP_number==1)]="singleton"
  full$result[which(full$signals==1 & full$SNP_number>1)]="GOOD"
  write.table(full, paste0("/Volumes/userdata/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS results/",POP,"_", s,"_sigregion_0.01clumping_resultsummary_anno.txt"), row.names = F, sep = "\t", na="", quote=F)
}
}

sexlist<-c("full","female", "male")
for (s in sexlist){
  for(POP in poplist){
    results_full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/data/summarystats/TAMA/tama_meta_", s,"_clean_rsid.nfiltered.biallelic.txt"), stringsAsFactors = FALSE, header = TRUE)
    full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/results/loci_clumping/TAMA/TAMA_", s,"_sigregion_0.01clumping_resultsummary.txt"))
    full<-merge(full, results_full, by = "SNP", all.x = T)
    full$result<-NA
    full$result[is.na(full$signals)]="sub significant"
    full$result[which(full$signals==1 & full$SNP_number==1)]="singleton"
    full$result[which(full$signals==1 & full$SNP_number>1)]="GOOD"
    write.table(full, paste0("/Volumes/userdata/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS results/TAMA_", s,"_sigregion_0.01clumping_resultsummary_anno.txt"), row.names = F, sep = "\t", na="", quote=F)
  }
}

sexlist<-c("full","female", "male")
for (s in sexlist){
  for(POP in poplist){
    results_full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/data/summarystats/LAT/LAT_meta_", s,"_clean_rsid.nfiltered.biallelic.txt"), stringsAsFactors = FALSE, header = TRUE)
    full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/results/loci_clumping/LAT/LAT_", s,"_TAMA-LD_sigregion_0.01clumping_resultsummary.txt"))
    full<-merge(full, results_full, by = "SNP", all.x = T)
    full$result<-NA
    full$result[is.na(full$signals)]="sub significant"
    full$result[which(full$signals==1 & full$SNP_number==1)]="singleton"
    full$result[which(full$signals==1 & full$SNP_number>1)]="GOOD"
    write.table(full, paste0("/Volumes/userdata/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS results/LAT_", s,"_TAMA-LD_sigregion_0.01clumping_resultsummary_anno.txt"), row.names = F, sep = "\t", na="", quote=F)
  }
}


##### MAGMA ####

#running magma is outlined in Ruth/2021_projects/post_GWAS/analysis/Magma/Magma_results_MajorGwas.Rmd
#merge in magma to loci file. 


#### Plotting ###
#Tanya does it for the paper. This is just code for me
####
sexlist<-c("full","female", "male")
for (s in sexlist){
for(POP in poplist){
  results_full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/data/summarystats",POP,"/",POP,"_",s,"_clean.nfiltered.biallelic.txt"), stringsAsFactors = FALSE, header = TRUE)  
  full<-read.delim(paste0("/Volumes/scratch/merrimanlab/tanya/GWAS_Paper/PostGWAS_Repo/results/loci_clumping/",POP,"/",POP,"_",s,"_sigregion_0.01clumping_resultsummary.txt"))
  full[full==""]=NA
  full<-full[is.na(full$SNP)==F,]
  colnames(results_full)[3]<-"BP"
  for (r in 1:nrow(full)){
    rown<-full[r,]
    locus.zoom(data = results_full,                                    
               region = c(rown$seqnames, rown$start, rown$end),
               snp = rown$SNP,
               ignore.lead = T,
               secondary.label = T,
               offset_bp = 0,                                                                                   
               genes.data = Unique.genes,                  
               plot.title = paste0(POP,"_full_",rown$locus, rown$SNP),     
               file.name = paste0(POP,"/",POP,"_full_",rown$locus,"_", rown$SNP,".jpeg"),
               population = POP,
               nominal = 4,
               significant = 6,
               sig.type = "BF", 
               plot.type = "jpg",
               colour.genes = FALSE
    )  
  }
}
}




### Make tables ###
#check LD between things - Tanya code


phenoscanner
#Run SNPs through phenoscanner
#install_github("phenoscanner/phenoscanner")
sexlist=c("full","male","female")
for (s in sexlist){

all_list<-read.delim(paste0("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/",s,"_indepSNP_summary_1jun2022.txt"))
d1<-split(all_list$SNP,ceiling(seq_along(all_list$SNP)/100))

library(phenoscanner)
total<-as.data.frame(matrix(ncol=22,nrow=0))
colnames(total)<-c("snp","rsid","hg19_coordinates","hg38_coordinates","a1","a2","trait","efo","study","pmid","ancestry","year","beta","se","p","direction","n","n_cases","n_controls","n_studies","unit","dataset")

for(d in seq(1,length(d1))){
res <- phenoscanner(snpquery=d1[[d]], pvalue = 0.00000005)
total<-rbind(total,res$results)
}
write.table(total, paste0("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/",s,"_phenoscanner_raw.txt"), quote = F, sep = "\t", na="", row.names = F)
}
### add in phenoscanner results ###
sexlist=c("full","male","female")
for (s in sexlist){
pheno<-read.delim(paste0("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/",s,"_phenoscanner_raw.txt"))
library(tidyverse)
pheno<-pheno %>% dplyr::select(snp,trait,study,pmid)
pheno2<-unique(pheno)
pheno2$combined<-NA
pheno2$combined<-paste(pheno2$trait, pheno2$study, pheno2$pmid, sep = ", ")
pheno2<-pheno2[,c(1,5)]
pheno2$combined2<-NA
pheno2$combined2<-paste(pheno2$snp, pheno2$combined, sep = ", ")

urate_gout<-pheno2[str_detect(pheno2$combined, "urate|uric|gout|Urate|Uric|Gout|allopurinol|Allopurinol"),]
phe2<-pheno2[grepl(pattern = "Hashimoto|Autoimmune thyroid disease|Type 1 diabetes|Addison|Adrenal insufficiency|Pernicious anemia|Acute disseminated encephalomyelitis|Myasthenia gravis|Uveitis|Goodpasture|Anti-Neutrophil Cytoplasmic Antibody associated vasculitis|Microscopic polyangiitis|Wegener|Giant cell arteritis|Takayasu|Autoimmune gastritis|Chronic atrophic gastritis|Crohn|Ulcerative colitis|Coeliac disease|Celiac disease|Primary biliary cholangitis|Primary biliary cirrhosis|Vitiligo|Pemphigus|Pemphigoid|Erythema nodosum|Psoriasis|Psoriatic arthritis|Sarcoidosis|Sjogren|Sicca syndrome|Scleroderma|Dermatomyositis|Polymyositis|Behcet|Rheumatoid arthritis|Ankylosing spondylitis|Basophil|Eosinophil|Granulocyte|Lymphocyte|Monocyte|Neutrophil|White blood cell|White cell|basophil|eosinophil|granulocyte|lymphocyte|monocyte|neutrophil|white blood cell|white cell", x= pheno2$combined),]
phe3<-pheno2[grepl(pattern = "hashimoto|autoimmune thyroid disease|type 1 diabetes|addison|adrenal insufficiency|pernicious anemia|acute disseminated encephalomyelitis|myasthenia gravis|uveitis|goodpasture|anti-neutrophil cytoplasmic antibody associated vasculitis|microscopic polyangiitis|wegener|giant cell arteritis|takayasu|autoimmune gastritis|chronic atrophic gastritis|crohn|ulcerative colitis|coeliac disease|celiac disease|primary biliary cholangitis|primary biliary cirrhosis|vitiligo|pemphigus|pemphigoid|erythema nodosum|psoriasis|psoriatic arthritis|sarcoidosis|sjogren|sicca syndrome|scleroderma|dermatomyositis|polymyositis|behcet|rheumatoid arthritis|ankylosing spondylitis", x= pheno2$combined),]
inflammatory_autoimmune<-unique(rbind(phe2,phe3))

pheno3<-pheno2[which(!pheno2$combined2 %in% urate_gout$combined2),]
pheno3<-pheno3[which(!pheno3$combined2 %in% inflammatory_autoimmune$combined2),]
urate_gout <-urate_gout %>% group_by(snp) %>% summarise(values = paste(combined, collapse =  "; "))
colnames(urate_gout)<-c("SNP","phenoscanner_urate_gout")
pheno3 <-pheno3 %>% group_by(snp) %>% summarise(values = paste(combined, collapse =  "; "))
colnames(pheno3)<-c("SNP","phenoscanner_other_phenotypes")
inflammatory_autoimmune <-inflammatory_autoimmune  %>% group_by(snp) %>% summarise(values = paste(combined, collapse =  "; "))
colnames(inflammatory_autoimmune )<-c("SNP","phenoscanner_inflammatory_autoimmune")

snps_all<-read.delim(paste0("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/",s,"_indepSNP_summary_1jun2022.txt"))
snps_all<-merge(snps_all,urate_gout, by="SNP", all.x=T)
snps_all<-merge(snps_all,inflammatory_autoimmune, by="SNP", all.x=T)
snps_all<-merge(snps_all,pheno3, by="SNP", all.x=T)
write.table(snps_all, paste0("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/",s,"_indepSNP_summary_1jun2022_phenoscanner.txt"), quote = F, sep = "\t", na="", row.names = F)
}



#Run SNPs through phenoscanner
#install_github("phenoscanner/phenoscanner")

  all_list<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/gene_prioritisation.unique.txt")
  protein<-all_list[which(all_list$coding=="proteincoding"),]
  psuedogene<-all_list[which(all_list$coding=="psuedogene"),]
  rna<-all_list[which(all_list$coding=="lncRNA"|all_list$coding=="lncRNA"),]
  
  genelist<-unique(protein$gene)
  genelist<-unique(psuedogene$gene)
  genelist<-unique(rna$gene)
  
  d1<-split(genelist,ceiling(seq_along(genelist)/10))
  
  library(phenoscanner)
  total<-as.data.frame(matrix(ncol=22,nrow=0))
  colnames(total)<-c("snp","rsid","hg19_coordinates","hg38_coordinates","a1","a2","trait","efo","study","pmid","ancestry","year","beta","se","p","direction","n","n_cases","n_controls","n_studies","unit","dataset")
  
  for(d in seq(1,length(d1))){
    res <- phenoscanner(genequery = d1[[d]], pvalue = 0.00000005)
    total<-rbind(total,res$results)
  }
  write.table(total,"/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/genes_phenoscanner_raw_protien.txt", quote = F, sep = "\t", na="", row.names = F)
  write.table(total,"/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/genes_phenoscanner_raw_psuedogene.txt", quote = F, sep = "\t", na="", row.names = F)
  write.table(total,"/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/genes_phenoscanner_raw_rna.txt", quote = F, sep = "\t", na="", row.names = F)
  totalrna<-total
  totalprotein<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/genes_phenoscanner_raw_protien.txt")
  totalpsuedogene<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/genes_phenoscanner_raw_psuedogene.txt")
  
  total<-rbind(totalprotein, totalpsuedogene, totalrna)
  total<-unique(total)
 check <-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/genes_phenoscanner_raw_part1.txt")
  genelist2<- genelist[!genelist %in% total$gene ]

  write.table(total,"/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/genes_phenoscanner_raw.txt", quote = F, sep = "\t", na="", row.names = F)
  
  # Filter for white blood cell traits:
  keywords = c('cell', 'count', 'cyte', 'phil')
  pattern = paste(keywords, collapse = '|')
  
  # The second filter is to remove red blood cell traits and other unwanted traits
  phewas = total %>% filter(grepl(pattern, trait, ignore.case = T)) %>% filter(!grepl('Red|scatter percentage|expression|cancer|carcinoma|death', trait)) 
  phewas$n_cases[which(phewas$n_cases==0)]=NA
  phewas$n_cases[which(phewas$n_cases=="-")]=NA
  
  phewas<-phewas[is.na(phewas$n_cases),]
  blood_cell_genes<-unique(phewas$gene)
  
  immune = c('cell', 'count', 'cyte', 'phil')
  
  
  
  total$n_cases[which(total$n_cases==0)]=NA
  total$n_cases[which(total$n_cases=="-")]=NA
  
  total$n_cases<-as.numeric(total$n_cases)
  total2<-total[which(total$n_cases>=100|is.na(total$n_cases)),]
  
  autoimmune<-c("alopecia areata", "ankylosing", "autoimmune", "behcets", "celiac", "inflammatory polyneuropathy", "crohn", "type 1 diabetes", "glomerulonephritis", "granulomatous?", "graves", "guillain", "thyroiditis", "iga nephropathy", "idiopathic arthritis", "menieres ", "multiple sclerosis", "myasthenia gravis", "pemphigus vulgaris", "polymyalgia", "primary biliary", "psoriasis", "psoriatic", "raynauds", "rheumatic fever", "rheumatoid", "sarcoidosis", "scleroderma", "sjogrens", "lupus erythematosus", "ulcerative", "vitiligo", "coeliac", "diabetes mellitus type 1", "malabsorption")
  pattern = paste(autoimmune, collapse = '|')
  autoimmune_associations<-total2 %>% filter(grepl(pattern, trait))
  
  autoimmune_associations_table<-as.data.frame(table(autoimmune_associations$trait,autoimmune_associations$gene))
  autoimmune_associations_table<-autoimmune_associations_table[which(autoimmune_associations_table$Freq>0),]
  
  ### add in phenoscanner results ###
sexlist=c("full","male","female")
for (s in sexlist){
  pheno<-read.delim(paste0("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/",s,"_phenoscanner_raw.txt"))
  library(tidyverse)
  pheno<-pheno %>% dplyr::select(snp,trait,study,pmid)
  pheno2<-unique(pheno)
  pheno2$combined<-NA
  pheno2$combined<-paste(pheno2$trait, pheno2$study, pheno2$pmid, sep = ", ")
  pheno2<-pheno2[,c(1,5)]
  pheno2$combined2<-NA
  pheno2$combined2<-paste(pheno2$snp, pheno2$combined, sep = ", ")
  
  urate_gout<-pheno2[str_detect(pheno2$combined, "urate|uric|gout|Urate|Uric|Gout|allopurinol|Allopurinol"),]
  phe2<-pheno2[grepl(pattern = "Hashimoto|Autoimmune thyroid disease|Type 1 diabetes|Addison|Adrenal insufficiency|Pernicious anemia|Acute disseminated encephalomyelitis|Myasthenia gravis|Uveitis|Goodpasture|Anti-Neutrophil Cytoplasmic Antibody associated vasculitis|Microscopic polyangiitis|Wegener|Giant cell arteritis|Takayasu|Autoimmune gastritis|Chronic atrophic gastritis|Crohn|Ulcerative colitis|Coeliac disease|Celiac disease|Primary biliary cholangitis|Primary biliary cirrhosis|Vitiligo|Pemphigus|Pemphigoid|Erythema nodosum|Psoriasis|Psoriatic arthritis|Sarcoidosis|Sjogren|Sicca syndrome|Scleroderma|Dermatomyositis|Polymyositis|Behcet|Rheumatoid arthritis|Ankylosing spondylitis|Basophil|Eosinophil|Granulocyte|Lymphocyte|Monocyte|Neutrophil|White blood cell|White cell|basophil|eosinophil|granulocyte|lymphocyte|monocyte|neutrophil|white blood cell|white cell", x= pheno2$combined),]
  phe3<-pheno2[grepl(pattern = "hashimoto|autoimmune thyroid disease|type 1 diabetes|addison|adrenal insufficiency|pernicious anemia|acute disseminated encephalomyelitis|myasthenia gravis|uveitis|goodpasture|anti-neutrophil cytoplasmic antibody associated vasculitis|microscopic polyangiitis|wegener|giant cell arteritis|takayasu|autoimmune gastritis|chronic atrophic gastritis|crohn|ulcerative colitis|coeliac disease|celiac disease|primary biliary cholangitis|primary biliary cirrhosis|vitiligo|pemphigus|pemphigoid|erythema nodosum|psoriasis|psoriatic arthritis|sarcoidosis|sjogren|sicca syndrome|scleroderma|dermatomyositis|polymyositis|behcet|rheumatoid arthritis|ankylosing spondylitis", x= pheno2$combined),]
  inflammatory_autoimmune<-unique(rbind(phe2,phe3))
  
  pheno3<-pheno2[which(!pheno2$combined2 %in% urate_gout$combined2),]
  pheno3<-pheno3[which(!pheno3$combined2 %in% inflammatory_autoimmune$combined2),]
  urate_gout <-urate_gout %>% group_by(snp) %>% summarise(values = paste(combined, collapse =  "; "))
  colnames(urate_gout)<-c("SNP","phenoscanner_urate_gout")
  pheno3 <-pheno3 %>% group_by(snp) %>% summarise(values = paste(combined, collapse =  "; "))
  colnames(pheno3)<-c("SNP","phenoscanner_other_phenotypes")
  inflammatory_autoimmune <-inflammatory_autoimmune  %>% group_by(snp) %>% summarise(values = paste(combined, collapse =  "; "))
  colnames(inflammatory_autoimmune )<-c("SNP","phenoscanner_inflammatory_autoimmune")
  
  snps_all<-read.delim(paste0("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/",s,"_indepSNP_summary_1jun2022.txt"))
  snps_all<-merge(snps_all,urate_gout, by="SNP", all.x=T)
  snps_all<-merge(snps_all,inflammatory_autoimmune, by="SNP", all.x=T)
  snps_all<-merge(snps_all,pheno3, by="SNP", all.x=T)
  write.table(snps_all, paste0("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/",s,"_indepSNP_summary_1jun2022_phenoscanner.txt"), quote = F, sep = "\t", na="", row.names = F)
}






#COJO
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")

library(GenomicRanges)
a<-list.files(path = "/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/COJO/results/", pattern ="_results.jma.cojo", full.names = T)
cojo<-map_dfr(a, read_tsv, col_types = cols(Chr = col_double(),  SNP = col_character(),  bp = col_double(),  refA = col_character(),  freq = col_double(),  b = col_double(),  se = col_double(),  p = col_double(), n = col_double(),  freq_geno = col_double(),  bJ = col_double(),  bJ_se = col_double(),  pJ = col_double(),  LD_r = col_double()))
write.table(cojo, "/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/COJO/results/All_european_cojo.txt", na="", row.names = F, quote = F, sep="\t")

a<-list.files(path = "/Volumes/scratch/merrimanlab/ruth/cojo_script/results/", pattern ="cojo_results_male.jma.cojo", full.names = T)
cojo<-map_dfr(a, read_tsv, col_types = cols(Chr = col_double(),  SNP = col_character(),  bp = col_double(),  refA = col_character(),  freq = col_double(),  b = col_double(),  se = col_double(),  p = col_double(), n = col_double(),  freq_geno = col_double(),  bJ = col_double(),  bJ_se = col_double(),  pJ = col_double(),  LD_r = col_double()))
write.table(cojo, "/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/COJO/results/All_european_cojo_male.txt", na="", row.names = F, quote = F, sep="\t")

a<-list.files(path = "/Volumes/scratch/merrimanlab/ruth/cojo_script/results/", pattern ="cojo_results_female.jma.cojo", full.names = T)
cojo<-map_dfr(a, read_tsv, col_types = cols(Chr = col_double(),  SNP = col_character(),  bp = col_double(),  refA = col_character(),  freq = col_double(),  b = col_double(),  se = col_double(),  p = col_double(), n = col_double(),  freq_geno = col_double(),  bJ = col_double(),  bJ_se = col_double(),  pJ = col_double(),  LD_r = col_double()))
write.table(cojo, "/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/COJO/results/All_european_cojo_female.txt", na="", row.names = F, quote = F, sep="\t")


EUR_cojo<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/COJO/results/All_european_cojo.txt")
loci<-read.delim("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/full_loci_summary_1jun2022.txt")
EUR_cojo[,"id"] <- row.names(EUR_cojo)
loci[,"id"]<- row.names(loci)
EUR_cojo_gr <- GRanges(seqnames = EUR_cojo$Chr, ranges = IRanges(start = EUR_cojo$bp, end = EUR_cojo$bp), SNP_name = EUR_cojo$SNP)
loci_all_gr <- GRanges(seqnames = loci$chr, ranges = IRanges(start = loci$start, end = loci$end), LOCI_name = loci$LOCI)
overlap_hits <- findOverlaps(loci_all_gr, EUR_cojo_gr)
combined_snp <- merge(merge(loci, overlap_hits, by.x = "id", by.y = "queryHits", all=T), EUR_cojo, by.x = "subjectHits", by.y = "id", all=T)
#remove extra cojo snp outside of loci
combined_snp<-combined_snp[is.na(combined_snp$locus)==FALSE,]
combined_snp2 <-combined_snp %>% arrange(p) %>% group_by(locus) %>% summarise(values = paste(SNP, collapse =  "; "))
colnames(combined_snp2)<-c("locus","EUR_COJO")
combined_snp2$EUR_COJO[which(combined_snp2$EUR_COJO=="NA")]=NA
loci<-loci[1:11]
loci<-merge(loci, combined_snp2, by="locus")
loci[loci==""]=NA
write.table(loci, "/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/full_loci_summary_1jun2022_cojo.txt", sep="\t", na="", row.names = F, quote=F)


EUR_cojo<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/COJO/results/All_european_cojo_female.txt")
loci<-read.delim("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/female_loci_summary_1jun2022.txt")
EUR_cojo[,"id"] <- row.names(EUR_cojo)
loci[,"id"]<- row.names(loci)
EUR_cojo_gr <- GRanges(seqnames = EUR_cojo$Chr, ranges = IRanges(start = EUR_cojo$bp, end = EUR_cojo$bp), SNP_name = EUR_cojo$SNP)
loci_all_gr <- GRanges(seqnames = loci$chr, ranges = IRanges(start = loci$start, end = loci$end), LOCI_name = loci$LOCI)
overlap_hits <- findOverlaps(loci_all_gr, EUR_cojo_gr)
combined_snp <- merge(merge(loci, overlap_hits, by.x = "id", by.y = "queryHits", all=T), EUR_cojo, by.x = "subjectHits", by.y = "id", all=T)
#remove extra cojo snp outside of loci
combined_snp<-combined_snp[is.na(combined_snp$locus)==FALSE,]
combined_snp2 <-combined_snp %>% arrange(p) %>% group_by(locus) %>% summarise(values = paste(SNP, collapse =  "; "))
colnames(combined_snp2)<-c("locus","EUR_COJO")
combined_snp2$EUR_COJO[which(combined_snp2$EUR_COJO=="NA")]=NA
loci<-loci[1:12]
loci<-merge(loci, combined_snp2, by="locus")
loci[loci==""]=NA
write.table(loci, "/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/male_loci_summary_1jun2022_cojo.txt", sep="\t", na="", row.names = F, quote=F)

EUR_cojo<-read.delim("/Volumes/biochemistry/Lab_Groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/COJO/results/All_european_cojo_male.txt")
loci<-read.delim("/Volumes/archive/merrimanlab/major_gwas_paper_archive/results_for_paper/loci_list/male_loci_summary_1jun2022.txt")
EUR_cojo[,"id"] <- row.names(EUR_cojo)
loci[,"id"]<- row.names(loci)
EUR_cojo_gr <- GRanges(seqnames = EUR_cojo$Chr, ranges = IRanges(start = EUR_cojo$bp, end = EUR_cojo$bp), SNP_name = EUR_cojo$SNP)
loci_all_gr <- GRanges(seqnames = loci$chr, ranges = IRanges(start = loci$start, end = loci$end), LOCI_name = loci$LOCI)
overlap_hits <- findOverlaps(loci_all_gr, EUR_cojo_gr)
combined_snp <- merge(merge(loci, overlap_hits, by.x = "id", by.y = "queryHits", all=T), EUR_cojo, by.x = "subjectHits", by.y = "id", all=T)
#remove extra cojo snp outside of loci
combined_snp<-combined_snp[is.na(combined_snp$locus)==FALSE,]
combined_snp2 <-combined_snp %>% arrange(p) %>% group_by(locus) %>% summarise(values = paste(SNP, collapse =  "; "))
colnames(combined_snp2)<-c("locus","EUR_COJO")
combined_snp2$EUR_COJO[which(combined_snp2$EUR_COJO=="NA")]=NA
loci<-loci[1:12]
loci<-merge(loci, combined_snp2, by="locus")
loci[loci==""]=NA
write.table(loci, "/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/male_loci_summary_1jun2022_cojo.txt", sep="\t", na="", row.names = F, quote=F)



# phecode PheWAS


# PRS
loci<-read.delim("/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/full_loci_summary_1jun2022_cojo.txt")
snps<-read.delim("/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/full_indepSNP_summary_1jun2022_phenoscanner.txt")
loci$snp_used<-NA
loci$snp_used[is.na(loci$EUR)==F]=loci$EUR[is.na(loci$EUR)==F]
loci[loci==""]=NA
loci$effect_from<-NA
loci$effect_from[is.na(loci$EUR)==F]="EUR"
loci[loci==""]=NA
loci$effect_from[which(is.na(loci$snp_used) & is.na(loci$TAMA)==F)]="TAMA"
loci[loci==""]=NA
loci$snp_used[which(is.na(loci$snp_used) & is.na(loci$TAMA)==F)]=loci$TAMA[which(is.na(loci$snp_used) & is.na(loci$TAMA)==F)]
loci<-merge(loci,snps[,c("SNP","MajorAllele","MinorAllele","OR_EUR" ,"OR_TAMA")], by.x="snp_used", by.y="SNP",all.x=T)
loci$effect_used_UKBB=NA
loci$effect_used_UKBB[which(loci$effect_from=="EUR")]=loci$OR_EUR[which(loci$effect_from=="EUR")]
loci$effect_used_UKBB[which(loci$effect_from=="TAMA")]=loci$OR_TAMA[which(loci$effect_from=="TAMA")]
loci$risk_allele<-NA
loci$risk_allele[which(loci$effect_used_UKBB>=1)]=loci$MinorAllele[which(loci$effect_used_UKBB>=1)]
loci$risk_allele[which(loci$effect_used_UKBB<=1)]=loci$MajorAllele[which(loci$effect_used_UKBB<=1)]
loci$oth_allele<-NA
loci$oth_allele[which(loci$effect_used_UKBB>=1)]=loci$MajorAllele[which(loci$effect_used_UKBB>=1)]
loci$oth_allele[which(loci$effect_used_UKBB<=1)]=loci$MinorAllele[which(loci$effect_used_UKBB<=1)]
loci$effect_used_UKBB2=NA
loci$effect_used_UKBB2[which(loci$effect_used_UKBB>=1)]=loci$effect_used_UKBB[which(loci$effect_used_UKBB>=1)]
loci$effect_used_UKBB2[which(loci$effect_used_UKBB<=1)]=1/loci$effect_used_UKBB[which(loci$effect_used_UKBB<=1)]


write.table(loci, "/Volumes/biochemistry/Lab_groups/merrimanlab/Merriman_Documents/Ruth/2021_projects/post_GWAS_analysis/GWAS_results/4_tables/full_loci_summary_1jun2022_forPRS.txt", sep="\t", na="", row.names = F, quote=F)

#then used PRS_prediction_EUROPEAN_2021.rmd


