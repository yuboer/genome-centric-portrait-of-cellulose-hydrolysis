####### to calculate frequencing of cazy modules co-occuring with each other
#######
path = as.character("D:/papers_in_progress/to finish/HKU_cellulolytic_genomes/BMC Genomics/Environmental microbiome_submission/Appendix files/")
setwd(path)

#check
getwd()

cazy_in_nc_list<-dir('dbCAN_results')
#CAZy harboring NCs' names, 
queryinfo_in_nc_list<-dir('gbk_info')
#all the extracted NCs' names
head(cazy_in_nc_list)
head(queryinfo_in_nc_list)
length(cazy_in_nc_list)
#[1] 3898
length(queryinfo_in_nc_list)
#[1] 5216

queryinfo_in_cazy_harboring_nc_list<-queryinfo_in_nc_list[queryinfo_in_nc_list %in% cazy_in_nc_list]

length(queryinfo_in_cazy_harboring_nc_list)
#[1] 3898

cazy_in_nc<-list()
for (i in 1:length(cazy_in_nc_list))
{
  tmp<-read.delim(file=paste("D:/papers_in_progress/to finish/HKU_cellulolytic_genomes/BMC Genomics/Environmental microbiome_submission/Appendix files/dbCAN_results/",cazy_in_nc_list[i],sep=""),sep='\t',header = F)
  #past the path ahead, since the working directory is not where the file locates
  colnames(tmp)<-c("cazy_module","HMM_length","Query_ID","Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  cazy_in_nc[[cazy_in_nc_list[i]]]<-tmp
  #cazy_in_nc$cazy_in_nc_list[i]<-temp
  #cazy_in_nc_list[i]=cazy_in_nc[[i]]
}


queryinfo_in_cazy_harboring_nc<-list()
for (i in 1:length(queryinfo_in_cazy_harboring_nc_list))
{
  tmp<-read.delim(file=paste('D:/papers_in_progress/to finish/HKU_cellulolytic_genomes/BMC Genomics/Environmental microbiome_submission/Appendix files/gbk_info/',queryinfo_in_cazy_harboring_nc_list[i],sep=''),sep='\t',header = F)
  colnames(tmp)<-c('ID', 'desc', 'product', 'locus', 'strand', 'start', 'end')
  queryinfo_in_cazy_harboring_nc[[queryinfo_in_cazy_harboring_nc_list[i]]]<-tmp
}


#######for each .gbk file, info on all the cazy modules indentified in genes and the gene info
complete_info_nc<-list()
for (i in 1:length(cazy_in_nc_list))
{
  tmp<-merge(cazy_in_nc[[cazy_in_nc_list[i]]],queryinfo_in_cazy_harboring_nc[[cazy_in_nc_list[i]]],by.x="Query_ID",by.y="ID",all.x=T,all.y=F)
  complete_info_nc[[cazy_in_nc_list[i]]]<-tmp
}

all_cazy_harboring_nc<-names(complete_info_nc)
######


###### CAZy module co-occuring with SLH module


all_cazy_in_slh_harboring_genes<-data.frame(Query_ID=c(NA),cazy_module=c(NA),HMM_length=c(NA),Query_length=c(NA),E_value=c(NA),HMM_start=c(NA),HMM_end=c(NA),
                                            Query_start=c(NA),Query_end=c(NA),Coverage=c(NA),desc=c(NA),product=c(NA),locus=c(NA),strand=c(NA),start=c(NA),end=c(NA))

for (i in 1:length(all_cazy_harboring_nc))
{
  slh_tmp<-complete_info_nc[[i]]
  slh_tmp$cazy_module<-gsub('.hmm','',slh_tmp$cazy_module)
  slh_harboring_genes<-slh_tmp$Query_ID[grep('SLH',slh_tmp$cazy_module)]
  slh_harboring_genes<-unique(slh_harboring_genes)
  all_cazy_in_slh_tmp<-slh_tmp[slh_tmp$Query_ID %in% slh_harboring_genes,]
  all_cazy_in_slh_harboring_genes<-rbind(all_cazy_in_slh_harboring_genes,all_cazy_in_slh_tmp)
}

colnames(all_cazy_in_slh_harboring_genes)

dim(all_cazy_in_slh_harboring_genes)
all_cazy_in_slh_harboring_genes<-all_cazy_in_slh_harboring_genes[-1,]


all_cazy_in_slh_harboring_genes_CBM<-all_cazy_in_slh_harboring_genes[grep('CBM',all_cazy_in_slh_harboring_genes$cazy_module),]
dim(all_cazy_in_slh_harboring_genes_CBM)
#[1] 504  16
all_slh_cbm_genes_id<-all_cazy_in_slh_harboring_genes_CBM[,'Query_ID']
#length(all_slh_cbm_genes_id)
#[1] 504
uniq_gene_id_slh_cbm<-unique(all_slh_cbm_genes_id)
#length(uniq_gene_id_slh_cbm)
#[1] 211
CBM_count_in_slh_genes<-as.data.frame(table(all_cazy_in_slh_harboring_genes_CBM$cazy_module))
colnames(CBM_count_in_slh_genes)<-c('CBM_module','Frequency')
write.table(CBM_count_in_slh_genes,file="CBM_count_in_slh_genes.tab",sep='\t',row.names = F, col.names = T)
cazy_count_in_slh_genes<-as.data.frame(table(all_cazy_in_slh_harboring_genes$cazy_module))
colnames(cazy_count_in_slh_genes)<-c('cazy_module','Frequency')
write.table(cazy_count_in_slh_genes,file="cazy_count_in_slh_genes.tab",sep='\t',row.names = F, col.names = T)

#
slh_harboring_gens<-all_cazy_in_slh_harboring_genes$Query_ID
uniq_slh_harboing_gene_id<-unique(slh_harboring_gens)
length(uniq_slh_harboing_gene_id)
#[1] 2757
cazy_in_slh_harboring_genes<-data.frame(cazy_module=c(NA))
for (i in 1:length(uniq_slh_harboing_gene_id))
  {
    tmp_cazy_slh<-all_cazy_in_slh_harboring_genes[all_cazy_in_slh_harboring_genes$Query_ID %in% uniq_slh_harboing_gene_id[i],]
    cazy_abun_in_slh_gene<-as.data.frame(table(tmp_cazy_slh$cazy_module))
    colnames(cazy_abun_in_slh_gene)<-c("cazy_module",uniq_slh_harboing_gene_id[i])
    cazy_in_slh_harboring_genes<-merge(cazy_in_slh_harboring_genes,cazy_abun_in_slh_gene,by.x = 'cazy_module',by.y = 'cazy_module',all=T)
}
dim(cazy_in_slh_harboring_genes)
#[1]   76 2758
cazy_in_slh_harboring_genes<-cazy_in_slh_harboring_genes[!is.na(cazy_in_slh_harboring_genes$cazy_module),]
rownames(cazy_in_slh_harboring_genes)<-cazy_in_slh_harboring_genes$cazy_module
cazy_in_slh_harboring_genes<-cazy_in_slh_harboring_genes[,-1]
cazy_in_slh_harboring_genes_t<-t(cazy_in_slh_harboring_genes)
write.table(cazy_in_slh_harboring_genes_t,file='cazy_in_slh_harboring_genes.tab',sep='\t',row.names = T, col.names =T)
######

all_cazy_in_cohesin_harboring_genes<-data.frame(Query_ID=c(NA),cazy_module=c(NA),HMM_length=c(NA),Query_length=c(NA),E_value=c(NA),HMM_start=c(NA),HMM_end=c(NA),
                                                Query_start=c(NA),Query_end=c(NA),Coverage=c(NA),desc=c(NA),product=c(NA),locus=c(NA),strand=c(NA),start=c(NA),end=c(NA))

for (i in 1:length(all_cazy_harboring_nc))
{
  cohesin_tmp<-complete_info_nc[[i]]
  cohesin_tmp$cazy_module<-gsub('.hmm','',cohesin_tmp$cazy_module)
  cohesin_harboring_genes<-cohesin_tmp$Query_ID[grep('cohesin',cohesin_tmp$cazy_module)]
  all_cazy_in_cohesin_tmp<-cohesin_tmp[cohesin_tmp$Query_ID %in% cohesin_harboring_genes,]
  all_cazy_in_cohesin_harboring_genes<-rbind(all_cazy_in_cohesin_harboring_genes,all_cazy_in_cohesin_tmp)
}
dim(all_cazy_in_cohesin_harboring_genes)
all_cazy_in_cohesin_harboring_genes<-all_cazy_in_cohesin_harboring_genes[-1,]
#dim(all_cazy_in_cohesin_harboring_genes)
#[1] 6696   16

all_cazy_in_cohesin_harboring_genes_CBM<-all_cazy_in_cohesin_harboring_genes[grep('CBM',all_cazy_in_cohesin_harboring_genes$cazy_module),]
dim(all_cazy_in_cohesin_harboring_genes_CBM)
#[1] 504  16
all_cohesin_cbm_genes_id<-all_cazy_in_cohesin_harboring_genes_CBM[,'Query_ID']
#length(all_cohesin_cbm_genes_id)
#[1] 504
uniq_gene_id_cohesin_cbm<-unique(all_cohesin_cbm_genes_id)
#length(uniq_gene_id_cohesin_cbm)
#[1] 211
CBM_count_in_cohesin_genes<-as.data.frame(table(all_cazy_in_cohesin_harboring_genes_CBM$cazy_module))
colnames(CBM_count_in_cohesin_genes)<-c('CBM_module','Frequency')
write.table(CBM_count_in_cohesin_genes,file="CBM_count_in_cohesin_genes.tab",sep='\t',row.names = F, col.names = T)
cazy_count_in_cohesin_genes<-as.data.frame(table(all_cazy_in_cohesin_harboring_genes$cazy_module))
colnames(cazy_count_in_cohesin_genes)<-c('cazy_module','Frequency')
write.table(cazy_count_in_cohesin_genes,file="cazy_count_in_cohesin_genes.tab",sep='\t',row.names = F, col.names = T)

#
cohesin_harboring_gens<-all_cazy_in_cohesin_harboring_genes$Query_ID
uniq_cohesin_harboing_gene_id<-unique(cohesin_harboring_gens)
length(uniq_cohesin_harboing_gene_id)
#[1] 2757
cazy_in_cohesin_harboring_genes<-data.frame(cazy_module=c(NA))
for (i in 1:length(uniq_cohesin_harboing_gene_id))
{
  tmp_cazy_cohesin<-all_cazy_in_cohesin_harboring_genes[all_cazy_in_cohesin_harboring_genes$Query_ID %in% uniq_cohesin_harboing_gene_id[i],]
  cazy_abun_in_cohesin_gene<-as.data.frame(table(tmp_cazy_cohesin$cazy_module))
  colnames(cazy_abun_in_cohesin_gene)<-c("cazy_module",uniq_cohesin_harboing_gene_id[i])
  cazy_in_cohesin_harboring_genes<-merge(cazy_in_cohesin_harboring_genes,cazy_abun_in_cohesin_gene,by.x = 'cazy_module',by.y = 'cazy_module',all=T)
}
dim(cazy_in_cohesin_harboring_genes)
#[1]   76 2758
cazy_in_cohesin_harboring_genes<-cazy_in_cohesin_harboring_genes[!is.na(cazy_in_cohesin_harboring_genes$cazy_module),]
rownames(cazy_in_cohesin_harboring_genes)<-cazy_in_cohesin_harboring_genes$cazy_module
cazy_in_cohesin_harboring_genes<-cazy_in_cohesin_harboring_genes[,-1]
cazy_in_cohesin_harboring_genes_t<-t(cazy_in_cohesin_harboring_genes)
write.table(cazy_in_cohesin_harboring_genes_t,file='cazy_in_cohesin_harboring_genes.tab',sep='\t',row.names = T, col.names =T)
######



all_cazy_in_dockerin_harboring_genes<-data.frame(Query_ID=c(NA),cazy_module=c(NA),HMM_length=c(NA),Query_length=c(NA),E_value=c(NA),HMM_start=c(NA),HMM_end=c(NA),
                                                 Query_start=c(NA),Query_end=c(NA),Coverage=c(NA),desc=c(NA),product=c(NA),locus=c(NA),strand=c(NA),start=c(NA),end=c(NA))

for (i in 1:length(all_cazy_harboring_nc))
{
  dockerin_tmp<-complete_info_nc[[i]]
  dockerin_tmp$cazy_module<-gsub('.hmm','',dockerin_tmp$cazy_module)
  dockerin_harboring_genes<-dockerin_tmp$Query_ID[grep('dockerin',dockerin_tmp$cazy_module)]
  dockerin_harboring_genes<-unique(dockerin_harboring_genes)
  all_cazy_in_dockerin_tmp<-dockerin_tmp[dockerin_tmp$Query_ID %in% dockerin_harboring_genes,]
  all_cazy_in_dockerin_harboring_genes<-rbind(all_cazy_in_dockerin_harboring_genes,all_cazy_in_dockerin_tmp)
}
dim(all_cazy_in_dockerin_harboring_genes)
all_cazy_in_dockerin_harboring_genes<-all_cazy_in_dockerin_harboring_genes[-1,]
#dim(all_cazy_in_dockerin_harboring_genes)
#[1] 6696   16

all_cazy_in_dockerin_harboring_genes_CBM<-all_cazy_in_dockerin_harboring_genes[grep('CBM',all_cazy_in_dockerin_harboring_genes$cazy_module),]
dim(all_cazy_in_dockerin_harboring_genes_CBM)
#[1] 504  16
all_dockerin_cbm_genes_id<-all_cazy_in_dockerin_harboring_genes_CBM[,'Query_ID']
#length(all_dockerin_cbm_genes_id)
#[1] 504
uniq_gene_id_dockerin_cbm<-unique(all_dockerin_cbm_genes_id)
#length(uniq_gene_id_dockerin_cbm)
#[1] 151
CBM_count_in_dockerin_genes<-as.data.frame(table(all_cazy_in_dockerin_harboring_genes_CBM$cazy_module))
colnames(CBM_count_in_dockerin_genes)<-c('CBM_module','Frequency')
write.table(CBM_count_in_dockerin_genes,file="CBM_count_in_dockerin_genes.tab",sep='\t',row.names = F, col.names = T)
cazy_count_in_dockerin_genes<-as.data.frame(table(all_cazy_in_dockerin_harboring_genes$cazy_module))
colnames(cazy_count_in_dockerin_genes)<-c('cazy_module','Frequency')
write.table(cazy_count_in_dockerin_genes,file="cazy_count_in_dockerin_genes.tab",sep='\t',row.names = F, col.names = T)

#
dockerin_harboring_gens<-all_cazy_in_dockerin_harboring_genes$Query_ID
uniq_dockerin_harboing_gene_id<-unique(dockerin_harboring_gens)
length(uniq_dockerin_harboing_gene_id)
#[1] 494
cazy_in_dockerin_harboring_genes<-data.frame(cazy_module=c(NA))
for (i in 1:length(uniq_dockerin_harboing_gene_id))
{
  tmp_cazy_dockerin<-all_cazy_in_dockerin_harboring_genes[all_cazy_in_dockerin_harboring_genes$Query_ID %in% uniq_dockerin_harboing_gene_id[i],]
  cazy_abun_in_dockerin_gene<-as.data.frame(table(tmp_cazy_dockerin$cazy_module))
  colnames(cazy_abun_in_dockerin_gene)<-c("cazy_module",uniq_dockerin_harboing_gene_id[i])
  cazy_in_dockerin_harboring_genes<-merge(cazy_in_dockerin_harboring_genes,cazy_abun_in_dockerin_gene,by.x = 'cazy_module',by.y = 'cazy_module',all=T)
}
dim(cazy_in_dockerin_harboring_genes)
#[1]   76 2758
cazy_in_dockerin_harboring_genes<-cazy_in_dockerin_harboring_genes[!is.na(cazy_in_dockerin_harboring_genes$cazy_module),]
rownames(cazy_in_dockerin_harboring_genes)<-cazy_in_dockerin_harboring_genes$cazy_module
cazy_in_dockerin_harboring_genes<-cazy_in_dockerin_harboring_genes[,-1]
cazy_in_dockerin_harboring_genes_t<-t(cazy_in_dockerin_harboring_genes)
write.table(cazy_in_dockerin_harboring_genes_t,file='cazy_in_dockerin_harboring_genes.tab',sep='\t',row.names = T, col.names =T)




#####
exo_cellulase_modules<-c('GH6','GH9','GH48')

endo_cellulase_modules<-c('GH5','GH8','GH12','GH44','GH45','GH51','GH74','GH124')
not_identified_endo_cellulase<-c('GH7','GH61')

hemicellulase<-c('GH10','GH11','GH16','GH26','GH28','GH53','GH74','GH113')
not_identified_hemicellulase_module<-c('GH134')


oligosaccharides<-c('GH1','GH2','GH3','GH29','GH30','GH35','GH38','GH39','GH42','GH43','GH47','GH52','GH59','GH92','GH94','GH116','GH120','GH125')
#ß_glucosidase<-c('GH1','GH2','GH3','GH30')
LPMO<-c('AA9','AA10','AA11')


for (j in 1:length(exo_cellulase_modules))
{
  all_cazy_in_harboring_genes<-data.frame(Query_ID=c(NA),cazy_module=c(NA),HMM_length=c(NA),Query_length=c(NA),E_value=c(NA),HMM_start=c(NA),HMM_end=c(NA),
                                          Query_start=c(NA),Query_end=c(NA),Coverage=c(NA),desc=c(NA),product=c(NA),locus=c(NA),strand=c(NA),start=c(NA),end=c(NA))
  for (i in 1:length(all_cazy_harboring_nc))
  {
    tmp<-complete_info_nc[[i]]
    tmp$cazy_module<-gsub('.hmm','',tmp$cazy_module)
    harboring_genes<-tmp$Query_ID[(tmp$cazy_module %in% exo_cellulase_modules[j])]
    harboring_genes<-unique(harboring_genes)
    all_cazy_in_tmp<-tmp[tmp$Query_ID %in% harboring_genes,]
    all_cazy_in_harboring_genes<-rbind(all_cazy_in_harboring_genes,all_cazy_in_tmp)
  }
  all_cazy_in_harboring_genes<-all_cazy_in_harboring_genes[-1,]
  cazy_count_in_genes<-as.data.frame(table(all_cazy_in_harboring_genes$cazy_module))
  colnames(cazy_count_in_genes)<-c('cazy_module','Frequency')
  write.table(cazy_count_in_genes,file=paste('cazy_count_in_',exo_cellulase_modules[j],'_genes.tab',sep = ''),sep='\t',row.names = F, col.names = T)
}
#####


######
for (j in 1:length(endo_cellulase_modules))
{
  all_cazy_in_harboring_genes<-data.frame(Query_ID=c(NA),cazy_module=c(NA),HMM_length=c(NA),Query_length=c(NA),E_value=c(NA),HMM_start=c(NA),HMM_end=c(NA),
                                          Query_start=c(NA),Query_end=c(NA),Coverage=c(NA),desc=c(NA),product=c(NA),locus=c(NA),strand=c(NA),start=c(NA),end=c(NA))
  for (i in 1:length(all_cazy_harboring_nc))
  {
    tmp<-complete_info_nc[[i]]
    tmp$cazy_module<-gsub('.hmm','',tmp$cazy_module)
    harboring_genes<-tmp$Query_ID[(tmp$cazy_module %in% endo_cellulase_modules[j])]
    harboring_genes<-unique(harboring_genes)
    all_cazy_in_tmp<-tmp[tmp$Query_ID %in% harboring_genes,]
    all_cazy_in_harboring_genes<-rbind(all_cazy_in_harboring_genes,all_cazy_in_tmp)
  }
  all_cazy_in_harboring_genes<-all_cazy_in_harboring_genes[-1,]
  cazy_count_in_genes<-as.data.frame(table(all_cazy_in_harboring_genes$cazy_module))
  colnames(cazy_count_in_genes)<-c('CBM_module','Frequency')
  write.table(cazy_count_in_genes,file=paste('cazy_count_in_',endo_cellulase_modules[j],'_genes.tab',sep = ''),sep='\t',row.names = F, col.names = T)
}
####



####
for (j in 1:length(hemicellulase))
{
  all_cazy_in_harboring_genes<-data.frame(Query_ID=c(NA),cazy_module=c(NA),HMM_length=c(NA),Query_length=c(NA),E_value=c(NA),HMM_start=c(NA),HMM_end=c(NA),
                                          Query_start=c(NA),Query_end=c(NA),Coverage=c(NA),desc=c(NA),product=c(NA),locus=c(NA),strand=c(NA),start=c(NA),end=c(NA))
  for (i in 1:length(all_cazy_harboring_nc))
  {
    tmp<-complete_info_nc[[i]]
    tmp$cazy_module<-gsub('.hmm','',tmp$cazy_module)
    harboring_genes<-tmp$Query_ID[(tmp$cazy_module %in% hemicellulase[j])]
    harboring_genes<-unique(harboring_genes)
    all_cazy_in_tmp<-tmp[tmp$Query_ID %in% harboring_genes,]
    all_cazy_in_harboring_genes<-rbind(all_cazy_in_harboring_genes,all_cazy_in_tmp)
  }
  all_cazy_in_harboring_genes<-all_cazy_in_harboring_genes[-1,]
  cazy_count_in_genes<-as.data.frame(table(all_cazy_in_harboring_genes$cazy_module))
  colnames(cazy_count_in_genes)<-c('cazy_module','Frequency')
  write.table(cazy_count_in_genes,file=paste('cazy_count_in_',hemicellulase[j],'_genes.tab',sep = ''),sep='\t',row.names = F, col.names = T)
}
######

for (j in 1:length(LPMO))
{
  all_cazy_in_harboring_genes<-data.frame(Query_ID=c(NA),cazy_module=c(NA),HMM_length=c(NA),Query_length=c(NA),E_value=c(NA),HMM_start=c(NA),HMM_end=c(NA),
                                          Query_start=c(NA),Query_end=c(NA),Coverage=c(NA),desc=c(NA),product=c(NA),locus=c(NA),strand=c(NA),start=c(NA),end=c(NA))
  for (i in 1:length(all_cazy_harboring_nc))
  {
    tmp<-complete_info_nc[[i]]
    tmp$cazy_module<-gsub('.hmm','',tmp$cazy_module)
    harboring_genes<-tmp$Query_ID[(tmp$cazy_module %in% LPMO[j])]
    harboring_genes<-unique(harboring_genes)
    all_cazy_in_tmp<-tmp[tmp$Query_ID %in% harboring_genes,]
    all_cazy_in_harboring_genes<-rbind(all_cazy_in_harboring_genes,all_cazy_in_tmp)
  }
  all_cazy_in_harboring_genes<-all_cazy_in_harboring_genes[-1,]
  cazy_count_in_genes<-as.data.frame(table(all_cazy_in_harboring_genes$cazy_module))
  colnames(cazy_count_in_genes)<-c('cazy_module','Frequency')
  write.table(cazy_count_in_genes,file=paste('cazy_count_in_',LPMO[j],'_genes.tab',sep = ''),sep='\t',row.names = F, col.names = T)
}
######
save.image(file="abundance_of_identified_cazy_in_slh_cellulase_hemicellulase_genes.RData")




