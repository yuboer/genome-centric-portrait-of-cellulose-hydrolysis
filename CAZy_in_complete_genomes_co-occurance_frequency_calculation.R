#this script was written to calculate the CAZy co-occuring frequencies along same genes in the 2786 complete genomes investigated
setwd("/workhome/yuboer_CAZy_in_2786complete_genome")
investigated_modules_file<-dir(path = "/workhome/yuboer_CAZy_in_2786complete_genome/co-occurance_gh_slh_cohe/cellulase_hemicellulase_cohesin_slh_dok/")

inves_gh_slh_coh_dok<-data.frame(cazy_module=c(NA))
for (i in 1:length(investigated_modules_file))
{
  tmp<-read.delim(file=paste("/workhome/yuboer_CAZy_in_2786complete_genome/co-occurance_gh_slh_cohe/cellulase_hemicellulase_cohesin_slh_dok/",investigated_modules_file[i],sep = ""),header = F,sep='\t' )
  colnames(tmp)<-c('cazy_module',strsplit(investigated_modules_file[i],"\\.")[[1]] [1])
  inves_gh_slh_coh_dok<-merge(inves_gh_slh_coh_dok,tmp,by.x ='cazy_module',by.y="cazy_module",all=T)
}

write.table(inves_gh_slh_coh_dok,file='/workhome/yuboer_CAZy_in_2786complete_genome/co-occurance_gh_slh_coh_dok.tab',sep='\t',row.names = F,col.names = T)

exo_cellulase_modules<-c('GH6','GH9','GH48')

endo_cellulase_modules<-c('GH5','GH7','GH8','GH12','GH44','GH45','GH51','GH61','GH74','GH124')
#not_identified_endo_cellulase<-c('GH7','GH61')

hemicellulase<-c('GH10','GH11','GH16','GH26','GH28','GH53','GH113','GH134')
#not_identified_hemicellulase_module<-c('GH134')
inves_gh_slh_coh_dok_modules<-c("cazy_module","slh",'cohesin','dockerin',exo_cellulase_modules,endo_cellulase_modules,hemicellulase)

inves_gh_slh_coh_dok<-inves_gh_slh_coh_dok[inves_gh_slh_coh_dok_modules[inves_gh_slh_coh_dok_modules %in% colnames(inves_gh_slh_coh_dok)]]

inves_gh_slh_coh_dok_t<-t(inves_gh_slh_coh_dok)

cm<-inves_gh_slh_coh_dok_t["cazy_module",]
colnames(inves_gh_slh_coh_dok_t)<-cm
inves_gh_slh_coh_dok_t<-inves_gh_slh_coh_dok_t[-1,]

colnames(inves_gh_slh_coh_dok_t)
inves_gh_slh_coh_dok_t<-inves_gh_slh_coh_dok_t[,!is.na(colnames(inves_gh_slh_coh_dok_t))]
celllulose_binding_CBM<-c("CBM1","CBM2","CBM3","CBM4","CBM6","CBM7","CBM8","CBM9","CBM10","CBM11","CBM16","CBM17","CBM28","CBM30","CBM37","CBM44","CBM46","CBM49","CBM63","CBM64")
hemicellulose_binding_CBM<-c('CBM13','CBM15','CBM22','CBM23','CBM24','CBM27','CBM29','CBM31','CBM35','CBM39','CBM41','CBM43','CBM52','CBM54','CBM56','CBM59','CBM60','CBM61','CBM62','CBM65','CBM66')
sel_cbm<-c(celllulose_binding_CBM,hemicellulose_binding_CBM)
co_cbm<-colnames(inves_gh_slh_coh_dok_t)[grep("CBM",colnames(inves_gh_slh_coh_dok_t))]

co_sel_cbm_in_seq<-sel_cbm[sel_cbm %in% co_cbm]
co_cbm_others<-co_cbm[!(co_cbm %in% co_sel_cbm_in_seq)]
co_cbm_inseq<-c(co_sel_cbm_in_seq,co_cbm_others)

co_others<-colnames(inves_gh_slh_coh_dok_t)[!(colnames(inves_gh_slh_coh_dok_t) %in% co_cbm )]
co_modules_inseq<-c(co_cbm_inseq,co_others)

inves_gh_slh_coh_dok_co_occur<-inves_gh_slh_coh_dok_t[,co_modules_inseq[co_modules_inseq %in% colnames(inves_gh_slh_coh_dok_t)]]
write.table(inves_gh_slh_coh_dok_co_occur,file="co-occurance_slh_coh_dok_gh.tab",row.names = T,col.names = T,sep='\t')

save.image(file="co-occurance_merge.RData")
