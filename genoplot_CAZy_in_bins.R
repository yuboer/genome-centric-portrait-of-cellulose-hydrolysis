####set the working directory to the folder where files of dbCAN annotation locates
setwd("/envbio/yuboer/yuboer_2nd_binning/cazy_annotation_2nd_binning")
cazy_in_bins_list_0<-dir('/envbio/yuboer/yuboer_2nd_binning/cazy_annotation_2nd_binning')
cazy_in_bins_list<-sapply(strsplit(cazy_in_bins_list_0,'.faa.out.dm.ps'),'[',1)
#cazy_in_bins_list<-gsub(".faa.out.dm.ps","",cazy_in_bins_list_0)
####read all the files and store them in a list, namely 'cazy_in_bins_hmm' here:
cazy_in_bins_hmm<-list()
for (i in 1: length(cazy_in_bins_list))
{
  temp<-read.delim(file = paste("/envbio/yuboer/yuboer_2nd_binning/cazy_annotation_2nd_binning/",cazy_in_bins_list_0[i],sep=""),header = F, sep='\t')
  colnames(temp)<-c("cazy_module","HMM_length","Query_ID","Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  temp$cazy_module<-gsub(".hmm","",temp$cazy_module)
  #temp$cazy_module<-sapply(strsplit(temp$cazy_module,".hmm") '[',1)
  temp$strand<-c(1)
  cazy_in_bins_hmm[[cazy_in_bins_list[i]]]<-temp
  #cazy_in_bins_hmm$cazy_in_bins_list_0[i]<-temp
  #cazy_in_bins_list_0[i]=cazy_in_bins_hmm[[i]]
}

####genoplot on arrangment of all the cazy modules along genes in the bins:
#cazy_in_bins_hmm_plotfrags<-list()
for (j in 1:length(cazy_in_bins_list))
{
  bins_tmp<-cazy_in_bins_hmm[[j]]
  query_id_in_bins_tmp<-bins_tmp[,'Query_ID']
  uni_query_id_in_bins_tmp<-unique(query_id_in_bins_tmp)
  uni_query_id_in_bins_tmp<-as.character(uni_query_id_in_bins_tmp)
  plotfrag<-list()
  for (i in 1:length(uni_query_id_in_bins_tmp))
  {
    tmp<-bins_tmp[bins_tmp$Query_ID %in% uni_query_id_in_bins_tmp[i],]
    tmp1<-tmp[c("cazy_module", "Query_start","Query_end","strand")]
    colnames(tmp1)<-c("name","start","end","strand")
    tmp1$col<-c("Sky Blue 1")
    tmp1$gene_type<-c("blocks")
    plotfrag[[uni_query_id_in_bins_tmp[i]]]<-as.dna_seg(tmp1)
    #plotfrag[[i]]<-as.dna_seg(tmp1)
  }
  #names(plotfrag)<-uni_query_id_in_bins_tmp
  allannots <- lapply(plotfrag, function(x)
  {mid <- middle(x)
  annot <- annotation(x1=mid, text=x$name, rot=40, col='Gray 41')})
  #cazy_in_bins_hmm_plotfrags[[cazy_in_bins_list_0[j]]]<-plotfrag
  #cazy_in_bins_hmm[[j]]<-plotfrag
  #pdf(file="test.pdf",height=100,width=20)
  pdf(file=paste(cazy_in_bins_list[j],".pdf",sep=""),height=130,width=30)
  plot_gene_map(plotfrag,annotations = allannots,annotation_cex =0.8,
                main=paste("Arrangement of the CAZy modules in",cazy_in_bins_list[j],sep=" "),cex.main=1.2,col.main="Gray 41",
                dna_seg_labels =uni_query_id_in_bins_tmp, dna_seg_label_cex = 0.8,dna_seg_label_col="Gray 41")
  #plot_gene_map(cazy_in_bins_hmm_plotfrags[[cazy_in_bins_list_0[j]]],annotations = allannots,annotation_cex =0.8,dna_seg_label_cex=1)
  #plot_gene_map(cazy_in_bins_hmm_plotfrags[j],annotations = allannots,annotation_cex =0.8,dna_seg_label_cex=1)
  dev.off()
}
