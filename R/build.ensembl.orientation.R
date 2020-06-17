build.ensembl.orientation <- function(attributes.name=c('ensembl_gene_id','ensembl_transcript_id','refseq_mrna','hgnc_symbol','strand')){
  require(biomaRt)
  mart=useMart("ensembl")
  mart=useDataset("hsapiens_gene_ensembl", mart = mart)
  ensembl.orientation <- getBM(attributes=attributes.name, mart = mart)
  ensembl.orientation$strand <- ifelse(ensembl.orientation$strand==1,"+","-")
  write.csv(ensembl.orientation[!duplicated(ensembl.orientation$ensembl_transcript_id),],"ENSEMBL.orientation.txt",row.names = FALSE,quote = FALSE,na = "")
}
