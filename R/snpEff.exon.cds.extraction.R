
#' exon.cds.extraction
#'
#' Extracting exon or cds informtaion of primary annotation from snpEff.
#'
#' @param input annotation file outputted by snpEff.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' exon.cds.extraction(input = "mid.anno")
#'
#' exon.cds.extraction(input = "__breakpoint.annotated")
#'
snpEff.exon.cds.extraction <- function(input){
library(data.table)
library(parallel)
options(width=204)

in.anno = data.frame(fread(input)) ###Modified by Baifeng###
#in.anno.exon = commandArgs(TRUE)[2]
#configFile = commandArgs(TRUE)[2] ###Modified by Baifeng###
#source(configFile)
#cds <- data.frame(fread(in.anno.cds))
#exon <- data.frame(fread(in.anno.exon))
anno <- mclapply(1:nrow(in.anno),function(x) {
	y = strsplit(in.anno[x,8],"\\|")
	#	c(paste0(in.anno[x,1],"_",in.anno[x,2]),y[[1]][4],"NA",y[[1]][8],y[[1]][2],y[[1]][7],paste0(ifelse(grepl("p",y[[1]][11]),"exon","intron"),sub("/\\d+","",y[[1]][9])),sub("[c|n].*?(\\d+).*","\\1",y[[1]][10]))
	c(paste0(in.anno[x,1],"_",in.anno[x,2]),ifelse(y[[1]][4]=="","NA",y[[1]][4]),"NA",ifelse(y[[1]][8]=="",NA,y[[1]][8]),ifelse(y[[1]][2]=="",NA,y[[1]][2]),ifelse(y[[1]][7]=="",NA,y[[1]][7]),paste0(ifelse(y[[1]][11]!="","exon","intron"),ifelse(y[[1]][9]!="",sub("/\\d+","",y[[1]][9]),"NA")),ifelse(y[[1]][13]!="",sub("/\\d+","",y[[1]][13]),"NA"))
	})

write.table(data.frame(do.call(rbind,anno)),paste0(input,".ext0"),col.names=F, row.names=F, quote=F, sep='\t')
}

