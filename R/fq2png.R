#' fq2png
#'
#' Converted fastq sequnce into png file.
#'
#' @param seq raw fastq suquence of breakpoint
#' @param aln Transformed cigar records of breakpoint reads
#' @param nam prefixed name of output file
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' fq2png(seq = "example.EML4_intron6---ALK_exon20.txt", aln = "breakpoint.reads"
#' ,nam = "example.EML4_intron6---ALK_exon20")
#'
fq2png <- function(seq,aln,nam){

  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library(data.table)
  fa <- data.frame(fread(seq,col.names = c("RD", "flag", "seq")))

### remove secondary alignment ###
#fa <- fa %>% dplyr::filter(flag <1000)
  fa <- subset(fa,flag < 1000)

### look for read 1 and read 2? ###
  fa$read <- sapply(fa$flag, function(x) ifelse(intToBits(x)[8]==01, 2, 1))
  fa$len <- sapply(fa$seq, function(x) nchar(x))
  fa$RD <- paste(fa$RD,fa$read,sep = "/")
  fa <- fa[fa$read==2,]
### too slow, use grep in sh ###
### read in partition ###
#L<-readLines("/home/a/Documents/Split_fusion_all/visualisation/SR_info")
#L<-gsub(" ", "\t", L)
#L<-gsub("/", "\t", L)
#bk<-read.table(text=L, sep="\t", as.is = T)

  bk <- read.csv(aln,sep = "",header = F)
### mk alignment table ###
  colnames(bk)[1:3]<-c("part_order", "RD", "len")
  bk$RD <- as.character(bk$RD)


  BK_table <- right_join(bk,fa, by=c("RD"))
### plot... table? ###
  BK_table <- BK_table %>% dplyr::mutate(part_seq=substring(seq, V10, V11)) %>% select(RD, part_order, seq, V10, V11)

### convert it to a table
  seq_order=c()
  RD_col=c()
  part_order=c()
  base=c()
  location=c()
  for (i in 1:nrow(BK_table)){
    for (j in BK_table[i, "V10"]:BK_table[i, "V11"]){
      seq_order<-c(seq_order, i)
      base<-c(base, substring(BK_table[i, "seq"], j, j))
      RD_col<-c(RD_col, BK_table[i, "RD"])
      location<-c(location, j)
      part_order<-c(part_order, BK_table[i, "part_order"])
    }
  }

  base_table<-data.frame(seq_order, RD=RD_col, part_order, base, location, stringsAsFactors = F)

  title_name <- nam
  library(ggsci)
#g<-base_table %>% ggplot(aes(y=RD, x=location, label=base, fill=part_order, ymin=RD, xmin=location, ymax=RD, xmax=location)) +
#  geom_text(aes(y=RD, x=location, label=base), colour="black")+
#  geom_rect(aes(ymin=0, ymax=RD, xmin=location, fill=part_order, xmax=location+0.9), alpha=0.2)+
#  facet_grid(RD ~., scale="free_y")+
#  ylab("read coordinate")+
#  ggtitle(title)+
#  scale_fill_npg()+
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank(),
#        legend.position = "none")
#ggsave(g, file="/home/a/Documents/Split_fusion_all/visualisation/test_plot.pdf")

  base_table$part_order <- factor(base_table$part_order,levels = sort(unique(base_table$part_order)))

  g <- ggplot(base_table)+geom_text(aes(y=RD,x=location,label=base,color=part_order))+
  #geom_point(aes(y=RD,x=location,fill=part_order),size=3,shape=21,alpha=0.5)+
    ggtitle(title_name)+
    theme(axis.title = element_text(face = "bold"),axis.text = element_text(face = "bold"),panel.grid = element_blank(),panel.background = element_blank())
  ggsave(g, file=paste0(title_name,".png"),height = 15, width = 15)
}



