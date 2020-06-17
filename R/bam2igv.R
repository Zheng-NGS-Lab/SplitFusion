#setwd("/Users/baizha/Desktop/KIproject/Project/Fusion/test/SplitFusion/inst/data/example_data/result/example/")


#fq2igv <- function(config){
  
  #system(paste0(bwa," mem -T 18 -t ",thread," ",database_dir,"/",refGenome," ",fastq_dir,"/",r1filename," |",samtools," view -@ ",thread," -bS - > ",paste0(fastq_dir,"/",r1filename,".bam")))
  
  bam2igv <- function(bamfile,outname = bamfile, win.size=100){
    
    #### New version ###
  
    require("igvR")
    #outname <- "L18-00527T.MET_exon13---MET_exon15"
    igv <- igvR()
    #setBrowserWindowTitle(igv,outname)
    setGenome(igv, "hg19")
    #Sys.sleep(5)   # wait a few seconds before zooming into MEF2C
  
  
    x <- readGAlignments(bamfile, use.names=TRUE)
    track <- GenomicAlignmentTrack(outname, x, trackHeight = 300)
    displayTrack(igv, track)
  
    for (i in unique(seqnames(x)@values)) {
      x.p <- c(start(x[which(strand(x)=="-" & seqnames(x)==i)]),end(x[which(strand(x)=="-" & seqnames(x)==i)]),
               start(x[which(strand(x)=="+" & seqnames(x)==i)]),end(x[which(strand(x)=="+" & seqnames(x)==i)]))
      x.p.start <- min(x.p)
      x.p.end <- max(x.p)
    
      if (x.p.end - x.p.start < win.size){
      
        showGenomicRegion(igv, list(chrom=i, start=round((x.p.start+x.p.end)/2)-round(win.size/2), end=round((x.p.start+x.p.end)/2)+round(win.size/2)))
        Sys.sleep(5)
        saveToSVG(obj = igv,filename = paste0(outname,".",i,".",x.p.start,".",x.p.end,".svg"))
      
      }else{
      
        showGenomicRegion(igv, list(chrom=i, start=x.p.start-round(win.size/2), end=x.p.start+round(win.size/2)))
        Sys.sleep(5)
        saveToSVG(obj = igv,filename = paste0(outname,".",i,".",x.p.start,".svg"))
    
        showGenomicRegion(igv, list(chrom=i, start=x.p.end-round(win.size/2), end=x.p.end+round(win.size/2)))
        Sys.sleep(5)
        saveToSVG(obj = igv,filename = paste0(outname,".",i,".",x.p.end,".svg"))
      
      }
    }
  }
  
#  igv.bed.show(pos = 42492091,name = "EML4_intron6",chr = 2,bed = "breakpoint.reads")
#library(parallel)
#mclapply(list.files("./"),function(x){system(paste0(bwa," mem -T 18 -t ",thread," ",database_dir,"/",refGenome," ",fastq_dir,"/",x," |",samtools," view -@ ",thread," -bS - > ",paste0(fastq_dir,"/",x,".bam")))})
#mclapply(list.files("./",pattern = ".bam"),function(x){igv.bam.show(bamfile = x)})

#setwd("/Users/baizha/Dropbox/SplitFusion/Picture/SplitFusionManuscript/Fig4/fusion.sequences")
#for (i in list.files("./",pattern = "Lib*.*.bam$")){
#  igv.bam.show(bamfile = i)
#  #Sys.sleep(5)
#  #print(i)
#}

