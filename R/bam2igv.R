  
  bam2igv <- function(bamfile,outname = bamfile, genomeVer="hg38", win.size=100,chr=NULL,breakpoint=NULL,wait.time=10){
    
    #### New version ###
  
    require("igvR")
    #outname <- "L18-00527T.MET_exon13---MET_exon15"
    saveToSVG_source <-  function(obj, filename){
      send(obj, list(cmd="getSVG", callback="handleResponse", status="request", payload=""))
      while (!browserResponseReady(obj)){
        service(100)
      }
      svgText <- getBrowserResponse(obj);
      file <- file(filename)
      write(svgText, file)
      close(file)
      return(sprintf("%d characters written to %s", nchar(svgText), filename))
    }
      
    igv <- igvR()
    #setBrowserWindowTitle(igv,outname)
    setGenome(igv, genomeVer) 
    #Sys.sleep(5)   # wait a few seconds before zooming into MEF2C
  
  
    x <- readGAlignments(bamfile, use.names=TRUE)
    track <- GenomicAlignmentTrack(outname, x, trackHeight = 300)
    displayTrack(igv, track)
  
    if (is.null(chr) | is.null(breakpoint)){
      for (i in as.character(unique(seqnames(x)@values))) {
        x.p <- c(start(x[strand(x) == "-" & seqnames(x) == i]),end(x[strand(x)=="-" & seqnames(x)==i]),
                 start(x[strand(x) == "+" & seqnames(x) == i]),end(x[strand(x)=="+" & seqnames(x)==i]))
        x.p.start <- min(x.p)
        x.p.end <- max(x.p)
    
        if (x.p.end - x.p.start < win.size){
      
          showGenomicRegion(igv, list(chrom=i, start=round((x.p.start+x.p.end)/2)-round(win.size/2), end=round((x.p.start+x.p.end)/2)+round(win.size/2)))
          Sys.sleep(wait.time)
          saveToSVG_source(obj = igv,filename = paste0(outname,".",i,".",x.p.start,".",x.p.end,".svg"))
        
        }else{
      
          showGenomicRegion(igv, list(chrom=i, start=x.p.start-round(win.size/2), end=x.p.start+round(win.size/2)))
          Sys.sleep(wait.time)
          saveToSVG_source(obj = igv,filename = paste0(outname,".",i,".",x.p.start,".svg"))
        
          showGenomicRegion(igv, list(chrom=i, start=x.p.end-round(win.size/2), end=x.p.end+round(win.size/2)))
          Sys.sleep(wait.time)
          saveToSVG_source(obj = igv,filename = paste0(outname,".",i,".",x.p.end,".svg"))
        
        }
      }
    }else{
      
      showGenomicRegion(igv, list(chrom=chr, start=breakpoint-round(win.size/2), end=breakpoint+round(win.size/2)))
      Sys.sleep(wait.time)
      saveToSVG_source(obj = igv,filename = paste0(outname,".",chr,".",breakpoint,".svg"))
      
    }
    
  }

  
