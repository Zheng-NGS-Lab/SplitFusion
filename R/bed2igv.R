bed2igv <- function(pos,bed,name,chr,n.reads=10){
    igv <- igvR()

    setBrowserWindowTitle(igv,name)
    setGenome(igv, "hg19")

    Sys.sleep(5)   # wait a few seconds before zooming into MEF2C
    #showGenomicRegion(igv, "EML4")
    showGenomicRegion(igv, list(chrom=chr, start=pos-75, end=pos+75))

    #bam.data <- data.frame(do.call(rbind,mclapply(system(paste0("samtools view ",bam), intern = TRUE),function(x){str_split(x,"\t")[[1]]})))
    #bam.gr <- GRanges(seqnames = as.character(bam.data$X3),
    #                  IRanges(start = ifelse(as.numeric(as.character(bam.data$X4))<as.numeric(as.character(bam.data$X8)),as.numeric(as.character(bam.data$X4)),as.numeric(as.character(bam.data$X8))),
    #                          end = ifelse(as.numeric(as.character(bam.data$X4))<as.numeric(as.character(bam.data$X8)),as.numeric(as.character(bam.data$X8)),as.numeric(as.character(bam.data$X4)))),
    #                  name=as.character(bam.data$X1))

    bam.data <- data.frame(fread(bed))
    bam.data <- subset(bam.data,V3 == chr & V4 > pos-150 & V5 < pos+150)

    #bam.bed <- (bam.data[,c("X3","X4","X8","X1")])

    #bam.data$X1 <- as.character(bam.data$X1)

    #colors <- terrain.colors(10)
    #color.index <- 0

    ### read orientation ###
    bam.data[,"read_ori"] <- gsub(".+umi:(.+)/(\\d)","\\2",bam.data$V1)
    ### read part ###
    bam.data[,"read_part"] <- as.numeric(gsub("(^\\d*) .+umi:(.+/\\d)","\\1",bam.data$V1))

    #### read name ###
    #bam.data$V1 <- gsub(".+umi:(.+)/(\\d)","\\1",bam.data$V1)
    bam.data$V1 <- gsub(".+umi:(.+/\\d)","\\1",bam.data$V1)

    for (i in head(unique(bam.data$V1),n.reads)){
      #color.index <- color.index + 1
      dd <- bam.data[bam.data$V1==i,]
      for (j in 1:nrow(dd)){

        track.read <- DataFrameAnnotationTrack(i, dd[j,c("V3","V4","V5","V1")],
                                             #trackHeight=20, color=ifelse(dd[j,"strand"]=="1","red","green"))#color=colors[color.index])
                                             trackHeight=20, color=ifelse(dd[j,"read_ori"]=="1",
                                                                          colorRampPalette(c("red","black"))(ifelse(nrow(dd) < dd[j,"read_part"],dd[j,"read_part"],nrow(dd)))[dd[j,"read_part"]],
                                                                          colorRampPalette(c("green","black"))(ifelse(nrow(dd) < dd[j,"read_part"],dd[j,"read_part"],nrow(dd)))[dd[j,"read_part"]]))
        displayTrack(igv, track.read)
        Sys.sleep(5)
      }
    }
  }

