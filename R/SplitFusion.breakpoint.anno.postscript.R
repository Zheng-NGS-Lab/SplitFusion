# breakpoint post annotation processing
# 
# required: 
#	- fu.anno.right
#	- fu.anno.left
#	- fu.breakpointID.stgLeft
#	- fu.breakpointID.stgRight
#
SplitFusion.breakpoint.anno.postscript = function(){

library(data.table)
library(plyr)
options(width=204)
options(scipen=999)
(SampleId = sub('.*/', '', getwd()))
source('config.txt')
source(paste0(SFpath,"/R/SplitFusion.breakpoint.anno.postscript.mid.anno.R"))
source(paste0(SFpath,"/R/SplitFusion.breakpoint.anno.postscript.direction.R"))
source(paste0(SFpath,"/R/SplitFusion.breakpoint.anno.postscript.direction.sub.R"))

##====1: left-right annotations
    colnames = c('readID', 'chrorp_L', 'chr_L', 'orp_L', 'pos_L', 'strand_L'
		, 'num_unique_molecules', 'num_partner_ends', 'num_partner_ends2', 'breakpoint', 'overlap'
		, 'gene_L', 'geneStrand_L', 'inEx_L', 'functiontype_L', 'nm_L', 'exon_L', 'cdna_L'
		, 'chrorp_R', 'chr_R', 'orp_R', 'pos_R', 'strand_R', 'gene_R', 'geneStrand_R', 'inEx_R', 'functiontype_R', 'nm_R', 'exon_R', 'cdna_R')
    lr1 = fread('anno.left.right', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnames)
    lr1$cdna_L = suppressWarnings(as.numeric(lr1$cdna_L))
    lr1$cdna_R = suppressWarnings(as.numeric(lr1$cdna_R))
    lr1$pos_L = suppressWarnings(as.numeric(lr1$pos_L))
    lr1$pos_R = suppressWarnings(as.numeric(lr1$pos_R))
    #head(lr1)
    n.lr1 = nrow(lr1)
    lr1$exonn_L = suppressWarnings(as.numeric(sub('exon','',lr1$exon_L)))
    lr1$exonn_R = suppressWarnings(as.numeric(sub('exon','',lr1$exon_R)))

    ## kepp original cdna pos for later frameness calculation
    lr1$cdna_L0 = lr1$cdna_L
    lr1$cdna_R0 = lr1$cdna_R

##==== 2: connect left-mid-right
##==== For breakpoints with middle split: to correct breakpoint exon number, cdna, gdna by add/minus middle split size
lr2 = SplitFusion.breakpoint.anno.postscript.mid.anno(lr1)

##==== 3: sorting direction
lr3 = SplitFusion.breakpoint.anno.postscript.direction(lr2)
SplitFusion.breakpoint.anno.postscript.direction.sub(lr3)
}
