##==== 4: sorting direction
##	1) reverse read2 Left/Right to get the same as Read1
SplitFusion.breakpoint.anno.postscript.direction = function(lr2b){
r1 = lr2b[!(grepl("/2", lr2b$readID)),]
r2 = lr2b[grepl("/2", lr2b$readID),]

    if (nrow(r2)>0){
	r2r = r2
	    tmp.n = gsub('_L', '_Tmp000R', names(r2))
	    tmp.n2 = gsub('_R', '_Tmp000L', tmp.n)
	    tmp.n3 = gsub('_Tmp000', '_', tmp.n2)
	    names(r2r) = tmp.n3
	    r2r$strand_L[r2r$strand_L == '+'] = 'm'
	    r2r$strand_L[r2r$strand_L == '-'] = '+'
	    r2r$strand_L[r2r$strand_L == 'm'] = '-'
	    r2r$strand_R[r2r$strand_R == '+'] = 'm'
	    r2r$strand_R[r2r$strand_R == '-'] = '+'
	    r2r$strand_R[r2r$strand_R == 'm'] = '-'
	lr2c = rbind(r1, r2r)
    } else {lr2c = lr2b}

##	2) make Left as 5' and Right as 3'
# tab.sen = ddply(lr2c, .(strand_L, geneStrand_L, strand_R, geneStrand_R), summarize, n=length(readID))
sense = subset(lr2c, (strand_L == geneStrand_L & strand_R == geneStrand_R)
			| (is.na(geneStrand_L) & strand_R == geneStrand_R)
			| (is.na(geneStrand_R) & strand_L == geneStrand_L)
)
antisense = subset(lr2c, (strand_L != geneStrand_L & strand_R != geneStrand_R)
			| (is.na(geneStrand_L) & strand_R != geneStrand_R)
			| (is.na(geneStrand_R) & strand_L != geneStrand_L)
)
nosense = subset(lr2c, is.na(geneStrand_L) & is.na(geneStrand_R)
			| (strand_L == geneStrand_L & strand_R != geneStrand_R)
			| (strand_L != geneStrand_L & strand_R == geneStrand_R)
			)
#nrow(sense); nrow(antisense); nrow(nosense)
## reverse antisend L/R to get the same as sense
if (nrow(antisense)>0){
	anti2 = antisense
		tmp.n = gsub('_L', '_Tmp000R', names(antisense))
        	tmp.n2 = gsub('_R', '_Tmp000L', tmp.n)
        	tmp.n3 = gsub('_Tmp000', '_', tmp.n2)
        names(anti2) = tmp.n3
    lr2d = rbind(sense, anti2, nosense)
} else {lr2d = lr2c}
#nrow(lr2b); nrow(lr2c); nrow(lr2d)

n.lr2d = nrow(lr2d)
return(lr2d)
}
