##==== 3: connect left-mid-right
##==== For breakpoints with middle split: to correct breakpoint exon number, cdna, gdna by add/minus middle split size
SplitFusion.breakpoint.anno.postscript.mid.anno = function(lr2){
if (file.exists('anno.mid')){
    colnamesM = c('annoPos', 'readID', 'chr_M', 'start_M', 'end_M', 'gene_M', 'geneStrand_M', 'inEx_M', 'functiontype_M', 'nm_M', 'exon_M', 'cdna_M')
    #mid = read.table('anno.mid', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnamesM)
    mid = fread('anno.mid', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnamesM)
	if (nrow(mid)>0){
	mid$exonn_M = suppressWarnings(as.numeric(sub('exon', '', mid$exon_M)))
	head(mid)

	lmr = merge(lr2, mid, by="readID")
	if (nrow(lmr) >0){
		lmr_L = subset(lmr, nm_L==nm_M)
		lmr_R = subset(lmr, nm_R==nm_M & nm_L != nm_M)

		# mid belong to Left
		if (nrow(lmr_L)>0){
		    	lmr_L$exon_L = lmr_L$exon_M
		    	lmr_L$exonn_L = lmr_L$exonn_M
			lmr_L$cdna_L = lmr_L$cdna_L + (1 - (lmr_L$exonn_L > lmr_L$exonn_M)*2) * abs(lmr_L$overlap)
			
			lmr_L$absMstart = abs(lmr_L$start_M - lmr_L$orp_L)
			lmr_L$absMend = abs(lmr_L$end_M - lmr_L$orp_L)
			lmr_L$pos_M = ifelse(
				lmr_L$absMstart > lmr_L$absMend
				, lmr_L$start_M
				, lmr_L$end_M
				)
			lmr_L$chrpos_M = paste(lmr_L$chr_L, lmr_L$pos_M, sep='_')
			lmr_L$chrpos_Partner = paste(lmr_L$chr_R, lmr_L$pos_R, sep='_')
			lmr_L$breakpoint = ifelse(
				lmr_L$chrpos_M < lmr_L$chrpos_Partner
				, paste(lmr_L$chrpos_M, lmr_L$chrpos_Partner, sep='__')
				, paste(lmr_L$chrpos_Partner, lmr_L$chrpos_M, sep='__')
				)
			}

		# mid belong to Right
		if (nrow(lmr_R)>0){
		    	lmr_R$exon_R = lmr_R$exon_M
		    	lmr_R$exonn_R = lmr_R$exonn_M
			lmr_R$cdna_R = lmr_R$cdna_R + (1 - (lmr_R$exonn_R > lmr_R$exonn_M)*2) * abs(lmr_R$overlap)
			
			lmr_R$absMstart = abs(lmr_R$start_M - lmr_R$orp_R)
			lmr_R$absMend = abs(lmr_R$end_M - lmr_R$orp_R)
			lmr_R$pos_M = ifelse(
				lmr_R$absMstart > lmr_R$absMend
				, lmr_R$start_M
				, lmr_R$end_M
				)
			lmr_R$chrpos_M = paste(lmr_R$chr_R, lmr_R$pos_M, sep='_')
			lmr_R$chrpos_Partner = paste(lmr_R$chr_L, lmr_R$pos_L, sep='_')
			lmr_R$breakpoint = ifelse(
				lmr_R$chrpos_M < lmr_R$chrpos_Partner
				, paste(lmr_R$chrpos_M, lmr_R$chrpos_Partner, sep='__')
				, paste(lmr_R$chrpos_Partner, lmr_R$chrpos_M, sep='__')
				)
			}
		   
		lmr2 = rbind(lmr_L, lmr_R, fill=TRUE)
		#lmr2.keep = lmr2[, c(names(lr2))]
		lmr2.keep = lmr2[, c(names(lr2)), with=FALSE]
		    lr0 = lr2[!(lr2$readID %in% lmr2.keep$readID),]
		    nrow(lr0); nrow(lmr2.keep); nrow(lr2)
		lr2.update = rbind(lmr2.keep, lr0)
		} else {lr2.update = lr2}
		} else {lr2.update = lr2}
} else {lr2.update = lr2}

	##==== update num_start_site by bk
	lr2.update$stst = sub("[^0-9].*", "", sub(".*P", "", lr2.update$readID))
	nssbk = ddply(lr2.update, .(breakpoint), summarize, 'nss.update' = length(unique(stst)), 'umi.update' = length(unique(readID)))
	lr2b = merge(lr2.update, nssbk, by='breakpoint')
	lr2b$num_start_site = lr2b$nss.update
	lr2b$num_unique_molecules = lr2b$umi.update

return(lr2b)
}
