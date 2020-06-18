#=====================================================================
#===== Final analysis: frame status, filters, exon-junction, output
#=====================================================================

SplitFusion.breakpoint.anno.postscript.direction.sub = function(lr3){
n.lr3 = nrow(lr3)
sampleID = sub(".*\\/", "", getwd())

if (n.lr3 >0){
	## add the 6 base back to cDNA position
	lr3$cdna_L = lr3$cdna_L0 + 6
	lr3$cdna_R = lr3$cdna_R0 - 6

	## intra-gene
	lr3$intragene = ifelse(lr3$gene_L == lr3$gene_R, 1, 0)
	lr3$exon = substr(lr3$exon_L,1,4)
	lr3$exonD = 0
	lr3$exonD = ifelse(lr3$gene_L == lr3$gene_R, lr3$exonn_R - lr3$exonn_L, lr3$exonD)
	lr3$gec_L = paste0(lr3$gene_L, ' ', lr3$exon_L, ' c.', lr3$cdna_L, ' (', lr3$nm_L, ')')
	lr3$gec_R = paste0(lr3$gene_R, ' ', lr3$exon_R, ' c.', lr3$cdna_R, ' (', lr3$nm_R, ')')
	lr3$gec_LR = paste(lr3$gec_L, lr3$gec_R, sep='---')
	lr3$gg = paste0(lr3$gene_L, '_', lr3$gene_R)
	lr3$ge1 = paste0(lr3$gene_L, '_', lr3$exon_L)
	lr3$ge2 = paste0(lr3$gene_R, '_', lr3$exon_R)
	lr3$ge1ge2 = paste0(lr3$ge1,'---', lr3$ge2)

	##======================
	## frame status
	##======================
	    lr3$frameD = (lr3$cdna_R0 %% 3 - (lr3$cdna_L0 - lr3$overlap) %% 3) %% 3
	    lr3$frame = ifelse(lr3$frameD ==1, 'in-frame', 'out-frame')
	    lr3$frame[is.na(lr3$frame)] = 'N.A.'
		#head(lr3)
		## expect mostly in-frame or _NA_. Otherwise, might need inspection
		#table(lr3$frame,useNA='always')
				     
	    lr3$neighb = ifelse(lr3$exonD==1, 1, 0)
	    lr3$neighb[is.na(lr3$neighb)] = 0

		## order
		lr3 = lr3[order(lr3$intragene, lr3$frame, lr3$gec_LR),]

	##===========================================================
	##=== Optional: Backend database supported target output or filtering
	##===========================================================
	    ##==== Targeted output known significant fusions or splicing isoforms (e.g. MET ex14 skipping)
		known.gg=NA; known.ge=NA; known.partner=NA; known.3UTR=NA;
	    known.gg.file = paste0(panel_dir, '/known.gene-gene.txt'); if (file.exists(known.gg.file)){known.gg = readLines(known.gg.file)}
	    known.ge.file = paste0(panel_dir, '/known.gene-exon.txt'); if (file.exists(known.ge.file)){known.ge = readLines(known.ge.file)}
	    known.partner.file = paste0(panel_dir, '/known.partners.txt'); if (file.exists(known.partner.file)){known.partner = readLines(known.partner.file)}
	    known.3UTR.truncation.file = paste0(panel_dir, '/known.3UTR.truncation.txt'); if (file.exists(known.3UTR.truncation.file)){known.3UTR = readLines(known.3UTR.truncation.file)} 

	    ##==== Targeted remove unwanted fusions (e.g. involving homologous genes or recurrent falsed positives)
		filter.gg=NA; filter.ge=NA; 
	    filter.gg.file = paste0(panel_dir, '/filter.gene-gene.txt'); if (file.exists(filter.gg.file)){filter.gg0 = read.table(filter.gg.file, header=F, stringsAsFactors=F);
				 filter.gg = c(paste(filter.gg0$V1, filter.gg0$V2, sep='_'), paste(filter.gg0$V2, filter.gg0$V1, sep='_'))} 
	    filter.ge.file = paste0(panel_dir, '/filter.gene-exon.txt'); if (file.exists(filter.ge.file)){filter.ge = readLines(filter.ge.file)}

	    lr3$known=0
	    lr3$known[(lr3$gene_T == lr3$gene_L & lr3$gene_R %in% known.partner & lr3$intragene==0)] = 1
	    lr3$known[(lr3$gene_T == lr3$gene_R & lr3$gene_L %in% known.partner & lr3$intragene==0)] = 1
	    lr3$known[(lr3$ge1ge2 %in% known.ge)] = 1
	    lr3$known[(lr3$ge1 %in% known.3UTR & lr3$intragene==0)] = 1
	    lr32 = subset(lr3, (known==1 | (num_partner_ends >= as.numeric(FusionMinStartSite) & intragene==0)) & !(gg %in% filter.gg | ge1 %in% filter.ge | ge2 %in% filter.ge))
		    # nrow(lr3); nrow(lr32)
		    # output for furture research
		    write.table(lr32, paste(sampleID, '.fusion.list.pre-processing.txt', sep=''), row.names=F, quote=F, sep='\t')

	##======================
	## Further filters
	##======================
	## To remove:
	#	- read-through (defined as within 100K distance on chrosome)
	#	- homologous genes
	#	- neighbour exons of same gene
	#	- out-frame (keep out-frame of different genes. FGFRs fusion can be functional even out-frame)
	#	- same exon
	#
	## To calculate:
	#	- exon boundary/junction
	if (nrow(lr32)>0){
		lr32$read.through = 0
		lr32$read.through[lr32$gene_L != lr32$gene_R 
					& lr32$chr_L == lr32$chr_R 
					& abs(lr32$pos_L - lr32$pos_R)<100000
					]=TRUE

		lr32$gene_Lchr = gsub('[0-9]', '', lr32$gene_L)
		lr32$gene_Rchr = gsub('[0-9]', '', lr32$gene_R)
		homolog.match = function(row){
				h1 = pmatch(row[1], row[2])
				h2 = pmatch(row[2], row[1])
				h3 = (h1 || h2)
              			return(h3)
			}
		lr32$homolog = apply(lr32[,c('gene_Lchr', 'gene_Rchr')], 1, homolog.match)
		lr32$homolog[is.na(lr32$homolog)] = 0
		lr32$homolog[lr32$known==1] = 0

		lr40 = subset(lr32, !(read.through | homolog | neighb==1 | (intragene==1 & exon_L==exon_R)
					| (intragene==1 & frame=='out-frame')
					))
		lr4 = lr40[!duplicated(lr40$readID),]
		n.lr40 = nrow(lr40)
		n.lr4 = nrow(lr4)
		
		if (n.lr4>0){
			lr4$SampleID = sampleID
			# remove illegal symbol
			lr4$gec_LR = gsub('\\(', '.', lr4$gec_LR)
			lr4$gec_LR = gsub('\\)', '.', lr4$gec_LR)
			lr4$ge1ge2 = gsub('\\(', '.', lr4$ge1ge2)
			lr4$ge1ge2 = gsub('\\)', '.', lr4$ge1ge2)

			##======================
			## exon junction support based on RefGene exon starts/ends
			##======================
				bk = lr4[, c('breakpoint', 'gene_L', 'pos_L', 'gene_R', 'pos_R')]
					# most representative pos
				bk2 = ddply(bk, .(breakpoint, gene_L, gene_R), summarize
					, 'pos_L' = names(sort(table(pos_L), decreasing = TRUE)[1])
					, 'pos_R' = names(sort(table(pos_R), decreasing = TRUE)[1]))
					
					bk.genes = unique(c(bk2$gene_L, bk2$gene_R))
				refGene = read.table(paste0(path.package("SplitFusion"), '/data/refGene0.txt'), header=T, sep='\t')
				refGene2 = subset(refGene, name2 %in% bk.genes)[,c('name2', 'exonStarts', 'exonEnds')]

					toStarts = function(row){
						dist = as.numeric(unlist(strsplit(as.character(row[1]),','))) - as.numeric(row[2])
						if (any(dist %in% 0:2)) {junction=1} else {junction=0}
						return(junction)
						}

					fromEnds = function(row){
						dist = as.numeric(unlist(strsplit(as.character(row[1]),','))) - as.numeric(row[2])
						if (any(dist %in% -2:0)) {junction=1} else {junction=0}
						return(junction)
						}
					
				junc.L = unique(merge(bk2[, c('breakpoint', 'gene_L', 'pos_L')], refGene2, by.x='gene_L', by.y='name2', all.x=T))
					junc.L$distLstart = apply(junc.L[,c('exonStarts', 'pos_L')], 1, toStarts)
					junc.L$distLend = apply(junc.L[,c('exonEnds', 'pos_L')], 1, fromEnds)

				junc.R = unique(merge(bk2[, c('breakpoint', 'gene_R', 'pos_R')], refGene2, by.x='gene_R', by.y='name2', all.x=T))
					junc.R$distRstart = apply(junc.R[,c('exonStarts', 'pos_R')], 1, toStarts)
					junc.R$distRend = apply(junc.R[,c('exonEnds', 'pos_R')], 1, fromEnds)

				junc.LR = merge(junc.L, junc.R, by='breakpoint')[,c('breakpoint', 'distLstart', 'distLend', 'distRstart', 'distRend')]
				junc.LR$exon.junction = junc.LR$distLstart + junc.LR$distLend + junc.LR$distRstart + junc.LR$distRend

				#lr5 = merge(lr4, unique(junc.LR[,c('breakpoint', 'exon.junction')]), by='breakpoint')
				lr5 = merge(lr4, unique(junc.LR[,c('breakpoint', 'exon.junction')]), by='breakpoint', allow.cartesian=TRUE)

		## output
			lr5 = lr5[order(-lr5$exon.junction, -lr5$known, -lr5$num_partner_ends, lr5$intragene, lr5$frame),]
			fusion.list = subset(lr5, select=c('SampleID', 'ge1ge2', 'frame','gec_LR','chr_L','pos_L','inEx_L','gene_L','nm_L','exon_L','cdna_L'
						, 'overlap','chr_R','pos_R','inEx_R','gene_R','nm_R','exon_R','cdna_R','intragene','readID'
						, 'breakpoint', 'num_partner_ends', 'known', 'exon.junction'))
				write.table(fusion.list, paste(sampleID, '.fusion.list.post-processing.txt', sep=''), row.names=F, quote=F, sep='\t')

			## consolidate breakpoint
				fusion.list[is.na(fusion.list)]='-'
				fusion.list$readID2 = sub(':umi:.*', '', fusion.list$readID)
			fusion.tab1 = ddply(fusion.list, .(breakpoint, SampleID, ge1ge2, frame, nm_L, inEx_L, gene_L, cdna_L, nm_R, inEx_R, gene_R, cdna_R, intragene, overlap, num_partner_ends, known), summarize
							, 'num_unique_reads'=length(unique(readID2)), 'exon.junction' = max(exon.junction)
							)
				fusion.tab1$transcript_L = fusion.tab1$nm_L
				fusion.tab1$transcript_R = fusion.tab1$nm_R
				fusion.tab1$function_L = fusion.tab1$inEx_L
				fusion.tab1$function_R = fusion.tab1$inEx_R
			fusion.tab1 = fusion.tab1[order(-fusion.tab1$exon.junction, -fusion.tab1$known, -fusion.tab1$num_partner_ends, fusion.tab1$intragene, fusion.tab1$frame, -fusion.tab1$num_unique_reads),]

			## keep first of same breakpoint (exon.junction, most abundant nss, reads)
			fusion.tab2 = fusion.tab1[!duplicated(fusion.tab1$breakpoint),]
				head(fusion.tab2)

		out.names = c("SampleID", "ge1ge2", "frame", "num_partner_ends", "num_unique_reads", "exon.junction", "breakpoint"
				, "transcript_L", "transcript_R", "function_L", "function_R", "gene_L", "cdna_L", "gene_R", "cdna_R", "intragene", "known")
		fusion.table = fusion.tab2[, out.names]
		fusion.table = fusion.table[order(-fusion.tab2$exon.junction, -fusion.table$known, fusion.table$frame, fusion.table$intragene, -fusion.table$num_partner_ends),]
				name1 = names(fusion.table)
				name2 = sub('_L', "_5", name1)
				name3 = sub('_R', "_3", name2)
				name4 = sub('ge1ge2', "GeneExon5_GeneExon3", name3)
				names(fusion.table) = name4
				fusion.table$exon.junction[fusion.table$exon.junction==2] = 'Both'
				fusion.table$exon.junction[fusion.table$exon.junction==1] = 'One'
			write.table(fusion.table, paste(sampleID, '.fusion.table.NoFilter.txt', sep=''), row.names=F, quote=F, sep='\t')

	##=======================================
	## Filter: exon boundary
	##=======================================
		    # inter-gene, in-frame
		    # or known
			fusion.table$g5g3 = paste(fusion.table$gene_5, fusion.table$gene_3, sep='_')
			minPartnerEnds_BothExonJunction = as.numeric(minPartnerEnds_BothExonJunction)
			minPartnerEnds_OneExonJunction = as.numeric(minPartnerEnds_OneExonJunction)
				
			ex = subset(fusion.table, (known | (g5g3 %in% known.gg & (exon.junction != '0') | frame =='gDNA') ## for future compatibility with gDNA reads
							 | (exon.junction == 'Both' 
							    & (known ==1 | (frame =='in-frame' & num_partner_ends >= minPartnerEnds_BothExonJunction 
										& (num_partner_ends >= num_unique_reads | num_partner_ends >= minPartnerEnds_OneExonJunction))))
							 | (exon.junction == 'One' 
							    & frame != 'out-frame'
							    & num_partner_ends >= minPartnerEnds_OneExonJunction
							    & (known ==1 | (intragene==0 & frame =='in-frame')))
					))[, 1:15]

			# Output 1st of records with same GeneExon5_GeneExon3
			ex1 = ex[order(ex$frame, -ex$num_partner_ends, -ex$num_unique_reads),]
				write.table(ex1, paste(sampleID, '.brief.summary', sep=''), row.names=F, quote=F, sep='\t')
			
			
	##=======================================
	## export max 10 example fusion reads
	##=======================================
			ex1 = read.table(paste(sampleID, '.brief.summary', sep=''), sep='\t', stringsAsFactors=F, header=T)
			ex2 = ex1[!duplicated(ex1$"GeneExon5_GeneExon3"),]
			# ex2 = subset(ex2, frame != 'gDNA')
			(nex = min(10, nrow(ex2)))
			if (nex >0){
				## show max 25 example reads 
				for (i in 1:nex){
					## get read ID
					(gei = ex2$"GeneExon5_GeneExon3"[i])
					(ngeci = ex2$num_unique_reads[i])
					readids = unique(lr5$readID[lr5$ge1ge2 == gei])
					(sample.n = min(10, ngeci, length(readids)))
					(readid = sample(unique(lr5$readID[lr5$ge1ge2 == gei]), sample.n))
					writeLines(readid, 'tmp.readid')
					system('sed -e "s/:umi.*//" -e "s:/1.*::" -e "s:/2.*::" tmp.readid | sort -u > tmp.readid2')

					## get fa
					if (file.exists(paste0(sampleID, '.consolidated.bam'))){
						bamfile=paste0(sampleID, '.consolidated.bam')
						}else{ bamfile=paste0(bam_dir, "/", sampleID, ".consolidated.bam")}
					system(paste0(samtools, " view ", bamfile, " | grep -f tmp.readid2 | cut -f1,10 | sed 's/^/>/' |\
							 awk '{$3=length($2); print $0}' | sort -k1,1b -k3,3nr | sort -k1,1b -u |\
							 cut -d ' ' -f1,2 | tr ' ' '\n' | sed 's/:umi.*//' > "
							, sampleID, '.', gei, ".txt", sep=''))
					system(paste0(samtools, " view ", bamfile, " | grep -f tmp.readid2 | ", paste0(samtools, " view -T "), refGenome, " -bS - > ", sampleID, '.', gei, ".bam", sep=''))
					system(paste0(samtools, " index ",
                                                        sampleID, '.', gei, ".bam"))
					}
				system('rm tmp.readid*')
				}
		       }
		}
}

##==== make empty table for negative samples
if (!file.exists(paste(sampleID, '.brief.summary', sep=''))){
	file.create(paste(sampleID, '.brief.summary', sep=''))
}
}
