
## For genes annotated with multiple transcripts by ANNOVAR, keep one of them (NM, exon#, cdna#) -- here keep the refGene transript with most exons

annovar.exon.cds.extraction <- function(input){
library("SplitFusion")
source('config.txt')
in.anno = input
refNM = read.table(paste(path.package("SplitFusion"), '/data/', genomeVer, '_refGene.mostExon.NM', sep=''), sep=' ', header=T, stringsAsFactors=F)
refNM2 = refNM[,c('name2', 'name')]
names(refNM2) = c('Gene.refGene', 'refNM')

d1 = read.table(in.anno, sep='\t', header=T, stringsAsFactors=F)
        d1$Gene.refGene = sub(';.*', '', d1$Gene.refGene)
        d1$Gene.refGene = sub(',.*', '', d1$Gene.refGene)
        d2 = merge(d1, refNM2, by='Gene.refGene', all.x=T)
        d2$tmp1 = mapply(sub, pattern=paste('.*', d2$refNM, ':', sep=''), replacement='', x=d2$AAChange.refGene)

            # when the two refGenes (ANNOVAR.refGenes and refGene.mostExon.NM) do not match:
            d2$s1.4 = substr(d2$tmp1, 1, 4)
            d2$notMatch = 0
            d2$notMatch[d2$Func.refGene == 'exonic' & d2$s1.4 != 'exon'] = 1
                    #table(d2$notMatch)
            d2$tmp1[d2$notMatch==1] = sub('.*NM_', 'NM_', d2$AAChange.refGene[d2$notMatch==1])
                    #head(d2[d2$notMatch==1,])
            d2$refNM[d2$notMatch==1] = sub(':.*', '', d2$tmp1[d2$notMatch==1])
            d2$tmp1[d2$notMatch==1] = sub('.*:exon', 'exon', d2$tmp1[d2$notMatch==1])
                    #head(subset(d2, notMatch==1))

        d2$tmp2 = sub(',.*', '', d2$tmp1)
        d2$exon = sub(':.*', '', d2$tmp2)
        d2$tmp3 = sub('.*c.[ACTG]', '', d2$tmp2)
        d2$cdna = sub('[ACTG].*', '', d2$tmp3)
        d2$tmp4 = sub(':p.*', '', d2$tmp3)
        d2$tmp5 = sub("\\d+", '', d2$tmp4)
        d2$geneStrand = ifelse(d2$tmp5 == 'A', '+', '-')
        d2$chrpos = paste(d2$Chr, d2$Start, sep='_')
        d2 = d2[order(d2$chrpos),]
        d2$ExonicFunc.refGene = sub(';.*', '', gsub(' ', '_', d2$ExonicFunc.refGene))
        d2$Func.refGene = sub(';.*', '', d2$Func.refGene)

                #table(d2$Func, d2$exon, useNA='always')
        d2$exon[is.na(d2$exon)] = d2$Func.refGene[is.na(d2$exon)]
        d3 = subset(d2, select=c('chrpos', 'Gene.refGene', 'geneStrand', 'Func.refGene','ExonicFunc.refGene', 'refNM', 'exon', 'cdna'))
                #head(d3)

## same name replace
        oldnames = names(d3)
        d3$Gene.refGene2[d3$Gene.refGene=='MIR548F1'] = 'TPR'
        d3$Func.refGene2[d3$Gene.refGene=='MIR548F1'] = 'intronic'

        d3$Gene.refGene[!is.na(d3$Gene.refGene2)]=d3$Gene.refGene2[!is.na(d3$Gene.refGene2)]
        d3$Func.refGene[!is.na(d3$Func.refGene2)]=d3$Func.refGene2[!is.na(d3$Func.refGene2)]

write.table(d3[,oldnames], paste0(in.anno, '.ext0'), col.names=F, row.names=F, quote=F, sep='\t')

}
