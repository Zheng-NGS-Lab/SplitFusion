#!/bin/bash

. config.txt
SampleId=$( pwd | sed "s:.*/::")

	#=== [Kickstart mode] If specify bam_dir, start from bam ===
	if [ "$bam_dir" != "" ]; then
		if [ -s $bam_dir/$SampleId.bam ]; then
		    $samtools view -@ $thread $bam_dir/$SampleId.bam > _raw.sam
		elif [ -s $bam_dir/$SampleId.consolidated.bam ]; then
		    $samtools view -@ $thread $bam_dir/$SampleId.consolidated.bam > _raw.sam
		fi

	#=== If specify fastq_dir, do fasq to bam ===
	elif [ "$fastq_dir" != "" ]; then
		if [ "$r1filename" != "" ]; then
			$bwa mem -T 18 -t $thread $refGenome $fastq_dir/$r1filename $fastq_dir/$r2filename > _raw.sam 2> bwa.log
		else	
			$bwa mem -T 18 -t $thread $refGenome $fastq_dir/$SampleId.R1.fq $fastq_dir/$SampleId.R2.fq > _raw.sam 2> bwa.log
		fi
	else 
		echo "Must specify fastq_dir or bam_dir"
		exit
	fi

	head -n 10000 _raw.sam | grep ^@ > header

#== Consolidate reads based on unique UMI

	#==== if no umi, add umi:A for compatability
        hasUmi=$(grep -v ^@ _raw.sam | head -n 1 | cut -f 1 | grep umi: | wc -l)
		if [ $hasUmi -eq 0 ]; then
			grep -v ^@ _raw.sam | sed 's/\t/:umi:A\t/' > _raw.sam2
			mv _raw.sam2 _raw.sam
		fi

	$samtools view -@ $thread -T $refGenome -bS _raw.sam > _raw.bam

	#==== Tag chr.pos (at ligation site) to umi
                $bedtools bamtobed -cigar -i _raw.bam > _raw.bed

		#==== Reformat if not already in the required ligate UMI format, which 
		#====  has chr.pos (of ligation site) integrated in UMI (some 
		#====  read with suppl alignments need correction for ligate site later)
		ligateUmiFormat=$(head -n 1 _raw.bed | cut -f 4 | grep umi: | sed 's:.*umi:umi:' | grep C | grep P | wc -l)
		    if [ $ligateUmiFormat -eq 0 ]; then

			sed -e 's/:umi:/\tumi\t/' _raw.bed |\
			gawk '{OFS="\t";
				if ($4 != preID){
					if ($8 == "+"){
						posC = 100000001 + $2
					} else {
						posC = 100000001 + $3
					}
					umi="C"$1"P"posC
				} else {
					umi=preUmi
				};
	 
				preUmi=umi;
				preID=$4;
	 
				if ($6 ~ /N/){
					print $4":umi",umi"-"$6 > "_id.N.umi"
				} else {
					print $4":umi",umi"-"$6 > "_id.ligateUmi"
				}
			}' 2>/dev/null
		    else
			   cut -f4 _raw.bed | sed 's/umi:/umi\t/' > _id.ligateUmi 
		    fi
     
	#==== consolidation
	sed 's:/[12]$::' _id.ligateUmi >  _id.ligateUmi1
	sort --parallel=$thread -k2,2b -u _id.ligateUmi1 > uniq.ligateUmi
	sort --parallel=$thread -k1,1b -u uniq.ligateUmi > _consolidated.readID

	#=== 	join raw sam with consolidated ID ===
		sed -e 's:\t\t:\t*\t:g' _raw.sam | sed -e 's/umi:/umi\t/' > _raw.samC
		sort --parallel=$thread -k1,1b _raw.samC > _raw.sam.s
			    rm _raw.samC
     
	join _consolidated.readID _raw.sam.s > _consolidated.sam0
		cp header consolidated.sam
		sed -e 's/ /:/' _consolidated.sam0 | tr ' ' '\t' | cut -f1,3- >> consolidated.sam
		    rm _raw.sam.s
	
rm _*
grep -P '\tSA:Z:' consolidated.sam > _sa.sam
 
$samtools view -@ $thread -T $refGenome -bS consolidated.sam > _consolidated.bam
$samtools sort -@ $thread _consolidated.bam -o $SampleId.consolidated.bam
	rm consolidated.sam _consolidated.bam
$samtools index $SampleId.consolidated.bam 

