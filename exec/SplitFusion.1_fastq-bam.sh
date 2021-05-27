#!/bin/bash

. config.txt
SampleId=$( pwd | sed "s:.*/::")

	#=== [Kickstart mode] If specify bam_file, start from bam ===
	if [ "$bam_file" != "" ]; then
		    $samtools view -@ $thread $bam_file > _raw.sam
		#=== If specify fastq_file1, then start from fastq to bam ===
		#=== align to genome using bwa mem -q, correct primary alignment mapping quality with secondary alignment ===
	elif [ "$fastq_file1" != "" ]; then
		if [ "$fastq_file2" != "" ]; then
			$bwa mem -T 18 -q -t $thread $refGenome $fastq_file1 $fastq_file2 2> bwa.log | awk 'BEGIN{OFS="\t"}{if($2<2048&&$5==0){for(i=1;i<=NF;i++){if($i~/^SA:Z/){split($i,a,",");if(length(a)==6&&a[5]>0){$5=int(a[5]/2);break}}}}print}' > _raw.sam
		else	
			$bwa mem -T 18 -q -t $thread $refGenome $fastq_file1 2> bwa.log | awk 'BEGIN{OFS="\t"}{if($2<2048&&$5==0){for(i=1;i<=NF;i++){if($i~/^SA:Z/){split($i,a,",");if(length(a)==6&&a[5]>0){$5=int(a[5]/2);break}}}}print}' > _raw.sam
		fi
	else 
		echo "Must specify fastq_file or bam_file"
		exit
	fi

	head -n 10000 _raw.sam | grep ^@ > header

#== Consolidate reads based on unique UMI

	#==== if no umi, add umi:A## for compatability, where ## is the read pair number starting from 1.
        hasUmi=$(grep -v ^@ _raw.sam | head -n 1 | cut -f 1 | grep umi: | wc -l)
		if [ $hasUmi -eq 0 ]; then
			grep -v ^@ _raw.sam | awk '{OFS="\t";if(preId!=$1){n++}; preId=$1; $1=$1":umi:A"n; print $0}' > _raw.sam2
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
	cp header consolidated.sam
		sed -e 's:\t\t:\t*\t:g' _raw.sam | sed -e 's/umi:/umi\t/' | sort --parallel=$thread -k1,1b | \
		join _consolidated.readID - | sed -e 's/ /:/' | tr ' ' '\t' | cut -f1,3- >> consolidated.sam
	
rm _*
grep -P '\tSA:Z:' consolidated.sam > _sa.sam

$samtools view -@ $thread -T $refGenome -bS consolidated.sam | $samtools sort -@ $thread -o $SampleId.consolidated.bam -
$samtools index $SampleId.consolidated.bam 
