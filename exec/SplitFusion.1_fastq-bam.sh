#!/bin/bash

. config.txt
SampleId=$( pwd | sed "s:.*/::")
memG=$(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)))G
awk=$(which mawk)
mawkver=$($awk -W version |& head -n 1|cut -d ' ' -f2)
if [[ "$mawkver" != "1.3.4" || "$awk" == "" ]]; then
  awk=$(which gawk)
  if [ "$awk" == "" ]; then
    echo "Neither mawk 1.3.4 nor gawk is installed in the system!"
    exit
  fi
fi

	#=== [Kickstart mode] If specify bam_file, start from bam ===
	if [ "$bam_file" != "" ]; then
                    hasUmi=$($(samtools) view $bam_file | head -n 1 | cut -f 1 | grep umi: | wc -l)
		    $samtools view -@ $thread $bam_file > _raw.sam
		#=== If specify fastq_file1, then start from fastq to bam ===
		#=== align to genome using bwa mem -q, such that secondary alignment will display original mapping quality ===
	elif [ "$fastq_file1" != "" ]; then
		if [[ "$fastq_file1" == *.gz ]]; then
                	hasUmi=$(zcat $fastq_file1 | head -n 1 | cut -f 1 | grep umi: | wc -l)
		else
			hasUmi=$(cat $fastq_file1 | head -n 1 | cut -f 1 | grep umi: | wc -l)
		fi
		if [ "$fastq_file2" != "" ]; then
	#==== if no umi, add umi:A for compatability
			if [ $hasUmi -eq 0 ]; then
				$bwa mem -T 18 -q -K 10000000 -t $thread $refGenome $fastq_file1 $fastq_file2 2> bwa.log | $awk '/^@/{print}!/^@/{OFS="\t";if(preId!=$1){n++}; preId=$1; $1=$1":umi:A"n; print $0}' | $samtools view -@ $thread -bS - > _raw.bam
			else
				$bwa mem -T 18 -q -K 10000000 -t $thread $refGenome $fastq_file1 $fastq_file2 2> bwa.log | $samtools view -@ $thread -bS - > _raw.bam
			fi
		else	
			if [ $hasUmi -eq 0 ]; then
				$bwa mem -T 18 -q -K 10000000 -t $thread $refGenome $fastq_file1 2> bwa.log | $awk '/^@/{print}!/^@/{OFS="\t";if(preId!=$1){n++}; preId=$1; $1=$1":umi:A"n; print $0}' | $samtools view -@ $thread -bS - > _raw.bam
			else
				$bwa mem -T 18 -q -K 10000000 -t $thread $refGenome $fastq_file1 2> bwa.log | $samtools view -@ $thread -bS - > _raw.bam

			fi
		fi
	else 
		echo "Must specify fastq_file or bam_file"
		exit
	fi

#	$samtools view -H _raw.bam > consolidated.sam

#== Consolidate reads based on unique UMI

	ligateUmiFormat=$($samtools view _raw.bam | head -n 1 | cut -f 1 | grep umi: | sed 's:.*umi:umi:' | grep C | grep P | wc -l)

	#==== Tag chr.pos (at ligation site) to umi
#                $bedtools bamtobed -cigar -i _raw.bam > _raw.bed

		#==== Reformat if not already in the required ligate UMI format, which 
		#====  has chr.pos (of ligation site) integrated in UMI (some 
		#====  read with suppl alignments need correction for ligate site later)
		    if [ $ligateUmiFormat -eq 0 ]; then

			$bedtools bamtobed -cigar -i _raw.bam | sed -e 's/:umi:/\tumi\t/' |\
			$awk '{OFS="\t";
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
	 
				if ($6 !~ /N/){
					print $4":umi",umi"-"$6
				}
			}' 2>/dev/null | sed 's:/[12]$::' | sort --parallel=$thread -k2,2b -u -S $memG > uniq.ligateUmi
		    else
			   $bedtools bamtobed -cigar -i _raw.bam | cut -f4 | sed 's/umi:/umi\t/' | sed 's:/[12]$::' | sort --parallel=$thread -k2,2b -u -S $memG > uniq.ligateUmi
		    fi
     
	#==== consolidation
	#=== 	join raw sam with consolidated ID ===
#		$samtools view -@ $thread _raw.bam | sed -e 's:\t\t:\t*\t:g' | sed -e 's/umi:/umi\t/' | sort --parallel=$thread -k1,1b -S $memG | \
#		join <(sort --parallel=$thread -k1,1b -u -S $memG uniq.ligateUmi) - | sed -e 's/ /:/' | tr ' ' '\t' | cut -f1,3- | $samtools view -@ $thread -T $refGenome -bS - | $samtools sort -@ $thread -m $sortmemG -o $SampleId.consolidated.bam -

# Calculate number of threads and RAM per thread for samtools sort
# Assume 2GB RAM is used for OS and other programs
# Roughly, RAM around the size of _raw.bam is used for the awk step
# First, calculate number of KB RAM per thread after subtracting the
# above two numbers.
# If the RAM per thread calculated is larger than default 768MB,
# then use the calculated value and also set thread to original value
# Otherwise, reduce number of threads to the point that RAM per
# thread is larger than 768MB
totalmem=$(($(getconf _PHYS_PAGES)*$(getconf PAGE_SIZE)))
sortmem=$((($totalmem - $(du -k _raw.bam | cut -f1)*1024 - 2147483648)/ (1024 * $(echo $thread))))
if [[ $sortmem -gt 768000 ]]; then
	sortmem=$(echo "-@ $thread -m ${sortmem}K")
else
        newthread=$((($totalmem - $(du -k _raw.bam | cut -f1)*1024 - 2147483648)/ (1024 * 768000)))
        sortmem=$((($totalmem - $(du -k _raw.bam | cut -f1)*1024 - 2147483648)/ (1024 * $(echo $newthread))))
        sortmem=$(echo "-@ $newthread -m ${sortmem}K")
fi

	$samtools view -h -@ $thread _raw.bam | $awk -F"\t" 'BEGIN{OFS="\t";while(getline<"uniq.ligateUmi"){a[$1]=$2}}{if(/^@/){print}else{sub(":umi.*",":umi",$1);if(a[$1]!=""){$1=sprintf("%s:%s",$1,a[$1]);print}}}' | $samtools sort $sortmem -o $SampleId.consolidated.bam -	
rm _raw.bam
#grep -P '\tSA:Z:' consolidated.sam > _sa.sam

#$samtools view -@ $thread -T $refGenome -bS consolidated.sam | $samtools sort -@ $thread -o $SampleId.consolidated.bam -
$samtools index $SampleId.consolidated.bam 
#rm consolidated.sam
