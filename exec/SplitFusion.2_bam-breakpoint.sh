#!/bin/bash

. config.txt
SampleId=$( pwd | sed "s:.*/::")

##==== 1.1. get reads with SA
	if [ ! -s _sa.sam ]; then
		$samtools view -@ $thread $SampleId.consolidated.bam | awk -F "\t" 'BEGIN{OFS="\t"}/\tSA:Z/{for(i=1;i<=NF;i++){if($i~/^SA:Z/){n=split($i,a,",");if(n==6){if(a[5]>$5){$5=int(($5+a[5])/2)}print;break}}}}' > _sa.sam
	fi

	$samtools view -@ $thread -T $refGenome -bS _sa.sam > _sa.bam
	$bedtools bamtobed -cigar -i _sa.bam > _sa.bed0

	## bedtools uses 0-base, change to 1-base:
	awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed 

if [ -s _sa.bed ]; then

##==== 1.2. get read length
        awk '{n=length($10); print $1,n}' _sa.sam > _sa.len

##==== 2. Join bed and read length, and keep the longest
	# By default, number and order of reads in files _sa.sam, _sa.bam and _sa.bed0 are identical. So, can use paste.
    	paste _sa.len _sa.bed | tr ' ' '\t' | cut -f2- > _sa.len.bed
		
	cut -f1,5 _sa.len.bed > _len.readID1
	sort --parallel=$thread -k2,2b -k1,1nr _len.readID1 > _len.readID2
	sort --parallel=$thread -k2,2b -u _len.readID2 > _len.readID

	sort --parallel=$thread -k2,2b _len.readID > _len.readID.s
	sort --parallel=$thread -k4,4b _sa.bed > _sa.bed.s
	join -1 2 -2 4 _len.readID.s _sa.bed.s > _sa.len.bed

##==== 3. remove MQ low
	echo | awk -v minMQ=$minMQ '{OFS="\t"; if ($6 >= minMQ) print $0}' _sa.len.bed > _sa.len.bed.mq
#        awk -v minMQ=$minMQ '{OFS="\t"; for(i=1;i<=NF;i++){if($i~/^SA:Z/){split($i,a,",");if(length(a)==6&&a[5]>0){$5=int(a[5]/2);break}} if ($5 >= minMQ) print $0}' _sa.len.bed > _sa.len.bed.mq
##==== 4. calculate query start, end
	# corrected read length for small size clipping which is either
	#===============================#
	#=== transform CIGAR strings ===#
	#===============================#
	cut -f 2,7,8 -d ' ' _sa.len.bed.mq > _sa.SMH1

	sed -e 's/M/ M /g' -e 's/S/ S /g' -e 's/H/ H /g' _sa.SMH1 > _sa.SMH2

	    # correct D in CIGAR (1)
	    awk '{if ($5 ~ /1D/){
			    lenM = $3 + substr($5, 1, 1) + substr($5, 3, 3)
			    $1= $1 + substr($5, 1, 1); 
			    $3 = lenM; $5="_"; $6="_";
			    }
		    print $0
		    }' _sa.SMH2 | sed 's/_ //g' > _sa.SMH2a

	    # correct D in CIGAR (2)
	    awk '{if ($7 ~ /1D/){
			    lenM = $5 + substr($7, 1, 1) + substr($7, 3, 3)
			    $1= $1 + substr($7, 1, 1); 
			    $5 = lenM; $7="_"; $8="_";
			    }
		    print $0
		    }' _sa.SMH2a | sed 's/_ //g' > _sa.SMH2b

	    # correct I in CIGAR (1)
	    awk '{if ($5 ~ /1I/){
			    lenM = $3 - substr($5, 1, 1) + substr($5, 3, 3)
			    $1= $1 - substr($5, 1, 1);
			    $3 = lenM; $5="_"; $6="_";
			    }
		    print $0
		    }' _sa.SMH2b | sed 's/_ //g' > _sa.SMH2c

	    # correct I in CIGAR (2)
	    awk '{if ($7 ~ /1I/){
			    lenM = $5 - substr($7, 1, 1) + substr($7, 3, 3)
			    $1= $1 - substr($7, 1, 1); 
			    $5 = lenM; $7="_"; $8="_";
			    }
		    print $0
		    }' _sa.SMH2c | sed 's/_ //g' > _sa.SMH2d


	awk '{if ($4=="M"){
			if ($2=="+"){
				start=1; end=$3
			}else if ($2=="-"){
				start=$1-$3+1; end=$1
			}
		} else if ($6=="M"){
			if ($2=="+"){
				start=$3+1; end=$3+$5
			}else if ($2=="-"){
				start=$1-$3-$5+1; end=$1-$3
			}
		} else {start=0; end=0}; print start,end}' _sa.SMH2d > _sa.SMH3
	paste _sa.len.bed.mq _sa.SMH3 | tr ' ' '\t' > _sa.SMH4
	sort --parallel=$thread -k1,1b -k9,9n _sa.SMH4 > _sa.SMH4s

	#=== Correcting ligate.UMI based on Read1 head or, when Read1 is not mapped, Read2 tail
	awk '{n=split($1,a,"-");sub(/:umi:/,"\t",a[1]);printf "%s",a[1];for(i=2;i<n;i++){printf "-%s",a[i]}printf "\t%s",a[n];for(i=2;i<=NF;i++){printf "\t%s",$i}printf "\n"}' _sa.SMH4s > _corr.ligat1
	awk '{OFS="\t"; if ($3 ~ /\/1/) {order=$11} else {order=-$12}; print $0,order}' _corr.ligat1 > _corr.ligat1b
	sort --parallel=$thread -k1,1b -k3,3b -k13,13n _corr.ligat1b > _corr.ligat1s
	sort --parallel=$thread -k1,1b -u _corr.ligat1s > _corr.ligat2

	echo | awk -v minMapLength=$minMapLength '{OFS="\t"; 
		if ($3 ~ /\/1/){
			if ($11 > minMapLength) {$11=1}
		    	if ($9 == "+"){posC = 100000002 + $6 - $11
				}else{ posC = 100000000 + $7 + $11}
			}
		else { 	if ($9 == "-"){posC = 100000002 + $6 - $11
				}else{ posC = 100000000 + $7 + $11}
			}
                       
		umi="C"$5"P"posC;
		$2=umi;
		print $1,$2
		}' _corr.ligat2 > corr.ligat3

	join corr.ligat3 _corr.ligat1s | tr ' ' '\t' | cut -f1,2,4-13 | sed -e "s/\t/:umi:/" -e "s/\t/-/" > _sa.SMH.corr
	sort -k1,1b -k9,9n _sa.SMH.corr > _sa.SMH.corr.srt

##==== 5. separate left and right split alignments
	awk '{if ($1==pre1){
			n=n+1
		} else {n=1};

		if (n==1) {
			print $0 > "_left";
			if (pren !=1) {print pre0 > "_right"}
		} else if (n>2){
			print pre0 > "split.mid"
		};
		pre1=$1; pren=n; pre0=$0;
		print n,$0 > "_sa.SMH4sn"
	}' _sa.SMH.corr.srt

	# remove ligation.site not near Read 1 _Left
	grep "/1" _left | sed "s/:umi:C/\t/" | sed -e "s/P/\t/" -e "s/-/\t/" | awk '{OFS="\t";
			diff = $3 - 100000000 - $7;
			if (diff <0) {diff = -diff};
			if ($2 != $6 || diff > 750000){
				print $1":umi:C"$2"P"$3"-"$4 > "_diff_ligate_left"
			}
		}'
	
		if [ -s _diff_ligate_left ]; then
			join -v 1 _left _diff_ligate_left > _left2
		else 
			cp _left _left2
		fi; 

##==== 6. get the 4 positions on a SA query read (after left and right alignments are merged by ReadID):
		## skematic drawing:
		##  left.query.start....left.query.end-[breakpoint]-right.query.start....right.query.end
		##         1                   2      -[breakpoint]-        3                  4
			## 2: left.query.end [corresponding mapping position in the genome is breakpoint 1]
			## 3: right.query.start [corresponding mapping position in the genome is breakpoint 2]
	        ## add breadkpoint (chr and pos)
            awk '{if ($7 =="+") {
                                print $0,$3"_"$5
                        } else {
                                mstart=$5; mend=$4; $4=mstart; $5=mend;
                                print $0,$3"_"$5
                        }
                }' _left2 | tr ' ' '\t' > _lefts

            awk '{if ($7 =="+") {
                                print $0,$3"_"$4
                        } else {
                                mstart=$5; mend=$4; $4=mstart; $5=mend;
                                print $0,$3"_"$4
                        }
                }' _right | tr ' ' '\t' | sed 1d > _rights

	##=== join leftmost and rightmost
		sort --parallel=$thread -k1,1b _lefts > _lefts.s
		sort --parallel=$thread -k1,1b _rights > _rights.s
        join _lefts.s _rights.s > _left_right.j

##==== 7. breakpoint.candidates.preFilter
	## left: 2-11; right 12-21
	awk '{if ($11 < $21){
			pp=$11"__"$21
		} else {pp=$21"__"$11
		}; 
		print $0,pp
	}' _left_right.j > _breakpoint.noFilter1
	sort --parallel=$thread -k1,1b _breakpoint.noFilter1 > _breakpoint.noFilter2
else
	touch breakpoint.candidates.preFilter
fi
##==== 8. correct breakpoint for those containing mid
if [ -f split.mid ]; then
	cut -f1 split.mid | sort --parallel=$thread -u > _mid.id
	sort --parallel=$thread -k2,2b _sa.SMH4sn > _sa.SMH4sns
	join -1 1 -2 2 _mid.id _sa.SMH4sns > split.mid.expanded

	if [ -s split.mid.expanded ]; then
		awk '{if ($1==pre1){
			diff = $5 - pre5;
			if ($4 != pre4 || diff > 100000 || diff < -100000){
				if (pre8=="+"){bkp1=pre4"_"pre6} else {bkp1=pre4"_"pre5};
				if ($8=="+"){bkp2=$4"_"$5} else {bkp2=$4"_"$6}
				q1=pre10; q2=pre11
			   };
			};
		if (bkp1 < bkp2){bkp = bkp1"__"bkp2};
		if (bkp1 > bkp2){bkp = bkp2"__"bkp1};
		if (bkp != ""){print $1,bkp,q1,q2 > "_sa.mid.bkp"}

		diff=""; bkp=""; bkp1=""; bkp2=""; q1=""; q2=""
		pre1=$1; pre4=$4; pre5=$5; pre6=$6; pre8=$8; pre10=$10; pre11=$11
		}' split.mid.expanded

		sort --parallel=$thread -k1,1b _sa.mid.bkp > _sa.mid.bkps
		join -a1 _breakpoint.noFilter2 _sa.mid.bkps > _breakpoint.noFilter3

		awk '{if (NF==25) {$9=$24; $10=$25; $22=$23}; print}' _breakpoint.noFilter3 | cut -d ' ' -f 1-22 > _breakpoint.noFilter.bkp.corrected

		join _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.candidates.preFilter.w.mid
		join -v 1 _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.candidates.preFilter
	else 
        	cp _breakpoint.noFilter2 breakpoint.candidates.preFilter
	fi
else
        cp _breakpoint.noFilter2 breakpoint.candidates.preFilter
fi

# DONE breakpoint candidate
