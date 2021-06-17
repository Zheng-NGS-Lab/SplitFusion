#!/bin/bash

. config.txt
SampleId=$( pwd | sed "s:.*/::")
memG=$(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)))G
awk=$(which mawk)
if [ "$awk" == "" ]; then
  awk=$(which gawk)
  if [ "$awk" == "" ]; then
    echo "Neither mawk nor gawk is installed in the system!"
    exit
  fi
fi

##==== 1.1. get reads with SA
##==== 1.2. average the mapping quality of primary and secondary alignment
##==== 1.3. print alignment length to _sa.len
##==== 1.4. bedtools bamtobed output to _sa.bed

	$samtools view -@ $thread $SampleId.consolidated.bam | $awk -F "\t" 'BEGIN{OFS="\t"}/\tSA:Z/{for(i=1;i<=NF;i++){if($i~/^SA:Z/){n=split($i,a,",");if(n==6){if(a[5]>$5){$5=int(($5+a[5])/2)}print;printf "%s %d\n",$1,length($10) > "_sa.len";break}}}}' | $samtools view -T $refGenome -bS - | $bedtools bamtobed -cigar -i stdin | $awk -F"\t" '{start = $2+1; $2=start; print}' > _sa.bed

if [ -s _sa.bed ]; then

##==== 2. Join bed and read length, and keep the longest
##==== 3. Filter by mapping quality
	paste -d ' ' _sa.len _sa.bed | cut -d ' ' -f2- | cut -d ' ' -f1,5 | sort --parallel=$thread -k2,2b -k1,1nr -S $memG | sort --parallel=$thread -k2,2b -u -S $memG | sort --parallel=$thread -k2,2b -S $memG | join -1 2 -2 4 - <(sort --parallel=$thread -k4,4b -S $memG _sa.bed) | $awk -v minMQ=$minMQ '{OFS="\t"; if ($6 >= minMQ) print $0}'> _sa.len.bed.mq

##==== 4. calculate query start, end

	# corrected read length for small size clipping which is either
	#===============================#
	#=== transform CIGAR strings ===#
	#===============================#
	cut -f 2,7,8 -d ' ' _sa.len.bed.mq | sed -e 's/M/ M /g' -e 's/S/ S /g' -e 's/H/ H /g' | \
	    # correct D in CIGAR (1)
	    $awk '{if ($5 ~ /1D/){
			    lenM = $3 + substr($5, 1, 1) + substr($5, 3, 3)
			    $1= $1 + substr($5, 1, 1); 
			    $3 = lenM; $5="_"; $6="_";
			    }
		    print $0
		    }' | sed 's/_ //g' | \
	    # correct D in CIGAR (2)
	    $awk '{if ($7 ~ /1D/){
			    lenM = $5 + substr($7, 1, 1) + substr($7, 3, 3)
			    $1= $1 + substr($7, 1, 1); 
			    $5 = lenM; $7="_"; $8="_";
			    }
		    print $0
		    }' | sed 's/_ //g' | \
	    # correct I in CIGAR (1)
	    $awk '{if ($5 ~ /1I/){
			    lenM = $3 - substr($5, 1, 1) + substr($5, 3, 3)
			    $1= $1 - substr($5, 1, 1);
			    $3 = lenM; $5="_"; $6="_";
			    }
		    print $0
		    }' | sed 's/_ //g' | \
	    # correct I in CIGAR (2)
	    $awk '{if ($7 ~ /1I/){
			    lenM = $5 - substr($7, 1, 1) + substr($7, 3, 3)
			    $1= $1 - substr($7, 1, 1); 
			    $5 = lenM; $7="_"; $8="_";
			    }
		    print $0
		    }' | sed 's/_ //g' | \
	$awk '{if ($4=="M"){
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
		} else {start=0; end=0}; print start,end}' | \
	#=== Correcting ligate.UMI based on Read1 head or, when Read1 is not mapped, Read2 tail
	paste _sa.len.bed.mq - | tr ' ' '\t' |sort --parallel=$thread -k1,1b -k9,9n -S $memG | $awk '{n=split($1,a,"-");sub(/:umi:/,"\t",a[1]);printf "%s",a[1];for(i=2;i<n;i++){sub(/:umi:/,"\t",a[i]);printf "-%s",a[i]}printf "\t%s",a[n];for(i=2;i<=NF;i++){printf "\t%s",$i}printf "\n"}' | $awk '{OFS="\t"; if ($3 ~ /\/1/) {order=$11} else {order=-$12}; print $0,order}' | sort --parallel=$thread -k1,1b -k3,3b -k13,13n -S $memG > _corr.ligat1s
	sort --parallel=$thread -k1,1b -u -S $memG _corr.ligat1s | $awk -v minMapLength=$minMapLength '{OFS="\t"; 
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
		}' > corr.ligat3

	join corr.ligat3 _corr.ligat1s | tr ' ' '\t' | cut -f1,2,4-13 | sed -e "s/\t/:umi:/" -e "s/\t/-/" | sort -k1,1b -k9,9n -S $memG | \
##==== 5. separate left and right split alignments
	$awk '{if ($1==pre1){
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
	}'

	# remove ligation.site not near Read 1 _Left
	grep "/1" _left | sed "s/:umi:C/\t/" | sed -e "s/P/\t/" -e "s/-/\t/" | $awk '{OFS="\t";
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
            $awk '{if ($7 =="+") {
                                print $0,$3"_"$5
                        } else {
                                mstart=$5; mend=$4; $4=mstart; $5=mend;
                                print $0,$3"_"$5
                        }
                }' _left2 | tr ' ' '\t' | sort --parallel=$thread -k1,1b -S $memG > _lefts.s

	##=== join leftmost and rightmost
##==== 7. breakpoint.candidates.preFilter
	## left: 2-11; right 12-21
            $awk '{if ($7 =="+") {
                                print $0,$3"_"$4
                        } else {
                                mstart=$5; mend=$4; $4=mstart; $5=mend;
                                print $0,$3"_"$4
                        }
                }' _right | tr ' ' '\t' | sed 1d | sort --parallel=$thread -k1,1b -S $memG | join _lefts.s - | \
	$awk '{if ($11 < $21){
			pp=$11"__"$21
		} else {pp=$21"__"$11
		}; 
		print $0,pp
	}' | sort --parallel=$thread -k1,1b -S $memG > _breakpoint.noFilter2
        rm _left _left2 _lefts.s _right _corr.ligat1s _sa.len.bed.mq _sa.bed _sa.len
else
	rm _sa.len
	touch breakpoint.candidates.preFilter
fi
##==== 8. correct breakpoint for those containing mid
if [ -f split.mid ]; then
	cut -f1 split.mid | sort --parallel=$thread -u -S $memG > _mid.id
	sort --parallel=$thread -k2,2b -S $memG _sa.SMH4sn | join -1 1 -2 2 _mid.id - > split.mid.expanded

	if [ -s split.mid.expanded ]; then
		$awk '{if ($1==pre1){
			diff = $5 - pre5;
			if ($4 != pre4 || diff > 100000 || diff < -100000){
				if (pre8=="+"){bkp1=pre4"_"pre6} else {bkp1=pre4"_"pre5};
				if ($8=="+"){bkp2=$4"_"$5} else {bkp2=$4"_"$6}
				q1=pre10; q2=pre11
			   };
			};
		if (bkp1 < bkp2){bkp = bkp1"__"bkp2};
		if (bkp1 > bkp2){bkp = bkp2"__"bkp1};
		if (bkp != ""){print $1,bkp,q1,q2}

		diff=""; bkp=""; bkp1=""; bkp2=""; q1=""; q2=""
		pre1=$1; pre4=$4; pre5=$5; pre6=$6; pre8=$8; pre10=$10; pre11=$11
		}' split.mid.expanded | sort --parallel=$thread -k1,1b -S $memG | join -a1 _breakpoint.noFilter2 - | $awk '{if (NF==25) {$9=$24; $10=$25; $22=$23}; print}' | cut -d ' ' -f 1-22 > _breakpoint.noFilter.bkp.corrected

		join _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.candidates.preFilter.w.mid
		join -v 1 _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.candidates.preFilter
		rm _breakpoint.noFilter.bkp.corrected _mid.id _breakpoint.noFilter2 _sa.SMH4sn
	else
		rm _mid.id _sa.SMH4sn
        	mv _breakpoint.noFilter2 breakpoint.candidates.preFilter
	fi
else
	rm _sa.SMH4sn
        mv _breakpoint.noFilter2 breakpoint.candidates.preFilter
fi

# DONE breakpoint candidate
