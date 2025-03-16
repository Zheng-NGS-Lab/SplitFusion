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

##========================================
## If exists TARGET panel information, use the SplitFusion-TARGET model
##========================================

if [ $panel != "NA" ]; then
	sed "s:umi:\tumi:" breakpoint.candidates.preFilter | sort --parallel=$thread -k1,1b -S $memG > _breakpoint.candidates.preFilter.srt
	$awk '{if ($1 ~ /\/2/ && $9 == 1) print}' breakpoint.candidates.preFilter | $awk '{OFS="\t"; if ($7 == "-"){start=$5; $5=$4; $4=start};print $3,$4,$5,$1,$6,$7}' | sort --parallel=$thread -k1,1n -k2,2n -S $memG | $bedtools intersect -wa -wb -a $panel_dir/${genomeVer}_${panel}.GSP2.bed -b - | $awk '{diff1 = $8 - $2; diff2 = $9 - $3;if ($6 == "-"){if (diff1 >= -1 && diff1 <=1) {if (diff2 >= 3) {print} else {print > "_preFilter.anchored.offTarget.bed"}} else {print > "_preFilter.nonAnchored.bed"}} else {if (diff2 >= -1 && diff2 <=1) {if (diff1 <= -3) {print} else {print > "_preFilter.anchored.offTarget.bed"}} else {print > "_preFilter.nonAnchored.bed"}}}' | cut -f 10 | sed "s:umi.*::" | sort --parallel=$thread -k1,1b -u -S $memG | join - _breakpoint.candidates.preFilter.srt | sed 's: ::' > breakpoint.candidates.preFilter
	rm _breakpoint.candidates.preFilter.srt _preFilter.anchored.offTarget.bed _preFilter.nonAnchored.bed


	##==== reads with middle split
	if [ -f breakpoint.candidates.preFilter.w.mid ]; then	
	    sed "s:umi:\tumi:" breakpoint.candidates.preFilter.w.mid | sort --parallel=$thread -k1,1b -S $memG > _breakpoint.candidates.preFilter.srt
	    $awk '{if ($1 ~ /\/2/ && $9 == 1) print}' breakpoint.candidates.preFilter.w.mid | $awk '{OFS="\t"; if ($7 == "-"){start=$5; $5=$4; $4=start};print $3,$4,$5,$1,$6,$7}' | sort --parallel=$thread -k1,1n -k2,2n -S $memG | $bedtools intersect -wa -wb -a $panel_dir/${genomeVer}_${panel}.GSP2.bed -b - | $awk '{diff1 = $8 - $2; diff2 = $9 - $3;if ($6 == "-"){if (diff1 >= -1 && diff1 <=1) {if (diff2 >= 3) {print} else {print > "_preFilter.anchored.offTarget.bed"}} else {print > "_preFilter.nonAnchored.bed"}} else {if (diff2 >= -1 && diff2 <=1) {if (diff1 <= -3) {print} else {print > "_preFilter.anchored.offTarget.bed"}} else {print > "_preFilter.nonAnchored.bed"}}}' | cut -f 10 | sed "s:umi.*::" | sort --parallel=$thread -k1,1b -u -S $memG | join - _breakpoint.candidates.preFilter.srt | sed 's: ::' > breakpoint.candidates.preFilter.w.mid
	    rm _breapoint.candidates.preFilter.srt _preFilter.anchored.offTarget.bed _preFilter.nonAnchored.bed
	fi

	# trun off minMapLength filters by set to 0
	minMapLength=0; minMapLength2=0; 
fi


##========================================
##==== Filter 1: minMapLength, minExclusive, maxQueryGap, maxOverlap ====:
##========================================
		## 1. at least minMapLength 
		## 2. at least minExclusive bp mapped exculsively to alignments 1 and 2
		## 3. no larger than maxQueryGap
		## 4. no larger than max overlapping length

		## reads without middle split
		    ## left:  $9-------$10
		    ## right:       $19------$20
		$awk -v minMapLength=$minMapLength -v minMapLength2=$minMapLength2 -v minExclusive=$minExclusive -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap -v minMQ1=$minMQ1 \
		    '{ gap = $19-$10-1; overlap = $10-$19+1;
			if ($1 ~ /\/1/)	{mapLen1 = $10-$9+1; mapLen2 = $20-$19+1; mq1=$6
				} else { mapLen2 = $10-$9+1; mapLen1 = $20-$19+1; mq1=$16};
			if ((mapLen1 >= minMapLength && mapLen2 >= minMapLength2 && mq1 >= minMQ1) \
				&& ($19-$9 >= minExclusive && $20-$10 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			    ) {print $0,overlap} else {print $1 > "/dev/null"}
			}' breakpoint.candidates.preFilter > _sa.fu01

		## For reads with middle split, turn off maxQueryGap by let maxQueryGap=1000
		    ## left:  $9-------$10
		    ## right:		       $19------$20
		if [ -f breakpoint.candidates.preFilter.w.mid ]; then
			touch _filter2
		$awk -v minMapLength=$minMapLength -v minMapLength2=$minMapLength2 -v minExclusive=$minExclusive -v maxQueryGap=1000 \
		    '{ gap = $19-$10-1; overlap = $10-$19+1;
			if ($1 ~ /\/1/) {mapLen1 = $10-$9+1; mapLen2 = $20-$19+1
                                } else { mapLen2 = $10-$9+1; mapLen1 = $20-$19+1};

			if (  (mapLen1 >= minMapLength && mapLen2 >= minMapLength2) \
					&& ($19-$9 >= minExclusive && $20-$10 >= minExclusive && gap <= maxQueryGap) \
			    ) {print $0,overlap} else {print $1 > "_filter2"}
			}' breakpoint.candidates.preFilter.w.mid > _sa.fu02

			join -v 1 _sa.fu02 _filter2 > _sa.fu02f

			# filter gap and overlap
			$awk -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap \
			'{ if ($1==pre1){
					gap = $10 - pre11 -1; overlap = pre11 - $10 + 1
					if (gap > maxQueryGap || overlap > maxOverlap) {print $1}
				} else {gap=0; overlap=0};
				pre1=$1; pre11=$11;
			}' split.mid.expanded | sort --parallel=$thread -u -S $memG | join -v 1 _sa.fu02f - >> _sa.fu01
		        mv _sa.fu01 _sa.fu0
                        rm _sa.fu02 _sa.fu02f _filter2
		else mv _sa.fu01 _sa.fu0
		fi

##==== Filter 2: FusionMinStartSite
	## breakpoint ($22 now, later $23) and separate start site (chr+pos)
	# sort by breakpoint and  start.site.umi
	## breakpoint stats: num_unique_molecule (numi), num_start_site (nss), num_start_site2 (diff by at least 2, nss2), average MQ (avgMQ)
	sed 's/:umi:/\t/' _sa.fu0 | tr ' ' '\t' | $awk '{OFS="\t"; print $23,$2,$0}' | sed -e 's/C\([^\t]\+\)P\([0-9]\+\)-/\1\t\2\t/' -e 's:/[12]::' | sort --parallel=$thread -k1,1b -k28,28n -k6,6b -S $memG | $awk '{OFS="\t";
		if ($1 == pre1 && $2 == pre2 && $28 == pre28){
			i += 1;
			diff = $3-pre3;
			if (diff < 750000){
				siteID = preSiteID
				if (diff==0){
					if ($4 != pre4){
						numi += 1
					}
				} else {
					numi += 1;
					nss += 1;
					if (diff >1){nss2 += 1}
				}
			} else {
				siteID = NR
				numi=1; nss=1; nss2=1
			};
		} else {
			i=1; numi=1; nss=1; nss2=1; siteID=NR
		};

		print siteID,numi,nss,nss2,$0;
     		pre1=$1; pre2=$2; pre3=$3; pre4=$4; pre28=$28; preSiteID=siteID
	}' > breakpoint.stats

	##====  Apply Filter2, and
	## min start site step size of 2 (Deprecated: set to 1)
		minStartStepSize=1
	cut -f1-4 breakpoint.stats | tac | sort --parallel=$thread -k1,1b -u -S $memG | $awk -v FusionMinStartSite=$FusionMinStartSite -v minStartStepSize=$minStartStepSize \
		'{OFS="\t"; if ($3 >= FusionMinStartSite && $4 >= minStartStepSize){
			print $1,$2,$3,$4
			}
	}' | join -1 1 -2 1 - <(sort --parallel=$thread -k1,1b -S $memG breakpoint.stats) | sed 's/ /:umi:/12' \
                | tr ' ' '\t' | cut -f 2,3,4,12- > breakpoint.candidates

rm _sa.fu0 

## End: breakpoint candidates
