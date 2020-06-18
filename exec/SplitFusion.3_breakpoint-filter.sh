#!/bin/bash

. config.txt
SampleId=$( pwd | sed "s:.*/::")

##========================================
## If exists TARGET panel information, use the SplitFusion-TARGET model
##========================================

if [ $panel != "NA" ]; then
	gawk '{if ($1 ~ /\/2/ && $9 == 1) print}' breakpoint.candidates.preFilter > _preFilter.r2
	gawk '{OFS="\t"; if ($7 == "-"){start=$5; $5=$4; $4=start};
		print $3,$4,$5,$1,$6,$7
		}' _preFilter.r2 > _preFilter.r2.bed

		sort --parallel=$thread -k1,1n -k2,2n _preFilter.r2.bed > _preFilter.r2.srt.bed
	$bedtools intersect -wa -wb -a $panel_dir/$panel.GSP2.bed -b _preFilter.r2.srt.bed -S > _preFilter.GSP2.bed
                # anchored
                awk '{diff1 = $8 - $2; diff2 = $9 - $3;
		if ($6 == "-"){
			if (diff1 >= -1 && diff1 <=1) {
					if (diff2 >= 3) {
						print > "_preFilter.anchored.onTarget.bed"
					} else {print > "_preFilter.anchored.offTarget.bed"}
				} else {print > "_preFilter.nonAnchored.bed"}
		} else {
			if (diff2 >= -1 && diff2 <=1) {
					if (diff1 <= -3) {
                        	        print > "_preFilter.anchored.onTarget.bed"
				} else {print > "_preFilter.anchored.offTarget.bed"}
			} else {print > "_preFilter.nonAnchored.bed"}
		}}' _preFilter.GSP2.bed 

	cut -f 10 _preFilter.anchored.onTarget.bed | sed "s:umi.*::" > _preFilter.anchored.onTarget.readID0
        sort --parallel=$thread -k1,1b -u _preFilter.anchored.onTarget.readID0 > _preFilter.anchored.onTarget.readID.srt
	
	sed "s:umi:\tumi:" breakpoint.candidates.preFilter > _breakpoint.candidates.preFilter
        sort --parallel=$thread -k1,1b _breakpoint.candidates.preFilter > _breakpoint.candidates.preFilter.srt

	join _preFilter.anchored.onTarget.readID.srt _breakpoint.candidates.preFilter.srt > _breakpoint.candidates.preFilter.target
	sed 's: ::' _breakpoint.candidates.preFilter.target > breakpoint.candidates.preFilter
	rm _*


	##==== reads with middle split
	if [ -f breakpoint.candidates.preFilter.w.mid ]; then	
	    gawk '{if ($1 ~ /\/2/ && $9 == 1) print}' breakpoint.candidates.preFilter.w.mid > _preFilter.r2
	    gawk '{OFS="\t"; if ($7 == "-"){start=$5; $5=$4; $4=start};
		    print $3,$4,$5,$1,$6,$7
		    }' _preFilter.r2 > _preFilter.r2.bed

		    sort --parallel=$thread -k1,1n -k2,2n _preFilter.r2.bed > _preFilter.r2.srt.bed
	    $bedtools intersect -wa -wb -a $panel_dir/$panel.GSP2.bed -b _preFilter.r2.srt.bed -S > _preFilter.GSP2.bed
		    # anchored
		    awk '{diff1 = $8 - $2; diff2 = $9 - $3;
		    if ($6 == "-"){
			    if (diff1 >= -1 && diff1 <=1) {
					    if (diff2 >= 3) {
						    print > "_preFilter.anchored.onTarget.bed"
					    } else {print > "_preFilter.anchored.offTarget.bed"}
				    } else {print > "_preFilter.nonAnchored.bed"}
		    } else {
			    if (diff2 >= -1 && diff2 <=1) {
					    if (diff1 <= -3) {
					    print > "_preFilter.anchored.onTarget.bed"
				    } else {print > "_preFilter.anchored.offTarget.bed"}
			    } else {print > "_preFilter.nonAnchored.bed"}
		    }}' _preFilter.GSP2.bed 

	    cut -f 10 _preFilter.anchored.onTarget.bed | sed "s:umi.*::" > _preFilter.anchored.onTarget.readID0
	    sort --parallel=$thread -k1,1b -u _preFilter.anchored.onTarget.readID0 > _preFilter.anchored.onTarget.readID.srt
	    
	    sed "s:umi:\tumi:" breakpoint.candidates.preFilter.w.mid > _breakpoint.candidates.preFilter
	    sort --parallel=$thread -k1,1b _breakpoint.candidates.preFilter > _breakpoint.candidates.preFilter.srt

	    join _preFilter.anchored.onTarget.readID.srt _breakpoint.candidates.preFilter.srt > _breakpoint.candidates.preFilter.target
	    sed 's: ::' _breakpoint.candidates.preFilter.target > breakpoint.candidates.preFilter.w.mid
	    rm _*
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
		echo | awk -v minMapLength=$minMapLength -v minMapLength2=$minMapLength2 -v minExclusive=$minExclusive -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap -v minMQ1=$minMQ1 \
		    '{ gap = $19-$10-1; overlap = $10-$19+1;
			if ($1 ~ /\/1/)	{mapLen1 = $10-$9+1; mapLen2 = $20-$19+1; mq1=$6
				} else { mapLen2 = $10-$9+1; mapLen1 = $20-$19+1; mq1=$16};
			if ((mapLen1 >= minMapLength && mapLen2 >= minMapLength2 && mq1 >= minMQ1) \
				&& ($19-$9 >= minExclusive && $20-$10 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			    ) {print $0,overlap} else {print $1 > "_filter1"}
			}' breakpoint.candidates.preFilter > _sa.fu01

		## For reads with middle split, turn off maxQueryGap by let maxQueryGap=1000
		    ## left:  $9-------$10
		    ## right:		       $19------$20
		if [ -f breakpoint.candidates.preFilter.w.mid ]; then
			touch _filter2
		 echo | awk -v minMapLength=$minMapLength -v minMapLength2=$minMapLength2 -v minExclusive=$minExclusive -v maxQueryGap=1000 \
		    '{ gap = $19-$10-1; overlap = $10-$19+1;
			if ($1 ~ /\/1/) {mapLen1 = $10-$9+1; mapLen2 = $20-$19+1
                                } else { mapLen2 = $10-$9+1; mapLen1 = $20-$19+1};

			if (  (mapLen1 >= minMapLength && mapLen2 >= minMapLength2) \
					&& ($19-$9 >= minExclusive && $20-$10 >= minExclusive && gap <= maxQueryGap) \
			    ) {print $0,overlap} else {print $1 > "_filter2"}
			}' breakpoint.candidates.preFilter.w.mid > _sa.fu02

			join -v 1 _sa.fu02 _filter2 > _sa.fu02f

			# filter gap and overlap
			echo | awk -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap \
			'{ if ($1==pre1){
					gap = $10 - pre11 -1; overlap = pre11 - $10 + 1
					if (gap > maxQueryGap || overlap > maxOverlap) {print $1}
				} else {gap=0; overlap=0};
				pre1=$1; pre11=$11;
			}' split.mid.expanded > _mid.gap.id
			sort --parallel=$thread -u _mid.gap.id > _mid.gap.id.u

			join -v 1 _sa.fu02f _mid.gap.id.u > _sa.fu02ff 
			cat _sa.fu01 _sa.fu02ff > _sa.fu0
		else cp _sa.fu01 _sa.fu0
		fi

##==== Filter 2: FusionMinStartSite
	## breakpoint ($22 now, later $23) and separate start site (chr+pos)
	sed 's/:umi:/\t/' _sa.fu0 | tr ' ' '\t' | awk '{OFS="\t"; print $23,$2,$0}' | sed -e 's/C\([^\t]\+\)P\([0-9]\+\)-/\1\t\2\t/' -e 's:/[12]::' > _sa.fu2

	# sort by breakpoint and  start.site.umi
	sort --parallel=$thread -k1,1b -k28,28n -k6,6b _sa.fu2 > _sa.fu3

	## breakpoint stats: num_unique_molecule (numi), num_start_site (nss), num_start_site2 (diff by at least 2, nss2), average MQ (avgMQ)
        awk '{OFS="\t";
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

		print siteID,numi,nss,nss2,$0 > "breakpoint.stats";
     		pre1=$1; pre2=$2; pre3=$3; pre4=$4; pre28=$28; preSiteID=siteID
	}' _sa.fu3

	##====  Apply Filter2, and
	## min start site step size of 2 (Deprecated: set to 1)
		minStartStepSize=1
	cut -f1-4 breakpoint.stats | tac > _breakpoint.stats.tac
	sort --parallel=$thread -k1,1b -u _breakpoint.stats.tac > _breakpoint.stats.u
	echo | awk -v FusionMinStartSite=$FusionMinStartSite -v minStartStepSize=$minStartStepSize \
		'{OFS="\t"; if ($3 >= FusionMinStartSite && $4 >= minStartStepSize){
			print $1,$2,$3,$4 > "breakpoint.siteID.stats.filtered"
			}
		}' _breakpoint.stats.u

	sort --parallel=$thread -k1,1b breakpoint.stats > _breakpoint.stats.s
	join -1 1 -2 1 breakpoint.siteID.stats.filtered  _breakpoint.stats.s \
		| sed 's/ /:umi:/12' \
		| tr ' ' '\t' | cut -f 2,3,4,12- > breakpoint.candidates

touch _0; rm _*

## End: breakpoint candidates
