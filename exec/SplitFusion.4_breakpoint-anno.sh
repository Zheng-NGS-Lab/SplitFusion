#!/bin/bash

. config.txt

##==== Annotate breakpoint gene, exon, cDNA position

## Adjust breakpoint for annotation
#	- introduce overlap-removed point (orp)
#	- to correctly annotate exon and cDNA position (some alignment ends in intron could be overlapped with,
#		 and should belong to, next (on the sense strand if on cDNA) split alignment)
#	- Add 6 bases to breakpoint to make orp, which needs to be accounted for in later frame determination
    awk '{OFS="\t"; if ($10=="+") {
		    orpL = $8-6
	    } else if ($10=="-") {
		    orpL = $8+6
	    };

	    if ($20=="+") {
		    orpR = $17+6
	    } else if ($20=="-") {
		    orpR = $17-6
	    };

	    print $6"_"orpL,$6,orpL,$8,$10,$16"_"orpR,$16,orpR,$17,$20,$4,$1,$2,$3,$25,$26
    }' breakpoint.candidates > _orp

## Make mock variants for annotation to get gene, exon and fucntion info
## Left fileds 1-5, Right fileds 6-10, readID 11
    cut -f 1-5,11- _orp > __orpLeft 
    cut -f 6-11 _orp > __orpRight
    
    cat __orpLeft __orpRight | cut -f 2,3 | sort --parallel=$thread -u > __breakpoint.for.anno0

## Annotate

    if [ $AnnotationMethod = "annovar" ]; then
	awk '{print $1,$2,$2,"A","A"}' __breakpoint.for.anno0 > __breakpoint.for.anno
	$perl $annovar/table_annovar.pl __breakpoint.for.anno $annovar/humandb/ -buildver $genomeVer -out __breakpoint.annotated -remove -protocol refGene -operation g -nastring NA > /dev/null 2>&1
	$R -e "library(SplitFusion);annovar.exon.cds.extraction(input = \"__breakpoint.annotated.${genomeVer}_multianno.txt\")" > /dev/null 2>&1
    fi


    #if [ $AnnotationMethod = "snpEff" ]; then
	# under dev...
    #fi

sort --parallel=$thread -k1,1b __breakpoint.annotated.${genomeVer}_multianno.txt.ext0 > __breakpoint.annotated.extr

## Merge annotation back to read
	sort --parallel=$thread -k1,1b __orpLeft > __orpLeft.s
	sort --parallel=$thread -k1,1b __orpRight > __orpRight.s
	join __orpLeft.s  __breakpoint.annotated.extr > __anno.left
	join __orpRight.s __breakpoint.annotated.extr > __anno.right

	sort --parallel=$thread -k6,6b __anno.left > _anno.left
	sort --parallel=$thread -k6,6b __anno.right > _anno.right
	join -1 6 -2 6 _anno.left _anno.right > anno.left.right

##=== anno middle split
	if [ -s split.mid ]; then
	    awk '{mid = $4 + 10; print $3,mid,mid,"A","A",$1,$3,$4,$5}' split.mid > _mid.for.anno0

		if [ $AnnotationMethod = "annovar" ]; then
		    tr ' ' '\t' < _mid.for.anno0 | cut -f1-5 | sort --parallel=$thread -u > _mid.for.anno
		    $perl $annovar/table_annovar.pl _mid.for.anno $annovar/humandb/ -buildver $genomeVer -out _mid.anno -remove -protocol refGene -operation g -nastring NA > /dev/null 2>&1
		    $R -e "library(SplitFusion);annovar.exon.cds.extraction(input = \"_mid.anno.${genomeVer}_multianno.txt\")" > /dev/null 2>&1
		fi

	    #if [ $AnnotationMethod = "snpEff" ]; then
		# under dev...
	    #fi

		tr ' ' '\t' < _mid.for.anno0 | sed 's:\t:_:' | sort --parallel=$thread -k1,1b > _mid.for.anno1
		sort --parallel=$thread -k1,1b _mid.anno.${genomeVer}_multianno.txt.ext0 > _mid.anno.ext
		join _mid.for.anno1 _mid.anno.ext | cut -d ' ' -f2,5- > anno.mid
		
	fi
rm _*
