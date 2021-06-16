#!/bin/bash

. config.txt
memG=$(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)))G
awk=$(which mawk)
if [ "$awk" == "" ]; then
  awk=$(which gawk)
  if [ "$awk" == "" ]; then
    echo "Neither mawk nor gawk is installed in the system!"
    exit
  fi
fi

##==== Annotate breakpoint gene, exon, cDNA position

## Adjust breakpoint for annotation
#	- introduce overlap-removed point (orp)
#	- to correctly annotate exon and cDNA position (some alignment ends in intron could be overlapped with,
#		 and should belong to, next (on the sense strand if on cDNA) split alignment)
#	- Add 6 bases to breakpoint to make orp, which needs to be accounted for in later frame determination
## Make mock variants for annotation to get gene, exon and fucntion info
## Left fileds 1-5, Right fileds 6-10, readID 11
#	    print $6"_"orpL,$6,orpL,$8,$10,$16"_"orpR,$16,orpR,$17,$20,$4,$1,$2,$3,$25,$26
    $awk '{OFS="\t"; if ($10=="+") {
		    orpL = $8-6
	    } else if ($10=="-") {
		    orpL = $8+6
	    };

	    if ($20=="+") {
		    orpR = $17+6
	    } else if ($20=="-") {
		    orpR = $17-6
	    };

	    print $6"_"orpL,$6,orpL,$8,$10,$4,$1,$2,$3,$25,$26 > "__orpLeft";
	    print $16"_"orpR,$16,orpR,$17,$20,$4 > "__orpRight";
	    printf "%s\t%d\n%s\t%d\n",$6,orpL,$16,orpR;
    }' breakpoint.candidates | sort --parallel=$thread -u -S $memG > __breakpoint.for.anno0

## Annotate

    if [ $AnnotationMethod = "annovar" ]; then
	$awk '{print $1,$2,$2,"A","A"}' __breakpoint.for.anno0 > __breakpoint.for.anno
	$perl $annovar/table_annovar.pl __breakpoint.for.anno $annovar/humandb/ -buildver $genomeVer -out __breakpoint.annotated -remove -protocol refGene -operation g -nastring NA > /dev/null 2>&1
	$R -e "source(\"${SFpath}/R/annovar.exon.cds.extraction.R\");annovar.exon.cds.extraction(input = \"__breakpoint.annotated.${genomeVer}_multianno.txt\")" > /dev/null 2>&1
    fi


    #if [ $AnnotationMethod = "snpEff" ]; then
	# under dev...
    #fi

sort --parallel=$thread -k1,1b -S $memG __breakpoint.annotated.${genomeVer}_multianno.txt.ext0 > __breakpoint.annotated.extr

## Merge annotation back to read
	sort --parallel=$thread -k1,1b -S $memG __orpLeft | join - __breakpoint.annotated.extr | sort --parallel=$thread -k6,6b -S $memG > _anno.left
	sort --parallel=$thread -k1,1b -S $memG __orpRight | join - __breakpoint.annotated.extr | sort --parallel=$thread -k6,6b -S $memG | join -1 6 -2 6 _anno.left - > anno.left.right

##=== anno middle split
	if [ -s split.mid ]; then
	    $awk '{mid = $4 + 10; print $3,mid,mid,"A","A",$1,$3,$4,$5}' split.mid > _mid.for.anno0

		if [ $AnnotationMethod = "annovar" ]; then
		    tr ' ' '\t' < _mid.for.anno0 | cut -f1-5 | sort --parallel=$thread -u > _mid.for.anno
		    $perl $annovar/table_annovar.pl _mid.for.anno $annovar/humandb/ -buildver $genomeVer -out _mid.anno -remove -protocol refGene -operation g -nastring NA > /dev/null 2>&1
		    $R -e "source(\"${SFpath}/R/annovar.exon.cds.extraction.R\");annovar.exon.cds.extraction(input = \"_mid.anno.${genomeVer}_multianno.txt\")" > /dev/null 2>&1
		fi

	    #if [ $AnnotationMethod = "snpEff" ]; then
		# under dev...
	    #fi

		tr ' ' '\t' < _mid.for.anno0 | sed 's:\t:_:' | sort --parallel=$thread -k1,1b > _mid.for.anno1
		sort --parallel=$thread -k1,1b _mid.anno.${genomeVer}_multianno.txt.ext0 | join _mid.for.anno1 - | cut -d ' ' -f2,5- > anno.mid
		rm _mid.for.anno0 _mid.for.anno _mid.anno.${genomeVer}_multianno.txt
	fi
rm __orpLeft __orpRight __breakpoint.for.anno __breakpoint.for.anno0 _anno.left __breakpoint.annotated.extr __breakpoint.annotated.${genomeVer}_multianno.txt __breakpoint.annotated.${genomeVer}_multianno.txt.ext0
