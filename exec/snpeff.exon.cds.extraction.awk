#!/usr/bin/mawk
BEGIN {
        while (getline<sprintf("%s/inst/data/%s_refGene.mostExon.NM",w,p)) {
                v[$12]=$1;
        }
}
!/^#/ {
	g="NA";t="NA";dir="NA";nm="NA";ex="NA";pos="NA";
	n=split($8,a,"|");
	for (i=0;i<n;i+=15) {
		gsub(/\..*$/,"",a[i+7]);
		if (a[i+7]==v[a[i+4]]) {break}
	}
	if (i+1>=n) {i=0}
#	if (a[2]~/^intergenic/){printf "XXX %s\t%d\t%d\n",$0,i,n}
	if (a[i+2]~/^frameshift/) {
		g=a[i+4];nm=a[i+7];t="exonic";
		g=a[i+4];nm=a[i+7];
		if(a[i+10]~/T$/){dir="-"}else{dir="+"}
		split(a[i+9],b,"/");
		ex=sprintf("exon%d",b[1]);
		split(a[i+13],c,"/");
		pos=c[1];
	} else if (a[i+2]~/_prime_UTR/) {
		if (a[i+2]~/^3/) {t="UTR3";} else {t="UTR5"}
		g=a[i+4];nm=a[i+7];ex=t;
		if(a[i+10]~/T$/){dir="-"}else{dir="+"}
	} else if (a[i+2]~/splice_donor/||a[i+2]~/splice_acceptor/) {
		g=a[i+4];nm=a[i+7];t="splicing";
		if(a[i+10]~/T$/){dir="-"}else{dir="+"}
	} else if (a[i+2]~/^intergenic/) {
		g=a[i+4];t="intergenic";ex=t;
	} else if (a[i+2]~/stream_/ || a[i+2]~/intron/) {
		if (a[i+2]~/^up/) {t="upstream"}
		else if (a[i+2]~/^down/) {t="downstream"}
		else if (a[i+7]~/^NR_/) {t="ncRNA_intronic"}
		else {t="intronic"}
		g=a[i+4];nm=a[i+7];ex=t;
	} else if (a[i+2]~/non_coding_transcript_exon/) {
		g=a[i+4];nm=a[i+7];t="ncRNA_exonic";ex=t;
	}
	OFS="\t";
	printf "%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,g,dir,t,"NA",nm,ex,pos;
}
