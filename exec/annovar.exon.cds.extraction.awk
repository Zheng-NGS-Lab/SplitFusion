#!/usr/bin/mawk
BEGIN {
	while (getline<sprintf("%s/inst/data/%s_refGene.mostExon.NM",w,p)) {
		a[$12]=$1;
	}
	FS="\t";
# exonic
# mostExon transcript > other transcript
	while (getline<sprintf("%s.refGene.exonic_variant_function",q)) {
		if ($3=="UNKNOWN") {
			continue;
		}
		gsub(/ /,"_",$2);
		split($4,c," ");
		n=split($3,d,",");
		for(i=1;i<n;++i) {
			split(d[i],e,":");
			if (!($4 in g) || (a[g[$4]]!=nm[$4] && a[e[1]]==e[2])) {
				g[$4]=e[1];
				nm[$4]=e[2];
				ex[$4]=e[3];
				if(e[4]~/A$/){dir[$4]="+"}else{dir[$4]="-"}
				pos[$4]=substr(e[4],4,length(e[4])-4);
			}
		}
		if (!($4 in b) || a[g[$4]]==nm[$4]) {
			b[$4]=sprintf("%s_%s\t%s\t%s\texonic\t%s\t%s\t%s\t%d",c[1],c[2],g[$4],dir[$4],$2,nm[$4],ex[$4],pos[$4]);
			bs[$4]=0;
		}
	}
# 0. exonic
# 1. UTR5
# 2. UTR3
# 3. splicing
# 4. ncRNA_exonic
# 5. ncRNA_splicing
# 6. upstream
# 7. downstream
# 8. intronic
# 9. ncRNA_intronic
# 10. intergenic
	nc["exonic"]=0;
	nc["UTR5"]=1;nc["UTR3"]=2;nc["splicing"]=3;nc["intronic"]=8;
	nc["upstream"]=6;nc["downstream"]=7;nc["ncRNA_exonic"]=4;
	nc["ncRNA_splicing"]=5;nc["ncRNA_intronic"]=9;nc["intergenic"]=10;
	while (getline<sprintf("%s.refGene.variant_function",q)) {
		if (!($1 in nc)) {
			nc[$1]=11;
		}
		if ($3 in b) {
		       	if (bs[$3]==0) {
				continue;
			} else if (nc[$1] < bs[$3]) {
				split($3,c," ");
				if (nc[$1]==0) {
					b[$3]=sprintf("%s_%s\t%s\tNA\texonic\tunknown\tUNKNOWN\tUNKNOWN\tUNKNOWN",c[1],c[2],$2);
				} else {
					split($2,d,"(");
				        gene = d[1];
				        dir[$3]="NA";
			        	if ($1~/UTR/ || $1~/splicing/) {
						n=split($2,e,"),");
						for(i=1;i<=n;i++) {
							split(e[i],f,"(");
							if (f[1] in a) {
								if (f[2]~/A>A/) {
									dir[$3]="+";
								} else{
									dir[$3]="-";
								}
								gene = f[1]
							}
						}
					}
					if (gene in a) {
						b[$3]=sprintf("%s_%s\t%s\t%s\t%s\tNA\t%s\t%s\tNA",c[1],c[2],gene,dir[$3],$1,a[d[1]],$1);
					} else {
						b[$3]=sprintf("%s_%s\t%s\tNA\t%s\tNA\tNA\t%s\tNA",c[1],c[2],gene,$1,$1);
					}
				}
				bs[$3]=nc[$1];
			}
		} else {
			split($3,c," ");
			if (nc[$1]==0) {
				b[$3]=sprintf("%s_%s\t%s\tNA\texonic\tunknown\tUNKNOWN\tUNKNOWN\tUNKNOWN",c[1],c[2],$2);
			} else {
				split($2,d,"(");
				gene = d[1];
				dir[$3]="NA";
				if ($1~/UTR/ || $1~/splicing/) {
					n=split($2,e,"),");
					for(i=1;i<=n;i++) {
						split(e[i],f,"(");
						if (f[1] in a) {
							if (f[2]~/A>A/) {
								dir[$3]="+";
							} else{
								dir[$3]="-";
							}
							gene = f[1]
						}
					}
				}
				if (gene in a) {
					b[$3]=sprintf("%s_%s\t%s\t%s\t%s\tNA\t%s\t%s\tNA",c[1],c[2],gene,dir[$3],$1,a[d[1]],$1);
				} else {
					b[$3]=sprintf("%s_%s\t%s\tNA\t%s\tNA\tNA\t%s\tNA",c[1],c[2],gene,$1,$1);
				}
			}
			bs[$3]=nc[$1];
		}
	}
}
{
	FS=" ";
 	if ($0 in b) {
		printf "%s\n",b[$0];
	} else {
		printf "%s_%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n",$1,$2;
	}
}
