BEGIN {
  p=0;
  while(getline<sprintf("%s/inst/data/%s_pseudogene.bed",w,v)) {
    if (!/^#/) {
      c[$4]=$1;
      s[$4]=$2;
      e[$4]=$3;
      a[++p]=$4;
    }
  }
}
{
  if (!/^@/&&$5==0&&/XA:Z:/) {
# First, search for XA:Z field
    for(xaz=1;xaz<=NF;++xaz){if($xaz~/^XA:Z:/){break}}
# first check if the primary mapping is already a representative of 
# a pseduogene family.
    found=0; 
    for(j=1;j<=p;++j) {
      i=a[j];
      if(c[i]==$3&&$4>s[i]&&$4<=e[i]) {
	found=j;
        break;
      }
    }
# if not, check XA:Z tag to see if any alternative alignment matches
# if so, loop over the list of alternative mappings
# if none of the alternative mapping on our pseduogene.bed, just
# print the original line
    if(found==0) {
      split($xaz,f,":");
      am=split(f[3],g,";");
      for(k=1;k<am;++k) {
	split(g[k],h,",");
	dir=substr(h[2],1,1);
	pos=substr(h[2],2);
        for(j=1;j<=p;++j) {
          i=a[j];
          if(c[i]==h[1]&&pos>s[i]&&pos<=e[i]) {
	    found=j;
	# if different direction
#            if(and($2,16)==0&&dir=="-"||and($2,16)!=0&&dir=="+") {
#	      $3=h[1];
#	      $4=pos;
#	      if($6~/H/){
#		sub(/S/,"H",h[3]);
#	      }
#              $6=h[3];
#	      if(and($2,16)==0){$2=$2+16}
#              else{$2=$2-16}
#	    } else {$3=h[1];$4=pos;}
	    $3=h[1];$4=pos;
	    break;
          }
	}
      }
    }
    if(found!=0) {
# remove XA:Z, print
      printf "%s",$1;
      for(j=2;j<=NF;++j) {
	if(j!=xaz){printf "\t%s",$j}
      }
      printf "\n";
    }else{print}
  }
  else {print}
}
