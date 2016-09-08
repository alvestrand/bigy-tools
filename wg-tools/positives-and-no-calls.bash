#!/bin/bash
#
# Analyze the .bed/.vcf files for positives and no-calls.
#

cp variant-list.txt variant-match.txt
for BEDFILE in ${@}; do
	VCFFILE=`echo "$BEDFILE" | sed 's/.bed/.vcf/'`
	gawk 'NR==FNR {d[NR]=$0;split($0,x,",");s[NR]=x[1];n++} \
	      NR!=FNR && $1=="chrY" && $7=="PASS" && $4!="." && $5!="." {for (i=1;i<=n;i++) if ($2==s[i]) p[i]=$2"."$4"."$5} \
		  END {for (i=1;i<=n;i++) print d[i]","p[i]}' variant-match.txt "$VCFFILE" > foo
	gawk 'NR==FNR {d[NR]=$0;split($0,x,",");s[NR]=x[1];n++;xi=1} \
	      NR!=FNR && $1=="chrY" {m++;a[m]=$2;b[m]=$3} \
		  END {for (i=1;i<=n;i++) {call=";nc"; for (j=xi;j<=m;j++) {if (a[j]<=s[i] && s[i]<=b[j]) {call=""; \
		          if (a[j]==s[i]) call=";cbl"; if (b[j]==s[i]) call=";cbu"; if (a[j]==s[i] && b[j]==s[i]) call=";cblu"; xi=j; break}; \
				  if (a[j]>s[i]) break}; print d[i]""call}}' foo "$BEDFILE" > variant-match.txt
	echo -n "."
done
echo ""
