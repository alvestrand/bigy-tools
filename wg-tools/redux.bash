# !/bin/bash
# Author: Iain McDonald
# Purpose: Reduction and comparison script for Family Tree DNA's BigY data
# Release: 0.2.2 2016/09/05
# For free distribution under the terms of the GNU General Public License, version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# Set up instructions:
#    Create sub-directory called zip, and insert BigY zip files into the folder with the format bigy-<NAME>-<KIT NUMBER>.zip
#    Create lists of implications, bad SNPs, recurrent SNPs and SNP names in the following files:
#		implications.txt badlist.txt recurrencies.txt
#    with the following formats:
#         7246726 > 23612197 : L48 > Z381
#    	  22436300 : E206
#         20723243 : F552
#	The first example implies that anyone who is L48+ is also Z381+. Only columns 1 and 3 are used.
#   The second example states that E206 is an inconsistent SNP. Only the first column is used.
#   The third example shows the (same) format for recurrent SNPs, indicating F552 is recurrent. Only the first column is used.
#
# Windows WSL BASH users:
#   Non-native package dependences:
#      unzip
#   These must be installed before the script is run, otherwise an error message will be generated. Try:
#      sudo apt-get install unzip

# For circa 700 kits, this process takes around 2.5 hours on a heavy-duty 2014 laptop.

# Change log:
# 0.2.2.20160905a - Fixed bug in shared SNP counting involving presumed positives being counted twice
# 0.2.1.20160901a - BETA RELEASE (716): Introduced short report format
#					Introduced basic age calculation (based solely on SNP counts)
#					Fixed bug involving the NFILES parameter not being set on an initial run
#                   Moved backup location to be more intuitive
# 0.2.0.20160831a - Private release to Harald A.
#					Introduced SNP counts
#                   Encoded MAKE* flags to begin age analysis
#                   Encoded clade identification, tree formation and basic SNP counting
#					Fixed bug for multiple SNP names
#					Discontinued use of forced SNP names in favour of the single YBrowse output file
# 0.1.1.20160822a - BETA RELEASE (710): allowed basic compilation of report.csv

# Wish list:
# Basic report: parallisation / further optimisation
# Age analysis: the ability to select allowed regions of the chromosome where SNPs should be counted for age analysis.
# Age analysis: implementation of causality
# Age analysis: uncertainty ranges based on Poisson statistics
# HTML report: in the style of Alex W's BigTree
# General: BigY+YElite+WGS crossover
# Basic report: sorting criterion to match FTDNA tree

# These flags define what actually gets done
# Unless MAKEREPORT is non-zero, the other reports will be based on the current report.csv
# MAKEREPORT = do the main reduction and create report.csv, else use an existing one
# MAKEAGES = perform age analysis
# MAKESHORTREPORT = make a shorter version of the report, collapsing the singleton section and removing the shared and inconsistent list of SNPs
# MAKEHTMLTREE = make an HTML version of the tree for easier visualisation
MAKEREPORT=1
MAKEAGES=1
MAKESHORTREPORT=1
MAKEHTMLREPORT=1

# Set SKIPZIP to >0 skip the labour-intensive parts of the script. Requires MAKEREPORT=TRUE
# This should only be done if you have run the script before on the same dataset!
# If set to 1, unzipping and VCF/BED file scanning will be ignored
#    Statistics generation, sorting of already-identified variants into good/bad/shared lists and clade sorting will still take place
#    This reduces the processing time by a factor of ~12x (e.g. ~2 hours to ~10 minutes) but requires a successful (but not necessarily sorted) previous run.
# If set to 2, statistics generation is skipped. The statistics INCLUDING KIT NAMES will be taken from the header of the previous run.
#    It also uses the ORDER of the previous run, taken from order.txt
#    This reduces the processing time by a further factor of ~2x (e.g. to ~5 minutes) but requires a complete successful previous run
# If set to 3, SNP names aren't updated either. The previous SNP names are used.
SKIPZIP="0"

# The following parameters set the mutation rate for the age analysis.
# The rate is set in SNP mutations per year per *billion* base pairs.
# RATE0 and RATE1 give the lower and upper bounds for the 95% confidence interval.
RATE=0.822
RATE0=0.762
RATE1=0.880

# That should be everything you need to set up.
# Let's do stuff.

T0=`date +%s.%N`

# ---------------------------------------
# Always make a backup

echo "CREATING BACKUP COPIES OF EXISTING FILES..."

if [ ! -d "autobackup" ]; then
mkdir autobackup
else
rm -rf autobackup/*
fi
cp variant-*.txt autobackup/
cp snps_hg19.csv autobackup/
cp snp-names.csv autobackup/
cp snp-used.csv autobackup/
cp report.csv autobackup/
cp short-report.csv autobackup/
cp clades.csv autobackup/
cp tree.txt autobackup/
cp raw-ages.txt autobackup/

if [ "$MAKEREPORT" -gt "0" ]; then

echo "MAKING REPORT..."
rm -f report.csv

if [ "$SKIPZIP" == "0" ]; then

# Check the folder containing the zip files exists
if [ ! -d "zip" ]; then
echo "Input zip folder does not appear to exist. Aborting."
exit 1
fi

# Check the zip working exists
# Empty it, otherwise make it
if [ ! -d "working" ]; then
mkdir working
else
rm -rf working/*
fi

# Check the folder for the output exists
# Empty it, otherwise make it
if [ ! -d "unzip" ]; then
mkdir unzip
else
rm -f unzip/*.bed unzip/*.vcf
fi

# Get the list of input files
FILES=(`ls zip/big*.zip`)

if [ ${#FILES[@]} == 0 ]; then
echo "No input files detected in zip folder. Aborting."
echo "Check name format: should be bigy-<NAME>-<NUMBER>.zip"
exit 1
else
echo ${#FILES[@]} "input files detected"
fi

# Check whether unzip is installed
command -v unzip >/dev/null 2>&1 || { echo >&2 "Unzip package not found. Aborting."; exit 1; }

# Unzip each folder in turn
echo "Unzipping..."
FILECOUNT=0
for ZIPFILE in ${FILES[@]}; do
	let FILECOUNT+=1
	PREFIX=`echo "$ZIPFILE" | gawk -F- '{print $2"-"$3}' | sed 's/.zip//'`
	#echo $FILECOUNT: $ZIPFILE : $PREFIX
	unzip -q $ZIPFILE -d working/
	if [ -s working/*.vcf ]; then mv working/*.vcf working/"$PREFIX".vcf; fi
	if [ -s working/*.bed ]; then mv working/*.bed working/"$PREFIX".bed; fi
	if [ -s working/*/variants.vcf ]; then mv working/*/variants.vcf working/"$PREFIX".vcf; fi
	if [ -s working/*/regions.bed ]; then mv working/*/regions.bed working/"$PREFIX".bed; fi
	if [ -s working/"$PREFIX".vcf ] && [ -s working/"$PREFIX".bed ]; then
		mv working/"$PREFIX".vcf unzip/;
		mv working/"$PREFIX".bed unzip/;
		else echo ""; echo "Warning: could not identify VCF and/or BED file for $PREFIX"; fi
	rm -r working; mkdir working
    echo -n "."
done
echo ""

# Close SKIPZIP if
fi

# Skip some more if SKIPZIP set
if [ "$SKIPZIP" -gt "1" ]; then
cp header.csv report.csv
NFILES=`head -1 header.csv | gawk -v FS=, '{print NF-17}'`
echo "... $NFILES results to be post-processed"
fi
if [ "$SKIPZIP" -le "1" ]; then
# Check number of BED = number of VCF files
if [ `ls unzip/*.bed | wc -l` != `ls unzip/*.vcf | wc -l` ]; then
echo "Number of BED files does not equal number of VCF files."
echo "This is an unexpected error. Aborting."
exit 1
fi

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Generate statistics from BED & VCF files
echo "Generating preliminary statistics..."
FILES=(`ls unzip/*.bed`)
echo "Total Kits:,,,,,"${#FILES[@]}',,,,,,,,,,,Kit' >> report.csv
echo ',,,,,,,,,,,,,,,,Coverage' >> report.csv
echo ',,,,,,,,,,,,,,,,Regions' >> report.csv
echo ',,,,,,,,,,,,,,,,Variants' >> report.csv
echo ',,,,,,,,,,,,,,,,Passed' >> report.csv
echo ',,,,,,,,,,,,,,,,Simple SNPs' >> report.csv
echo ',,,,,,,,,,,,,,,,SNPs under U106' >> report.csv
echo ',,,,,,,,,,,,,,,,Singleton SNPs' >> report.csv
echo ',,,,,,,,,,,,,,,,Indels' >> report.csv
echo ',,,,,,,,,,,,,,,,Indels under U106' >> report.csv
echo ',,,,,,,,,,,,,,,,Singleton Indels' >> report.csv
echo 'Non-shared SNPs' >> report.csv
echo 'GrCh37,Name(s),Ref,Alt,Type,N+,(?+),N-,nc,cbl+,cbl-,cbu+,cbu-,cblu+,cblu-,1stCol,Recur' >> report.csv
echo "Generating statistics for" ${#FILES[@]} "BED files..."
for BEDFILE in ${FILES[@]}; do
	VCFFILE=`echo "$BEDFILE" | sed 's/.bed/.vcf/'`
	KITNAME=`echo "$BEDFILE" | gawk -F/ '{print $2}' | sed 's/.bed//'`
	STATS=`gawk '{s+=$3-$2-1} END {print k,s,NR}' k="$KITNAME" "$BEDFILE"`
	STATS2=`gawk '$1=="chrY" {n++} $1=="chrY" && $7=="PASS" {v++; if ($4!="." && $5!=".") {if (length($4)==1 && length($5)==1) {s++} else {i++}}} END {print n,v,s,0,0,i,0,0}' "$VCFFILE"`
	STATS="$STATS $STATS2"
	gawk -v s="$STATS" 'NR==1 {split(s,stat," ")} {print $0","stat[NR]}' report.csv > foo
	mv foo report.csv
	echo -n "."
done
echo ""
cp report.csv header.csv

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Close SKIPZIP if
fi

# Skip some more if SKIPZIP set
if [ "$SKIPZIP" == "0" ]; then

# Identify list of variants
echo "Identifying list of variants..."
rm -f variant-list.txt
for BEDFILE in ${FILES[@]}; do
	VCFFILE=`echo "$BEDFILE" | sed 's/.bed/.vcf/'`
	gawk '$1=="chrY" && $7=="PASS" && $4!="." && $5!="." {print $2"\t"$4"\t"$5}' "$VCFFILE" >> variant-list.txt
	echo -n "."
done
echo ""

# Create a unique list of variants
sort -nk1 variant-list.txt | uniq -c | sort -nk2 | gawk '{n="SNP"} length($3)>1 || length($4)>1 {n="Indel"} {print $2",,"$3","$4","n","$1",,,,,,,,,,,"}' > foo; mv foo variant-list.txt

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Write out positive cases
echo "Identifying positives and no calls..."
cp variant-list.txt variant-match.txt
for BEDFILE in ${FILES[@]}; do
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

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Close SKIPZIP if statement
fi
cp variant-match.txt variant-output.txt

# Insert presumed positives
gawk -v FS=, -v OFS=, 'NR==FNR {split($0,u," ");a[NR]=u[1];b[NR]=u[3];n++} \
NR!=FNR {m=NF; for (i=1;i<=m;i++) d[FNR,i]=$i; s[FNR]=$1; l=FNR} \
END {for (i=1;i<=n;i++) {delete t; for (j=1;j<=l;j++) {if (s[j]==a[i]) {for (k=18;k<=m;k++) if (length(d[j,k])>5 && substr(d[j,k],1,1)!=";") {t[k]=1} else {t[k]=0}}}; \
                         for (j=1;j<=l;j++) {if (s[j]==b[i]) {for (k=18;k<=m;k++) if (length(d[j,k])<=5 && t[k]==1 && substr(d[j,k],1,1)==";") d[j,k]="(?+)"d[j,k]}}}; \
     for (i=1;i<=l;i++) {for (j=1;j<=m;j++)  {printf "%s,",d[i,j]}; printf "\n"}}' implications.txt variant-output.txt > foo
mv foo variant-output.txt

# Identify recurrent SNPs
gawk -v FS=, -v OFS=, 'NR==FNR {split($0,u," ");a[NR]=u[1];n++} NR!=FNR {for (j=1;j<=n;j++) if (a[j]==$1) {$2="(R);";for (i=18;i<=NF;i++) if ($i~$1 || $i~/(?+)/) $i="(R);"$i}; print}' recurrencies.txt variant-output.txt > foo
mv foo variant-output.txt
	 
# Generate SNP stats
echo "Generating stats and segregating for SNPs..."
NFILES=`head -1 header.csv | gawk -v FS=, '{print NF-17}'`
gawk -v FS=, -v f="$NFILES" '{cbl=cbu=cblu=presp=nn=nc=cblp=cbln=cbup=cbun=cblup=cblun=0; \
   for (i=18;i<=18+f-1;i++) {if ($i!~$1) nn++; if ($i ~ /(?+)/) presp++; if ($i ~ /nc/) nc++; \
                             if ($i ~ /cbl/) cbl++; if ($i ~ /cbu/) cbu++; if ($i ~ /cblu/) cblu++; if ($i ~ /cbl/ && $i~$1) cblp++; if ($i ~ /cbu/ && $i~$1) cbup++; if ($i ~ /cblu/ && $i~$1) cblup++}; \
							 cbln=cbl-cblp-cblu; cbun=cbu-cbup; cblun=cblu-cblup; $6+=presp; $7=presp+0; $8=nn+0; $9=nc+0; $10=cblp+0; $11=cbln+0; $12=cbup+0; $13=cbun+0; $14=cblup+0; $15=cblun+0; \
							 for (i=1;i<=NF;i++) printf "%s,",$i; printf "\n"}' variant-output.txt > foo
mv foo variant-output.txt

# Remove common SNPs
# SNPs are declared common if the number of (presumed) positives, no calls, and coverage boundaries for that SNP is greater than or equal to the number of input files (i.e. there are no definite negatives)
# The number should never be greater than f, except in circumstances where you delete one of the input files without re-running the entire process again!
gawk -v FS=, -v f="$NFILES" '$6+$9+$11+$13+$15>=f' variant-output.txt > variant-shared.txt
gawk -v FS=, -v f="$NFILES" '$6+$9+$11+$13+$15<f' variant-output.txt > variant-not-shared.txt

# Remove bad SNPs
gawk -v FS=, 'NR==FNR {grch[NR]=$1;n++} NR!=FNR {flag=0; for (i=1;i<=n;i++) if ($1==grch[i]) flag=1; if (flag==1) print}' badlist.txt variant-not-shared.txt > variant-bad.txt
gawk -v FS=, 'NR==FNR {grch[NR]=$1;n++} NR!=FNR {flag=0; for (i=1;i<=n;i++) if ($1==grch[i]) flag=1; if (flag==0) print}' badlist.txt variant-not-shared.txt > foo
mv foo variant-not-shared.txt

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Initial vertical sorting of SNPs
echo "Sorting rows/columns..."
sort -nk6,6r -nk8,8 -nk1,1 -t, variant-not-shared.txt > foo; mv foo variant-not-shared.txt

# Rinse and repeat...
# Increase the number of repetitions (e.g. seq 1 10) if clades fail to sort
REPS=`seq 1 2`
for REP in $REPS; do

echo "Re-sort $REP"

# Sorting SNPs horizontally
ORDER=`cat variant-not-shared.txt | gawk -v FS=, '{for (i=1;i<=NF-17;i++) {a[NR,i]=$(i+17); p[NR,i]=0; if (($(i+17)~$1 || $(i+17)~/(?+)/) && $(i+17)!~/(R)/) p[NR,i]++}} \
           END {s=NR;for (i=1;i<=NF-17;i++) {new[i]=i}; \
		          for (s=NR;s>=1;s--) {n=0; for (i=1;i<=NF-17;i++) {if (p[s,new[i]]==1) {n++;new2[n]=new[i]}}; \
				     for (i=1;i<=NF-17;i++) {if (p[s,new[i]]==0) {n++;new2[n]=new[i]}}; \
					 for (i=1;i<=NF-17;i++) new[i]=new2[i]}; for (i=1;i<=NF-17;i++) printf "%i ",new[i]}'`
gawk -v FS=, -v o="$ORDER" 'NR==1 {n=split(o,s," ")} {for (i=1;i<=17;i++) printf "%s,",$i; for (i=1;i<=n;i++) printf "%s,",$(s[i]+17);printf "\n"}' report.csv > foo; mv foo report.csv
gawk -v FS=, -v o="$ORDER" 'NR==1 {n=split(o,s," ")} {for (i=1;i<=5;i++) printf "%s,",$i; for (i=6;i<=17;i++) if ($i>0) {printf "%04i,",$i} else {printf "%s,",$i}; for (i=1;i<=n;i++) printf "%s,",$(s[i]+17);printf "\n"}' variant-not-shared.txt > foo; mv foo variant-not-shared.txt
gawk -v FS=, -v o="$ORDER" 'NR==1 {n=split(o,s," ")} {for (i=1;i<=5;i++) printf "%s,",$i; for (i=6;i<=17;i++) if ($i>0) {printf "%04i,",$i} else {printf "%s,",$i}; for (i=1;i<=n;i++) printf "%s,",$(s[i]+17);printf "\n"}' variant-shared.txt > foo; mv foo variant-shared.txt
gawk -v FS=, -v o="$ORDER" 'NR==1 {n=split(o,s," ")} {for (i=1;i<=5;i++) printf "%s,",$i; for (i=6;i<=17;i++) if ($i>0) {printf "%04i,",$i} else {printf "%s,",$i}; for (i=1;i<=n;i++) printf "%s,",$(s[i]+17);printf "\n"}' variant-bad.txt > foo; mv foo variant-bad.txt

# Re-sort SNPs vertically
gawk -v FS=, -v OFS=, '{$16=NF; for (i=NF;i>=18;i--) if ($i~$1 || $i~/(?+)/) $16=i; pruns=prun=0; for (i=18;i<=NF;i++) if ($i~$1 || $i ~ /(?+)/) {if (prun==0) {prun=1; pruns++}} else {if (length($i)<2) prun=0}; $17=pruns; print}' variant-not-shared.txt > foo
sort -n -nk6,6r -nk16,16 -nk1,1 -t, foo > variant-not-shared.txt

done

# Sort other files
gawk -v FS=, -v OFS=, '{pruns=prun=0; for (i=18;i<=NF;i++) if ($i~$1 || $i ~ /(?+)/) {if (prun==0) {prun=1; pruns++}} else {if (length($i)<2) prun=0}; $17=pruns; print}' variant-shared.txt > foo; mv foo variant-shared.txt
gawk -v FS=, -v OFS=, '{pruns=prun=0; for (i=18;i<=NF;i++) if ($i~$1 || $i ~ /(?+)/) {if (prun==0) {prun=1; pruns++}} else {if (length($i)<2) prun=0}; $17=pruns; print}' variant-bad.txt > foo
sort -nk17,17 -nk6,6r foo > variant-bad.txt

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"


# Replace SNP names
# Names are taken from snp-mapping.txt for preference, but snps_hg19.csv is used to fill in the blanks
# The latter can be updated from: http://ybrowse.org/gbrowse2/gff/
# The previous set of names is used if SKIPZIP is set to above 2
echo "Introducing SNP names..."
if [ "$SKIPZIP" -le "2" ]; then
	gawk -v FS=, -v OFS=, 'NR>1 {print $4"."$11"."$12,$9}' snps_hg19.csv | sed 's/"//g' | sort -nk1,1 -t. | gawk -v FS=, -v OFS=, '$1==a {b=b"/"$2; print "...",$0} $1!=a {print a,b;a=$1;b=$2} END {print}' > snp-names.csv
	#cat foo | gawk -v FS=. '{print $1,$0}' | sort -nrk1 | uniq -w12 | sort -nk1 | gawk '{print $2}' | uniq > snp-names.csv
	# snp-mapping.txt allows for the inclusion of SNPs not named in the snps_hg19.csv file
	#cat snp-mapping.txt foo | gawk -v FS=. '{print $1,$0}' | sort -nrk1 | uniq -w12 | sort -nk1 | gawk '{print $2}' | uniq > snp-names.csv
	gawk -v FS=, '{print $1,$1"."$3"."$4}' variant-output.txt | sort -nk1 | gawk '{print $2}' > foo
	gawk -v FS=, 'NR==FNR {a[FNR]=$1;na++} NR!=FNR {b[FNR]=$1;bb[FNR]=$0;nb++} END {ij=1; for (i=1;i<=na;i++) {for (j=ij;j<=nb;j++) {if (a[i]==b[j]) {print bb[j]; ij=j; break}; xa=substr(a[i],1,8); xb=substr(b[i],1,8); if (xb+0>xa+0) {break}}}}' foo snp-names.csv > snp-used.csv
fi
gawk -v FS=, -v OFS=, 'NR==FNR {s[NR,1]=$1;s[NR,2]=$2;split($2,arr," ");s[NR,3]=arr[1];n=NR} NR!=FNR {m=$1"."$3"."$4; for (i=1;i<=n;i++) if (s[i,1]==m) {$2=$2""s[i,2]; gsub(m,s[i,3],$0)}; print}' snp-used.csv variant-not-shared.txt > foo; mv foo variant-not-shared.txt
gawk -v FS=, -v OFS=, 'NR==FNR {s[NR,1]=$1;s[NR,2]=$2;split($2,arr," ");s[NR,3]=arr[1];n=NR} NR!=FNR {m=$1"."$3"."$4; for (i=1;i<=n;i++) if (s[i,1]==m) {$2=$2""s[i,2]; gsub(m,s[i,3],$0)}; print}' snp-used.csv variant-shared.txt > foo; mv foo variant-shared.txt
gawk -v FS=, -v OFS=, 'NR==FNR {s[NR,1]=$1;s[NR,2]=$2;split($2,arr," ");s[NR,3]=arr[1];n=NR} NR!=FNR {m=$1"."$3"."$4; for (i=1;i<=n;i++) if (s[i,1]==m) {$2=$2""s[i,2]; gsub(m,s[i,3],$0)}; print}' snp-used.csv variant-bad.txt > foo; mv foo variant-bad.txt

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Counting SNPs
echo "Creating SNP counts..."
U106SNPS=`gawk -v FS=, -v OFS=, '$5=="SNP" {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
U106INDELS=`gawk -v FS=, -v OFS=, '$5=="Indel" {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
SINGSNPS=`gawk -v FS=, -v OFS=, '$5=="SNP" && $6==1 {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
SINGINDELS=`gawk -v FS=, -v OFS=, '$5=="Indel" && $6==1 {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
gawk -v FS=, -v OFS=, -v us="$U106SNPS" -v ui="$U106INDELS" -v ss="$SINGSNPS" -v si="$SINGINDELS" '\
	NR==7 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,us} \
	NR==8 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,ss} \
	NR==10 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,ui} \
	NR==11 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,si} \
	NR<7 || NR>11 {print}' report.csv > foo; mv foo report.csv

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Collating file
echo "Collating file..."
cat variant-not-shared.txt >> report.csv
echo '' >> report.csv
head -1 report.csv | awk '{$6="";$1="Shared SNPs"; print}' > foo; cat foo >> report.csv
echo 'GrCh37,Name(s),Ref,Alt,Type,N+,(?+),N-,nc,cbl+,cbl-,cbu+,cbu-,cblu+,cblu-,1stCol,Recur' >> report.csv
cat variant-shared.txt >> report.csv
echo '' >> report.csv
head -1 report.csv | awk '{$6="";$1="Inconsistent SNPs"; print}' > foo; cat foo >> report.csv
echo 'GrCh37,Name(s),Ref,Alt,Type,N+,(?+),N-,nc,cbl+,cbl-,cbu+,cbu-,cblu+,cblu-,1stCol,Recur' >> report.csv
cat variant-bad.txt >> report.csv

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Make TSV copy from CSV and backup files
# Strictly speaking, this isn't necessary, as the CSV output satisfies most demands
# However, it acts as a useful backup for those "oh ****" moments and is needed if SKIPZIP is set to certain values
# Just do sed 's/\t/,/g' report.tsv > backup.csv to recreate your backup CSV copy
echo "Converting TSV -> CSV..."
rm report.tsv
sed 's/,/\t/g' report.csv > report.tsv
rm foo

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "Report made after $DT seconds"

# Close MAKEREPORT if statement
fi


# ---------------------------------------
# Calculate the tree strcture whatever happens, because it's quick and it's useful for later reports

echo "CALCULATING TREE STRUCTURE..."

# First we need to identify clades
echo "Identifying clades..."
gawk -v FS=, -v OFS=, '$1+0>0 && n==0 {n=1} $1+0==0 && n==1 {n=2} n==1 && $17==1 && $6>1 {print $1,$2,$5,$6,$16}' report.csv | \
gawk -v FS=, '{if (length($2)>0) {y=$2} else {y=$1}} $4==n && $5==f {x=x","y} ($4!=n || $5!=f) && $4>1 {print NR,n,f,n+f-1,x; x=y} {n=$4;f=$5} END {print NR,n,f,n+f-1,x}' > clades.csv

# Then we need to create the tree
echo "Forming tree..."
sort -nrk2 clades.csv | gawk '{a0[NR]=$3;a1[NR]=$4;split($5,b,",");a2[NR]=$6=b[1]; for (i=NR-1;i>=1;i--) {if ($3>=a0[i] && $4<=a1[i]) $6=a2[i]">"$6}; print}' | sort -nk3,3 -nk1 > tree.txt
sed 's/,/;/g' tree.txt | sed 's/ /,/g' > tree.csv


T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "Tree complete after $DT seconds"


# ---------------------------------------
# Make a short version of the report, compressing the information in singletons

if [ "$MAKESHORTREPORT" -gt "0" ]; then

echo "MAKING A SHORT REPORT..."

# Write report
head -12 report.csv > short-report.csv
gawk -v FS=, -v OFS=, '$6>1' variant-not-shared.txt >> short-report.csv
gawk -v FS=, -v OFS=, '$6==1' variant-not-shared.txt > variant-singletons.txt
echo '' >> short-report.csv
head -1 report.csv | awk '{$6="";$1="Singletons"; print}' > foo; cat foo >> short-report.csv
# This flags the singletons which are not called in adjecent tests
gawk -v FS=, -v OFS=, 'NR==FNR {start[FNR]=$3;end[FNR]=$4;n=FNR} NR!=FNR {for (i=1;i<=n;i++) if (start[i]<=$16 && end[i]>=$16) {s=start[i];e=end[i]}; neg=0; for (i=s;i<=e;i++) if (length($i)==0) neg++; flag=""; if (neg<e-s) flag="(s?)"; if (neg==0) flag="(s?!)"; $($16)=flag""$($16); print}' tree.csv variant-singletons.txt > foo
# This compresses them into a short format
gawk -v FS=, -v OFS=, -v maxn=0 '{n[$16]++;a[$16,n[$16]]=$($16); if (maxn<n[$16]) maxn=n[$16]} END {for (i=1;i<=maxn;i++) {for (j=1;j<=NF;j++) printf "%s,",a[j,i]; printf "\n"}}' foo >> short-report.csv
echo '' >> short-report.csv
head -1 report.csv | awk '{$6="";$1="Inconsistent SNPs"; print}' > foo; cat foo >> short-report.csv
echo 'GrCh37,Name(s),Ref,Alt,Type,N+,(?+),N-,nc,cbl+,cbl-,cbu+,cbu-,cblu+,cblu-,1stCol,Recur' >> short-report.csv
cat variant-bad.txt >> short-report.csv

# Close MAKESHORTREPORT if statement
fi

# ---------------------------------------
# Make a short version of the report, compressing the information in singletons

if [ "$MAKEHTMLREPORT" -gt "0" ]; then

echo "MAKING THE HTML REPORT..."

# Make HTML header
echo "<!DOCTYPE html>" > report.html
echo "<HTML>" >> report.html
echo "<BODY>" >> report.html

# Oopsie...
echo "<H1>THIS ISN'T FINISHED YET!</H1>" >> report.html


# Make table header
echo "<TABLE border=0 cellpadding=1>" >> report.html

# Insert table contents
 gawk -v FS=, '$2+0>0 {nsnps=split($5,snpnames,";"); print "<TD colspan="$2+0" rowspan="nsnps" bgcolor=\"#FFAAAA\" align=\"center\" alt=\""$3-17"\">"$5"</TD>"}' tree.csv | sed 's/;/<BR>/g' | grep 'alt="1"' | gawk '{print "<TR>"$0"</TR>"}' >> report.html

# Make table footer
echo "</TABLE>" >> report.html

# Make HTML footer
echo "</BODY>" >> report.html
echo "</HTML>" >> report.html

# Close MAKEHTMLREPORT if statement
fi

# ---------------------------------------
# Now let's move on to the age estimation

if [ "$MAKEAGES" -gt "0" ]; then

echo "CALCULATING AGES..."

# Then for each branch of the tree, we need to calculate the average number of SNPs, average coverage, and form a basic age from this
echo "Counting SNPs..."
gawk -v FS=, -v r="$RATE" 'NR==FNR {for (i=1;i<=6;i++) {c[FNR,i]=$i}; delete a; n=split($5,a,";"); s[FNR]=a[n]; n=FNR} \
			  NR!=FNR && $17=="Coverage" {for (i=18;i<=NF;i++) cov[i]=$i}
              NR!=FNR && $5=="SNP" {for (i=1;i<=n;i++) { \
			                if (flag[i]==1) {for (j=c[i,3];j<=c[i,4];j++) if (($j~$2 && length($2)>0) || ($j~$1 && $1+0>0)) {nsnp[i]++}};
							if ($1==s[i] || $2==s[i]) {for (j=c[i,3];j<=c[i,4];j++) coverage[i]+=cov[j]; flag[i]=1}}} \
              END {for (i=1;i<=n;i++) if (coverage[i]>0 && c[i,2]>0 && r>0) print i,nsnp[i],c[i,2]+0,nsnp[i]/c[i,2],nsnp[i]/(r*1e-9)/coverage[i],c[i,6]}' tree.csv header.csv variant-not-shared.txt > raw-ages.txt

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "Age analysis complete after $DT seconds"

# Close MAKEAGES if statement
fi

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "Done in $DT seconds"
