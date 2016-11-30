# !/bin/bash
# Author: Iain McDonald
# Contributors: Harald Alvestrand
# Purpose: Reduction and comparison script for Family Tree DNA's BigY data
# Release: 0.4.3 2016/10/27
# For free distribution under the terms of the GNU General Public License, version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# WARNING: This code is not yet ready for clades outside of U106. Please contact the author for details.

# Set up instructions:
#    Create sub-directory called zip, and insert BigY zip files into the folder with the format bigy-<NAME>-<KIT NUMBER>.zip
#    Create lists of implications, bad SNPs, recurrent SNPs and SNP names in the following files:
#		implications.txt badlist.txt recurrencies.txt cladenames.txt
#    The first three lists should obey the following formats:
#         7246726 > 23612197 : L48 > Z381
#    	  22436300 : E206
#         20723243 : F552
#	The first example implies that anyone who is L48+ is also Z381+. Only columns 1 and 3 are used.
#   The second example states that E206 is an inconsistent SNP. Only the first column is used.
#   The third example shows the (same) format for recurrent SNPs, indicating F552 is recurrent. Only the first column is used.
#	The final file, cladenames.txt, allows the substitution of automagically generated SNP names, e.g.:
#	
#
# Windows WSL BASH users:
#   Non-native package dependences:
#      unzip
#   These must be installed before the script is run, otherwise an error message will be generated. Try:
#      sudo apt-get install unzip
#
# Additional package dependencies (installed by default under Windows WSL):
#      gawk, sed, python
#
# Required internal libraries:
#      positives-and-no-calls.py : This code reads the call status of SNPs from the BED files (Courtesy H.A.)
#      snps_hg19.csv : Contains SNP names. Can be updated, e.g.: wget http://ybrowse.org/gbrowse2/gff/snps_hg19.csv
#      age.bed : Contains the regions used to count SNPs for age analysis.
#      poisson.tbl : A look-up table containing a matrix of Poisson functions
#      cpoisson.tbl : A look-up table containing 95% confidence intervals for cumulative Poisson statistics
#
# Create copy using:
#      zip redux.zip redux.bash positives-and-no-calls.py snps_hg19.csv age.bed poisson.tbl cpoisson.tbl implications.txt badlist.txt recurrencies.txt cladenames.txt
#
# For circa 700 kits, this process takes around 50 minutes on a heavy-duty 2014 laptop.

# Change log:
# 0.4.3.20161027a - Began encoding clade naming preferences.
# 0.4.2.20161025b - Completed final ("top-down") age analysis with uncertainties
# 0.4.1.20161025a - Completed "bottom-up" age analysis with uncertainties
# 0.4.0.20161002a - Began setup for uncertainty analysis in age calculations
# 0.3.1.20160929a - EARLY RELEASE 733
#					Replaced BED file scanning with python script from Harald A.: dramatic speed improvement
#                   Introduced date information to header
#					Introduced calculation for age analysis coverage, reported in header
#					Introduced coverage restriction for age calculation
#					Introduced full age analysis (no uncertainty ranges yet)
#					Introduced key at top-left
#					Introduced clade listings at top of tree (not currently part of output)
# 0.3.0.20160909a - EARLY RELEASE 728: Added second auto-backup for fresh runs of the script
#					Introduced number of unique SNPs into age calculation
#					Introduced "top-down" age normalisation to account for causality
# 0.2.2.20160905a - EARLY RELEASE 723: Fixed bug in shared SNP counting involving presumed positives being counted twice
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
# Basic report: include GrCh37<->GrCh38 conversion (using CrossMap?)
# Basic report: sorting criterion to match FTDNA tree
# Basic report: parallisation / further optimisation
# Basic report: test equivalencies of SNPs and indels to remove duplicates and test for any missing clades they form (e.g. FGC564 / 7775535)
# Basic report: create file to allow swaps of ancestral/derived values on reference chromosome
# Basic report: add support for SNPs which are derived in the reference sequence
# Basic report: artificially inserted SNPs to allow clade subdivision on the basis of other data (e.g. Z301, DF98)
# Basic report: additional columns:
#					Clade's primary SNP
#					Tree position
#				additional rows:
#					Tree position ("R1b1a1a2...")
#					Lowest shared clade ("terminal SNP") 
# HTML report: in the style of Alex W's BigTree
# General: BigY+YElite+WGS crossover

# This label will define what your top-level SNP is called out of the (usually large number of) possible shared SNPs.
TOPSNP="U106"

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
#    This requires a successful (but not necessarily sorted) previous run.
# If set to 2, statistics generation is skipped. The statistics INCLUDING KIT NAMES will be taken from the header of the previous run.
#    It also uses the ORDER of the previous run, taken from order.txt
#    This requires a complete successful previous run
# If set to 3, SNP names aren't updated either. The previous SNP names are used.
SKIPZIP="0"

# Display warning if inconsistencies detected in report (0=no; 1=yes)
WARNINGSNPS=1

# The following parameters set the mutation rate for the age analysis.
# The rate is set in SNP mutations per year per *billion* base pairs.
# RATE0 and RATE1 give the lower and upper bounds for the 95% confidence interval.
# The rate you should use varies depending on the region of the chromsome involved.
# Rates are numerically lower in palindromic regions (see Helgason et al. 2015 and our own analyses).

# Rates for entire BigY test
#RATE=0.751
#RATE0=0.688
#RATE1=0.814

# Rates for standardised age.bed restricted region
RATE=0.8181
RATE0=0.7586
RATE1=0.8767

# ZEROAGE gives the date from which ages are computed
# ERRZEROAGE gives the (95% c.i.) +/- variation in this date for individual testers
ZEROAGE=1950
ERRZEROAGE=16

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

# Make further backup copies when running the script from scratch
# This is useful when you want to make changes to the bad/inconsistent list, but still want to compare to the original previous run.
# For example:
# gawk 'NR==FNR {c[$5]++;next};c[$5]==0' tree.txt autobackup2/tree.txt
# will tell you the changes to the tree structure that have resulted from the addition of new kits between "from-scratch" runs.
if [ ! -d "autobackup2" ]; then
mkdir autobackup2
else
rm -rf autobackup2/*
fi
cp variant-*.txt autobackup2/
cp snps_hg19.csv autobackup2/
cp snp-names.csv autobackup2/
cp snp-used.csv autobackup2/
cp report.csv autobackup2/
cp short-report.csv autobackup2/
cp clades.csv autobackup2/
cp tree.txt autobackup2/
cp raw-ages.txt autobackup2/

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
echo 'KEY:,,,,,,,,,,,,,,,,Date' >> report.csv
echo 'N+/N-,Number of +/- calls,,,,,,,,,,,,,,,Coverage' >> report.csv
echo '(?+),Call uncertain but presumed positive,(position forced),,,,,,,,,,,,,,...for age analysis' >> report.csv
echo 'cbl,Occurs on lower boundary of coverage,(often problematic),,,,,,,,,,,,,,Regions' >> report.csv
echo 'cbu,Occurs on upper boundary of coverage,(usually ok),,,,,,,,,,,,,,Variants' >> report.csv
echo 'cblu,Occurs as a 1-base-pair region,,,,,,,,,,,,,,,Passed' >> report.csv
echo '1stCol,First column which is positive,,,,,,,,,,,,,,,Simple SNPs' >> report.csv
echo 'Recur,Recurrencies in tree,(check: 1 or (R)),,,,,,,,,,,,,,SNPs under U106' >> report.csv
echo '(s?),Questionable singleton,(not negative in some clademates),,,,,,,,,,,,,,Singleton SNPs' >> report.csv
echo '(s?!),Questionable singleton,(not negative in all clademates),,,,,,,,,,,,,,...for age analysis' >> report.csv
echo '(R),Allowed recurrency,,,,,,,,,,,,,,,Indels' >> report.csv
echo 'Blank,Securely called negative,,,,,,,,,,,,,,,Indels under U106' >> report.csv
echo 'Full report at:,www.jb.man.ac.uk/~mcdonald/genetics/report.csv,,,,,,,,,,,,,,,Singleton Indels' >> report.csv
echo 'Non-shared SNPs' >> report.csv
echo 'GrCh37,Name(s),Ref,Alt,Type,N+,(?+),N-,nc,cbl+,cbl-,cbu+,cbu-,cblu+,cblu-,1stCol,Recur' >> report.csv
echo "Generating statistics for" ${#FILES[@]} "BED files..."
for BEDFILE in ${FILES[@]}; do
	VCFFILE=`echo "$BEDFILE" | sed 's/.bed/.vcf/'`
	KITNAME=`echo "$BEDFILE" | gawk -F/ '{print $2}' | sed 's/.bed//'`
	KITDATE=`ls -l --time-style +%Y-%m-%d "$BEDFILE" | cut -d\  -f6`
	#STATS=`gawk '{s+=$3-$2-1} END {print k,s,NR}' k="$KITNAME" "$BEDFILE"`
	STATS=`gawk 'NR==FNR {a[NR]=$1;b[NR]=$2;n=NR} NR!=FNR {s+=$3-$2-$1; for (i=1;i<=n;i++) if ($2<=b[i] && $3>=a[i]) {x=($3>b[i]?b[i]:$3)-($2>a[i]?$2:a[i])+1; if (x<0) x=0; as+=x}} END {print s,as,FNR}' age.bed "$BEDFILE"`
	STATS2=`gawk '$1=="chrY" {n++} $1=="chrY" && $7=="PASS" {v++; if ($4!="." && $5!=".") {if (length($4)==1 && length($5)==1) {s++} else {i++}}} END {print n,v,s,0,0,0,i,0,0}' "$VCFFILE"`
	STATS="$KITNAME $KITDATE $STATS $STATS2"
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
# Include python script by Harald A.
./positives-and-no-calls.py ${FILES[@]} > variant-match.txt

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
                         for (j=1;j<=l;j++) {if (s[j]==b[i]) {for (k=18;k<=m;k++) if (length(d[j,k])<=5 && t[k]==1 && (substr(d[j,k],1,1)==";" || length(d[j,k]==0))) d[j,k]="(?+)"d[j,k]}}}; \
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
# The previous set of names is used if SKIPZIP is set to above 2
echo "Introducing SNP names..."
if [ "$SKIPZIP" -le "2" ]; then
	gawk -v FS=, -v OFS=, 'NR>1 {print $4"."$11"."$12,$9}' snps_hg19.csv | sed 's/"//g' | sort -nk1,1 -t. | gawk -v FS=, -v OFS=, '$1==a {b=b"/"$2; print "...",$0} $1!=a {print a,b;a=$1;b=$2} END {print}' > snp-names.csv
	gawk -v FS=, '{print $1,$1"."$3"."$4}' variant-output.txt | sort -nk1 | gawk '{print $2}' > foo
	gawk -v FS=, 'NR==FNR {a[FNR]=$1;na++} NR!=FNR {b[FNR]=$1;bb[FNR]=$0;nb++} END {ij=1; for (i=1;i<=na;i++) {for (j=ij;j<=nb;j++) {if (a[i]==b[j]) {print bb[j]; ij=j; break}; xa=substr(a[i],1,8); xb=substr(b[i],1,8); if (xb+0>xa+0) {break}}}}' foo snp-names.csv > snp-used.csv
fi
gawk -v FS=, -v OFS=, 'NR==FNR {s[NR,1]=$1;s[NR,2]=$2;split($2,arr," ");s[NR,3]=arr[1];n=NR} NR!=FNR {m=$1"."$3"."$4; for (i=1;i<=n;i++) if (s[i,1]==m) {$2=$2""s[i,2]; gsub(m,s[i,3],$0)}; print}' snp-used.csv variant-not-shared.txt > foo; mv foo variant-not-shared.txt
gawk -v FS=, -v OFS=, 'NR==FNR {s[NR,1]=$1;s[NR,2]=$2;split($2,arr," ");s[NR,3]=arr[1];n=NR} NR!=FNR {m=$1"."$3"."$4; for (i=1;i<=n;i++) if (s[i,1]==m) {$2=$2""s[i,2]; gsub(m,s[i,3],$0)}; print}' snp-used.csv variant-shared.txt > foo; mv foo variant-shared.txt
gawk -v FS=, -v OFS=, 'NR==FNR {s[NR,1]=$1;s[NR,2]=$2;split($2,arr," ");s[NR,3]=arr[1];n=NR} NR!=FNR {m=$1"."$3"."$4; for (i=1;i<=n;i++) if (s[i,1]==m) {$2=$2""s[i,2]; gsub(m,s[i,3],$0)}; print}' snp-used.csv variant-bad.txt > foo; mv foo variant-bad.txt

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

# Count the number of SNPs and indels, count the singletons, and count the singletons for any later age analysis
echo "Creating SNP counts..."
U106SNPS=`gawk -v FS=, -v OFS=, '$5=="SNP" {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
U106INDELS=`gawk -v FS=, -v OFS=, '$5=="Indel" {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
SINGSNPS=`gawk -v FS=, -v OFS=, '$5=="SNP" && $6==1 {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
SINGINDELS=`gawk -v FS=, -v OFS=, '$5=="Indel" && $6==1 {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {n[i]++}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' variant-not-shared.txt`
AGESINGSNPS=`gawk -v FS=, -v OFS=, 'NR==FNR {split($0,ax," ");a[NR]=ax[1];b[NR]=ax[2];na=NR} NR!=FNR && $5=="SNP" && $6==1 {s=$2; for (i=1;i<=NF-17;i++) if (($(i+17)~$2 && length($2)>0) || ($(i+17)~$1 && length($2)==0) || $(i+17)~/(R)/) {for (j=1;j<=na;j++) if ($1>=a[j] && $1<=b[j]) {n[i]++}}} END {for (i=1;i<=NF-17;i++) {x+=n[i]; printf "%i,",n[i]}}' age.bed variant-not-shared.txt`
gawk -v FS=, -v OFS=, -v us="$U106SNPS" -v ui="$U106INDELS" -v ss="$SINGSNPS" -v si="$SINGINDELS" -v as="$AGESINGSNPS" '\
	NR==9 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,us} \
	NR==10 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,ss} \
	NR==11 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,as} \
	NR==13 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,ui} \
	NR==14 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,si} \
	NR<9 || NR==12 || NR>14 {print}' report.csv > foo; mv foo report.csv

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
gawk -v FS=, -v OFS=, -v lead="$TOPSNP" 'NR==1 {print 1,lead,"SNP",$6,18} $17+0>0 && n==0 {n=1} $17+0==0 && n==1 {n=2} n==1 && $17==1 && $6>1 {print $1,$2,$5,$6,$16}' report.csv | \
gawk -v FS=, 'NR==1 {x=$2;w=0} {if (length($2)>0) {y=$2} else {y=$1}; z=$1} $4==n && $5==f {x=x","y;w=w","z} ($4!=n || $5!=f) && $4>1 && NR>1 {print NR,n,f,n+f-1,x,w; x=y; w=z} {n=$4;f=$5} END {print NR,n,f,n+f-1,x,w}' > clades.csv

# Now let's count the number of SNPs in those clades that are to be used for age analysis purposes
gawk 'NR==FNR {split($0,a," ");counta[FNR]=a[1];countb[FNR]=a[2];nc=FNR} NR!=FNR {n=split($6,s,","); nsnp=0; for (i=1;i<=n;i++) {for (j=1;j<=nc;j++) if (s[i]>=counta[j] && s[i]<=countb[j]) nsnp++}; print $0,nsnp}' age.bed clades.csv | sort -nrk2 > cladecounts.csv

# Then we need to create the tree
echo "Forming tree..."

# Clade name determination (first gawk command)
# First, replace the separator for alternative designations (/) with that for equivalent SNPs (,)
# Create arrays [b,ba,b0] consisting of the list of SNPs, their alphabetical parts and their numerical parts
# Take the first SNP ($6==b[1]) as an initial name. Then try to find a better one.
# Preference order (if a better SNP is found at a later step, it gets overwritten):
# 3. Take the first named SNP
# 2. Seperate the SNP names into classes:
#	(e) If there is a Y name, use the lowest numbered SNP;
#	(d) If there is a BY or A name, use the lowest numbered SNP;
#	(c) If there is an FGC, DF, PF, Z, CTS or S name, use the lowest numbered SNP;
#	(b) If there is an DF name, use the lowest numbered SNP;
#	(a) If there is a M, P, U, L name, use the lowest numbered SNP.
# 1. If there is a manual override requested, use that instead (a check is first done to see if it is actually in the list).
gawk 'NR==FNR {autoname[NR]=$1;newname[NR]=$2;nrename=NR} \
NR!=FNR {a=$5; gsub("/",",",a); aa=a0=a; gsub(/[0-9]/,"",aa); gsub(/[A-Za-z]/,"",a0); n=split(a,b,","); split(aa,ba,","); split(a0,b0,","); \
$6=b[1]; lownum=9e9; \
if ($6+0>0) for (i=1;i<=n;i++) if (b[i]>0) {$6=b[i]; break}; \
for (i=1;i<=n;i++) if (ba[i]=="Y" && b0[i]<lownum) {$6=b[i]; lownum=b0[i]}; \
foundclass=0; for (i=1;i<=n;i++) if (ba[i]=="BY" || ba[i]=="A") {foundclass=1}; \
	if (foundclass==1) {lownum=9e9; for (i=1;i<=n;i++) if ((ba[i]=="BY" || ba[i]=="A") && b0[i]<lownum) {$6=b[i]; lownum=b0[i]}}; \
foundclass=0; for (i=1;i<=n;i++) if (ba[i]=="FGC" || ba[i]=="PF" || ba[i]=="Z" || ba[i]=="CTS" || ba[i]=="S") {foundclass=1}; \
	if (foundclass==1) {lownum=9e9; for (i=1;i<=n;i++) if ((ba[i]=="FGC" || ba[i]=="PF" || ba[i]=="Z" || ba[i]=="CTS" || ba[i]=="S") && b0[i]<lownum) {$6=b[i]; lownum=b0[i]}}; \
foundclass=0; for (i=1;i<=n;i++) if (ba[i]=="DF") {foundclass=1}; \
	if (foundclass==1) {lownum=9e9; for (i=1;i<=n;i++) if ((ba[i]=="DF") && b0[i]<lownum) {$6=b[i]; lownum=b0[i]}}; \
foundclass=0; for (i=1;i<=n;i++) if (ba[i]=="M" || ba[i]=="P" || ba[i]=="U" || ba[i]=="L") {foundclass=1}; \
	if (foundclass==1) {lownum=9e9; for (i=1;i<=n;i++) if ((ba[i]=="M" || ba[i]=="P" || ba[i]=="U" || ba[i]=="L") && b0[i]<lownum) {$6=b[i]; lownum=b0[i]}}; \
for (i=1;i<=nrename;i++) if ($6==autoname[i]) {for (j=1;j<=n;j++) if (b[j]==newname[i]) $6=newname[i]}; \
print}' cladenames.txt cladecounts.csv | \
gawk -v lead="$TOPSNP" '{a0[NR]=$3;a1[NR]=$4;a2[NR]=$6; if (NR==1) a2[NR]=""; for (i=NR-1;i>=1;i--) {if ($3>=a0[i] && $4<=a1[i]) $6=a2[i]">"$6}; $6=lead""$6; if (NR==1) {$7=0;$5=lead;$6=lead}; print}' | \
sort -nk3,3 -nk1 | \
gawk 'NR==1 {n[NR]=1;sh[1]="0"} \
      NR>1 {n[NR]=split($6,tt,">"); for (i=1;i<=n[NR];i++) t[NR,i]=tt[i]; shn=0; sh[NR]=0; for (i=1;i<=NR;i++) if (t[NR,n[NR]-1]==t[i,n[i]]) {sh[NR]=sh[i]}; for (i=1;i<=NR;i++) if (t[NR,n[NR]-1]==t[i,n[i]-1]) {shn++}; sh[NR]=sh[NR]"."shn} \
	  {$6=n[NR]" "sh[NR]" "$6; print}'  > tree.txt
sed 's/,/;/g' tree.txt | sed 's/ /,/g' > tree.csv

#sort -nrk2 cladecounts.csv | gawk -v lead="$TOPSNP" '{a0[NR]=$3;a1[NR]=$4;split($5,b,",");a2[NR]=$6=b[1]; if (NR==1) a2[NR]=""; for (i=NR-1;i>=1;i--) {if ($3>=a0[i] && $4<=a1[i]) $6=a2[i]">"$6}; $6=lead""$6; if (NR==1) {$7=$5;$5=".";$6=lead}; print}' | 
#gawk -v lead="$TOPSNP" '{a0[NR]=$3;a1[NR]=$4;a2[NR]=$6; if (NR==1) a2[NR]=""; for (i=NR-1;i>=1;i--) {if ($3>=a0[i] && $4<=a1[i]) $6=a2[i]">"$6}; $6=lead""$6; if (NR==1) {$7=$5;$5=".";$6=lead}; print}' | \

# Performing consistency check
gawk -v FS=, '$17>1 && $2 !~ /(R)/ && $6>1 {print $1,$2,$5,$6,$17,"Inconsistent: conflicting calls"}' variant-not-shared.txt > warning-list.txt
gawk -v FS=, '$17>1 && $2 !~ /(R)/ && $6==1 {print $1,$2,$5,$6,$17,"Inconsistent: multiple calls as singleton"}' variant-not-shared.txt >> warning-list.txt
gawk -v FS=, '$6>1 && $17==1 {first=0; last=0; for (i=NF;i>=18;i--) if ($i~$1 || ($i~$2 && length($2)>0) || $i~/(?+)/) first=i; for (i=18;i<=NF;i++) if ($i~$1 || ($i~$2 && length($2)>0) || $i~/(?+)/) last=i; if (last-first+1>$6+0) print $1,$2,$6,$17,"Inconsistent: out of order"}' variant-not-shared.txt >> warning-list.txt
WARNINGSNPS=`wc -l warning-list.txt | gawk '{print $1}'`
if [ "$WARNINGSNPS" -gt "0" ]; then
echo "WARNING! $WARNINGSNPS mutations are not consistently placed in the tree."
echo "Set SKIPZIP=3 and fix them using either badlist.txt or implications.txt."
echo "The inconsistent mutations are recorded in warning-list.txt."
fi

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "Tree complete after $DT seconds"



# ---------------------------------------
# Insert tree into report (if being made)
#if [ "$MAKEREPORT" -gt "0" ]; then
#head -1 report.csv > foo
#gawk -v FS=, '{split($8,a,">");if (NR==1) {n1=$3; n2=$4}; max=$6>max?$6:max; for (i=$3;i<=$4;i++) c[i,$6]=a[$6]} END {printf "Clades"; for (y=1;y<=max;y++) {for (x=1;x<=n1-1;x++) printf ","; for (x=n1;x<=n2;x++) printf "[%s],",c[x,y]; print ""}}' tree.csv > tree-report.csv
#echo "" >> foo
#cat tree-report.csv >> foo
#echo "" >> foo
#gawk -v FS=, 'NR>1' report.csv >> foo
#mv foo report.csv
# Close second MAKEREPORT if statement
#fi


# ---------------------------------------
# Make a short version of the report, compressing the information in singletons

if [ "$MAKESHORTREPORT" -gt "0" ]; then

echo "MAKING A SHORT REPORT..."

# Write report
gawk -v FS=, 'n==0 {print} $1=="Non-shared SNPs" {n=1}' report.csv > short-report.csv
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
# Make an HTML version of the report

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
# XXX STILL TO DO XXX
 
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

# For each branch of the tree, we need to calculate:
# the average number of SNPs (nsnp[i]/c[i,2]),
# the total number of unique SNPs (snpsused[i]),
# average coverage (coverage[i]/c[i,2]),
# and form a basic age from this (nsnp[i]/<rate>/<coverage>)
# NB: flag[1]=1 sets the top-level SNP to always be counted
# Uncertainties: the lower and upper 95% confidence intervals are set by:
#  Lower: m/<rate>/<coverage> where int_0^n Poisson(n,m) = 0.025 for n observed SNPs and a theoretical mean of m
#  Upper: m/<rate>/<coverage> where int_0^n Poisson(n,m) = 0.975 for n observed SNPs and a theoretical mean of m
echo "Counting SNPs..."
gawk -v FS=, -v r="$RATE" 'NR==FNR {for (i=1;i<=9;i++) {c[FNR,i]=$i}; delete a; n=split($5,a,";"); natlevel[FNR]=n; s[FNR]=a[n]; n=FNR} \
			  FILENAME=="header.csv" && FNR==4 {for (i=18;i<=NF;i++) {cov[i]=$i; covtotal+=$i}} \
			  FILENAME=="age.bed" {split($0,a," ");counta[FNR]=a[1];countb[FNR]=a[2];nc=FNR} \
              NR!=FNR {flag[1]=1; for (i=1;i<=n;i++) { \
			           if (flag[i]==1 && $5=="SNP") {snpused=0; for (j=c[i,3];j<=c[i,4];j++) if (($j~$2 && length($2)>0) || ($j~$1 && $1+0>0)) {for (k=1;k<=nc;k++) if ($1>=counta[k]&&$1<=countb[k]) {nsnp[i]++; snpused=1}}; snpsused[i]+=snpused};
					   if ($1==s[i] || $2==s[i]) {for (j=c[i,3];j<=c[i,4];j++) coverage[i]+=cov[j]; flag[i]=1}}} \
              END {coverage[1]=covtotal; for (i=1;i<=n;i++) {if (coverage[i]>0 && c[i,2]>0 && r>0) print i,nsnp[i]+0,c[i,2]+0,snpsused[i],nsnp[i]/c[i,2],coverage[i]/c[i,2],natlevel[i],c[i,6],c[i,3],c[i,4],c[i,7],c[i,8],c[i,9]}}' tree.csv header.csv age.bed variant-not-shared.txt > raw-ages.txt

# The following line consults a Poisson statistics look-up table
# For the specified number of standard deviations (s) it will take the number
# of unique SNPs in that line ($2) as the mean and identify the Poisson interval
# corresponding to that number of standard deviations from the mean.
# It will then take that interval and multiply it by the marginalised rate with
# central value $RATE, and lower/upper confidence intervals at <s> standard 
# deviations $RATE0 .. $RATE1. The age and confidence interval are appended
# as three extra columns.
# CPoisson.tbl only covers numbers 0-100. For larger numbers a Gaussian is used
# to approximate the Poisson distribution, corrected by a small addition
# <gausscorn> and <gausscorp> listed in the last column of poisson.tbl.
gawk -v s=1.96 -v r="$RATE" -v r0="$RATE0" -v r1="$RATE1" ' \
	NR==FNR && ns<-s && $1>-s && NR>1 {numf=split(sd,negd1," "); split($0,negd2," "); comf=(-s-ns)/($1-ns); \
		for (i=1;i<=numf-1;i++) negd[i]=negd1[i]+(negd2[i]-negd1[i])*comf; gausscorn=negd1[numf]+(negd2[numf]-negd1[numf])*comf; negd0=negd1[3]+(negd2[3]-negd1[3])*comf} \
	NR==FNR && ns<s && $1>s && NR>1 {numf=split(sd,posd," "); split($0,posd2," "); comf=(s-ns)/($1-ns); \
		for (i=1;i<=numf;i++) posd[i]=posd1[i]+(posd2[i]-posd1[i])*comf; gausscorp=posd1[numf]+(posd2[numf]-posd1[numf])*comf; posd0=posd1[3]+(posd2[3]-posd1[3])*comf} \
	NR==FNR {ns=$1;sd=$0} \
	NR!=FNR {if ($2==0) {print $0,0.67/$3/r*1e9/$6,posd0/$3/r1*1e9/$6,negd0/$3/r0*1e9/$6} else \
	         if ($2>0 && $2<=numf-4) {print $0,$5/r*1e9/$6,$5/r1*1e9/$6*posd[$2+3]/$2,$5/r0*1e9/$6*negd[$2+3]/$2} else \
				{print $0,$5/r*1e9/$6,$5/r1*1e9/$6*($2-s*sqrt($2)+gausscorn)/$2,$5/r0*1e9/$6*($2+s*sqrt($2)+gausscorp)/$2}}' cpoisson.tbl raw-ages.txt > raw-ages-err.txt
			  
T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "...complete after $DT seconds"

echo "Consolidating ages..."

# This step is a complicated mess of stastistical trickery that aims to better account for two things:
# (1) Causality: sub-clades must have formed after their parent clades.
# (2) Normalisation: since SNPs occur randomly, some subclades have more or fewer than others.
#     If the only surviving line in a subclade has a larger/smaller number of SNPs than average, these will be passed onto ALL their descendants,
#     so the line will have a higher/lower average number of SNPs than its "brother" clades which share the same parent.
# There are two approaches to this. Both are used here.
# TOP-DOWN:  An age is computed for the top-level SNP, which should best represent the average computed mutation rate.
#            Clades below the top level have their mutation rates normalised so that they match up to this age.
# BOTTOM-UP: (cf. YFull's method). First, an age for each bottom-level clade is computed.
#            Ages for parent clades are computed from the weighted averaage of ages of their sub-clades.
#            In this case, the weighting is done using the square root of the number of unique SNPs in that sub-clade.
# An advantage of the top-down method are that causality is automatically implemented. A disadvantage is that ages can be biased if a single sub-clade dominates.
# An advantage of the bottom-up method are that random uncertainties can more readily accounted for in older clades.
# Disadvantages of the bottom-up methodology are that judgement calls must be made in places where causality fails, 
# and ages for smaller clades are less likely to be accurately determined as random variations in the SNP counts aren't normalised out.

# Get the lowest tree level
# and the coverage of each kit from the report
# and the number of age-countable singletons from the report
MAXDEPTH=`sort -t, -nrk6,6 tree.csv | head -1 | awk -v FS=, '{print $6}'`
COVERAGE=`gawk -v FS=, '$17=="...for age analysis" {n++} n==1 {print; n++}' report.csv`
SINGLES=`gawk -v FS=, '$17=="...for age analysis" {n++} n==2 {print; n++}' report.csv`
# Iterate up the tree, merging branches as I go
# In this case, the weights are simply the number of tests that are within that clade, plus a subtle shift (0.67 SNPs for each clade) to account for Poisson statistics

# The full calculation with uncertainties requires the full propagation of the probability distribution function (PDF).
# In short, the bottom-up methodology merges branches from the bottom to the top of the tree.
# The top-down methodology then restricts the PDF of sub-clades based on the PDF of the parent.
# For clades with no sub-clades:
# 	1. The number of singletons [S] represents a point on an inverse Poisson function.
# 	2. The probability of finding [X] SNPs, given a mean of [S] is =POISSON(X,S,0)
# 	3. An array for different values of [X] and [S] is provided in poisson.tbl, which can be interpolated to provide a PDF.
# 	4. Each PDF is convolved with the PDF of the SNP formation rate to identify an age.
# 	5. The convolved PDFs of each individual in the clade are multiplied together to form the global PDF.
# 	6. The cumulative (integral) of this PDF is calculated, and the 50th centile position is taken as the most likely age.
# 	7. The points at probabilities corresponding to -[s] and +[s] sigma define the uncertainty range.
# For clades with sub-clades, the time taken to produce that clade's set of SNPs must be added to the age of the sub-clade.
# 	1. The PDF already generated for the sub-clade is used as a starting point.
# 	2. This is additively convolved with the PDF arising from the number of shared SNPs at that clade level.
# 	3. Then the age for the parent clade is calculated as before.
# Poisson tables are read in from poisson.tbl. POISSON(X,0,0) and POISSON(X,1,0) are listed.
# 	POISSON(X,S,0) can then be calculated using the relationship:
# 	POISSON(X,S,0) / MAX[POISSON(_,S,0)] = POISSON(X/S,1,0)^S / MAX[POISSON(_,1,0)^S]
#	where MAX[POISSON(_,1,0)^S] = 0.3678794.
# This process is completed as follows:
# 	1. The number of singletons, coverage and mutation rates are passed from the command line as arrays.
# 	2. The Poisson table is read in from disc [poisson.tbl, NR==FNR].
# 	3. Then the tree structure is read in from disc [raw-ages-err.txt, NR!=FNR].
# 	4. Once all data is read in [END] then begin the bottom-up analysis.
#		Set some parameters for the top-level clade,
#		and loop over all tree levels [i] (bottom to top) and over all clades at that level [j]
#		(a) Reset arrays for sub-clades and reset ages for each clade.
#		(b) Loop over all members of that clade. Assign them either to:
#			[flag[k]==0] Not belonging to a sub-clade
#			[flag[k]==1] Belonging to a sub-clade, so ignored in the subsequent analysis
#			[flag[k]==2] Belonging to a sub-clade, and being the (arbitrary) representative for that sub-clade
# 			In cases where flag[k]==2, also register information for that sub-clade in the appropriate arrays
#			(i) Read in the already-computed PDF of that sub-clade
# 			(ii) Calculate the PDF for the SNPs between the sub-clade and the test clade,
# 			(iii) Perform the convolution (i) + (ii) to provide the PDF [kprob].
#			In cases where flag[k]==0, create a PDF [kprob] from Poisson statistics over times [t=0..maxt] in intervals [dt]
#			...based on three cases:
#				[single[k]==0], when [pois0] can be used,
#				[single[k]==1], when [pois1] can be used,
#				[single[k]>1], when the PDF must be calculated from [pois1].
#		(c) Multiply the resultant probabilities [kprob] together to get the PDF for the parent clade [subpdf]
#		(d) Normalise the PDF
#	5. Now begin the top-down analysis.
#		Loop over all tree levels [i] (top to bottom, from level 2) and over all clades at that level [j]
#		(a) Identify the parent and read in its PDF.
#		(b) Calculate the PDF for the number of SNPs in between the parent and the sub-clade.
#		(c) Perform the convolution (a) - (b)
#		(d) Multiply the sub-clade's PDF by that convolution
#	6. For each clade...
#		(a) Convolve with the rate uncertainty, according to [$RATE0] and [$RATE1] (*)
#		(b) Calculate the confidence interval according to [ciprob].
#		(c) Print out results for each clade.
#(*) Not done fully: age2n[j]/age2[j]/age2p[j] (2.5, 50, 97.5% confidence limits) reports this convolution
#	 but it is not applied to the listed PDFs [subpdf].
gawk -v rate="$RATE" -v r0="$RATE0" -v r1="$RATE1" -v max="$MAXDEPTH" -v singles="$SINGLES" -v coverages="$COVERAGE" -v ciprob=0.95 -v maxt=10000 -v dt=10 '\
	NR==1 {p0=(1-ciprob)/2; p1=1-p0}; \
	NR==1 {split(coverages,coverage,","); split(singles,single,",")} \
	NR==FNR {prob[NR]=$1; pois0[NR]=$2; pois1[NR]=$3/0.3678794}; \
	NR!=FNR {r[NR]=rate; for (i=1;i<=NF;i++) inp[NR,i]=$i; n=split(inp[NR,12],clades,">"); parent[NR]=clades[n-1]; clade[NR]=clades[n]} \
	END {inp[1,9]=18; inp[1,13]=0; parent[1]="."; for (i=max;i>=1;i--) for (j=1;j<=NR;j++) if (inp[j,8]==i) \
	         {nsubs=0; delete kprob; sage=sweight=0; age2n[j]=0; age2[j]=age2p[j]=maxt; \
			  for (k=inp[j,9];k<=inp[j,10];k++) {flag[k]=0; \
			     for (l=1;l<=NR;l++) {if (inp[l,9]>=inp[j,9] && inp[l,10]<=inp[j,10] && inp[l,8]+0==inp[j,8]+1 && k>=inp[l,9] && k<=inp[l,10]) \
				     {flag[k]=1; if (k==inp[l,9]) \
						{flag[k]=2; nsubs++; subatlevel=inp[l,13]+0; delete kprob0; delete kprob1; \
						nt=0; for (t=0;t<=maxt;t+=dt) {nt++; kprob0[nt]=subpdf[l,nt]}; \
						nt=0; for (t=0;t<=maxt;t+=dt) {nt++; nexpect=t*inp[l,6]*rate/1e9; \
							if (subatlevel==0) {plookup=int(nexpect*100)+1; if (plookup<1 || plookup>1000) {kprob1[nt]=0} else 
								{kprob1[nt]=(nexpect-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois0[plookup+1]-pois0[plookup])+pois0[plookup]}} \
							else if (subatlevel==1) {plookup=int(nexpect*100)+1; if (plookup<1 || plookup>1000) {kprob1[nt]=0} else 
								{kprob1[nt]=(nexpect-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois1[plookup+1]-pois1[plookup])+pois1[plookup]}} \
							else {plookup=int(nexpect/subatlevel*100)+1; if (plookup<1 || plookup>1000) {kprob1[nt]=0} else
								{kprob1[nt]=(nexpect/subatlevel-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois1[plookup+1]**subatlevel-pois1[plookup]**subatlevel)+pois1[plookup]**subatlevel}}}; \
						nt=0; for (t=0;t<=maxt;t+=dt) {nt++; nt1=0; if (kprob0[nt]!=0) for (t1=0;t1<=maxt-t;t1+=dt) {nt1++; nt0=nt+nt1-1; kprob[nsubs,nt0]+=kprob0[nt]*kprob1[nt1]}}}}}; \
				 if (k+0>0 && flag[k]==0) {nsubs++; \
					nt=0; for (t=0;t<=maxt;t+=dt) {nt++; nexpect=t*coverage[k]*rate/1e9; \
					if (single[k]==0) {plookup=int(nexpect*100)+1; if (plookup<1 || plookup>1000) {kprob[nsubs,nt]=0} else 
						{kprob[nsubs,nt]=(nexpect-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois0[plookup+1]-pois0[plookup])+pois0[plookup]}} \
					else if (single[k]==1) {plookup=int(nexpect*100)+1; if (plookup<1 || plookup>1000) {kprob[nsubs,nt]=0} else 
						{kprob[nsubs,nt]=(nexpect-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois1[plookup+1]-pois1[plookup])+pois1[plookup]}} \
					else {plookup=int(nexpect/single[k]*100)+1; if (plookup<1 || plookup>1000) {kprob[nsubs,nt]=0} else
						{kprob[nsubs,nt]=(nexpect/single[k]-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois1[plookup+1]**single[k]-pois1[plookup]**single[k])+pois1[plookup]**single[k]}}}}}; \
				 normprob=0; nt=0; for (t=0;t<=maxt;t+=dt) {nt++; subpdf[j,nt]=1; for (m=1;m<=nsubs;m++) {subpdf[j,nt]*=kprob[m,nt]}; normprob+=subpdf[j,nt]}; \
				 nt=0; for (t=0;t<=maxt;t+=dt) {nt++; subpdf[j,nt]/=normprob}; \
				 }; \
			 for (i=2;i<=max;i++) for (j=1;j<=NR;j++) if (inp[j,8]==i) {for (l=1;l<=NR;l++) if (clade[l]==parent[j]) jparent=l; \
				nt=0; for (t=0;t<=maxt;t+=dt) {nt++; parentpdf[nt]=subpdf[jparent,nt]}; \
				subatlevel=inp[j,13]+0; delete kprob1; nt=0; for (t=0;t<=maxt;t+=dt) {nt++; nexpect=t*inp[j,6]*rate/1e9; \
				if (subatlevel==0) {plookup=int(nexpect*100)+1; if (plookup<1 || plookup>1000) {kprob1[nt]=0} else 
					{kprob1[nt]=(nexpect-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois0[plookup+1]-pois0[plookup])+pois0[plookup]}} \
				else if (subatlevel==1) {plookup=int(nexpect*100)+1; if (plookup<1 || plookup>1000) {kprob1[nt]=0} else 
					{kprob1[nt]=(nexpect-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois1[plookup+1]-pois1[plookup])+pois1[plookup]}} \
				else {plookup=int(nexpect/subatlevel*100)+1; if (plookup<1 || plookup>1000) {kprob1[nt]=0} else
					{kprob1[nt]=(nexpect/subatlevel-prob[plookup])/(prob[plookup+1]-prob[plookup])*(pois1[plookup+1]**subatlevel-pois1[plookup]**subatlevel)+pois1[plookup]**subatlevel}}};
				delete tprob; nt=0; for (t=0;t<=maxt;t+=dt) {nt++; nt1=0; if (kprob0[nt]!=0) for (t1=0;t1<=t;t1+=dt) {nt1++; nt0=nt-nt1-1; if (nt0>0) tprob[nt0]+=parentpdf[nt]*kprob1[nt1]}}; \
				nt=0; for (t=0;t<=maxt;t+=dt) {nt++; subpdf[j,nt]*=tprob[nt]};
				normprob=0; nt=0; for (t=0;t<=maxt;t+=dt) {nt++; normprob+=subpdf[j,nt]}; \
				nt=0; for (t=0;t<=maxt;t+=dt) {nt++; subpdf[j,nt]/=normprob}}; \
			 for (i=1;i<=max;i++) for (j=1;j<=NR;j++) if (inp[j,8]==i) \
				{cumprob=0; nt=0; for (t=0;t<=maxt;t+=dt) {nt++; lastprob=cumprob; cumprob+=subpdf[j,nt]; \
					if (cumprob>p0 && lastprob<=p0) age2n[j]=(p0-lastprob)/(cumprob-lastprob)*dt+(t-dt); \
					if (cumprob>0.5 && lastprob<=0.5) age2[j]=(0.5-lastprob)/(cumprob-lastprob)*dt+(t-dt); \
					if (cumprob>p1 && lastprob<=p1) age2p[j]=(p1-lastprob)/(cumprob-lastprob)*dt+(t-dt)}; \
					age2n[j]*=rate/r1; age2p[j]*=rate/r0; \
					printf "%i ",i; if (length(parent[j])==0) parent[j]="."; for (n=1;n<=13;n++) printf "%s ",inp[j,n];; printf "%s %s %s %s %s %s | ",parent[j],clade[j],age[j],age2n[j],age2[j],age2p[j]; \
				nt=0; for (t=0;t<=maxt;t+=dt) {nt++; printf "%s ",subpdf[j,nt]}; printf "\n"}}' poisson.tbl raw-ages-err.txt > final-ages.txt

# Create the final HTML table
# 	Strip out the information, discard the PDF
#	Convert ages to dates
#	Correct for year zero problem, and replace +/- with AD/BC
#	Sort into "ISOGG" order
#	Encode HTML
gawk '{print $1}' FS="|" final-ages.txt | \
gawk -v z="$ZEROAGE" -v ez="$ERRZEROAGE" -v OFS="\t" '{print $12,$16,z+ez-int($17+.5),z-int($18+.5),z-ez-int($19+.5),$1}' | \
gawk '$3<1 {$3-=2; $3=-$3" BC"} $4<1 {$4-=2; $4=-$4" BC"} $5<1 {$5-=2; $5=-$5" BC"} $3!~"BC" {$3=$3" AD"} $4!~"BC" {$4=$4" AD"} $5!~"BC" {$5=$5" AD"} {print}' FS='\t' OFS='\t' | \
sort -rk1 | \
gawk 'NR==1 {print "<HTML><BODY><TABLE border=\"0\">"; print "<TR><TH>Clade\t<TH>Best guess<TH>95% confidence interval"} {printf "<TR><TD><SPAN style=\"display:inline-block; width:%ipx\"></SPAN>%s \t<TD>%s\t<TD>(%s &mdash; %s)\n",$6*10,$2,$4,$5,$3} END {print "</TABLE></BODY></HTML>"}' FS='\t' > table.html


T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "Age analysis complete after $DT seconds"

# Close MAKEAGES if statement
fi

T1=`date +%s.%N`
DT=`echo "$T1" "$T0" | gawk '{print $1-$2}'`
echo "Done in $DT seconds"
