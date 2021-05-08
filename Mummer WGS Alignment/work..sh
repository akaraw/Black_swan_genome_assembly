#Mummer3 version
#in the working directory

nucmer -b 200 muteswan.fa blackswan.fa
delta-filter -1 out.delta > out.1delta 
mummerplot --png --large out.1delta

grep -v \# out.fplot | sed 's/^$/NA NA 0/' > ali.fplot
grep -v \# out.rplot | sed 's/^$/NA NA 0/' > ali.rplot

show-coords -TH out.1delta | sort -k5,5nr | awk '{z+=1;print $8"\t"$1"\t"$2"\t"$9"\t"$3"\t"$4"\tz="z""}' > links.tsv

awk '{hue=sprintf("%03d", h); print "chr - "$1" "$1" 0 "$2" hue"hue; h+=40;}' muteswan.fasta.fai >  chr.kar
tac blackswan.fasta.fai | awk '{print "chr - "$1" "$1" 0 "$2" white"}' >> chr.kar
