nucmer -b 200 muteswan.fa blackswan.fa
delta-filter -1 out.delta > out.1delta 
mummerplot --png --large out.1delta

grep -v \# out.fplot | sed 's/^$/NA NA 0/' > ali.fplot
grep -v \# out.rplot | sed 's/^$/NA NA 0/' > ali.rplot

show-coords -TH out.1delta | sort -k5,5nr | awk '{z+=1;print $8"\t"$1"\t"$2"\t"$9"\t"$3"\t"$4"\tz="z""}' > links.tsv
