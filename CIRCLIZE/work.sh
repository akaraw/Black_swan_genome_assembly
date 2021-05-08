#GC content of the genomes
bedtools makewindows -g ${genome}.fa.fai -w 25000 -s 12500 > ${genome}.genome.windows
bedtools nuc -fi ${genome}.fa -bed ${genome}.genome.windows |cut -f1,2,3,5 > ${genome}.genome.gc
