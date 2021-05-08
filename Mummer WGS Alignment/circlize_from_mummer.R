#Here we only selected the largest 35 chromsomes from the mute swan assemvly and the largest 36 Hi-C scaffolds from the black swan genome
library(circlize)
library(dplyr)
setwd('D:/05.OneDrive/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/')
kar <- read.table('chr.kar', header = F, sep = ' ')
head(kar)
kar_bs = kar[37:604,]
kar_ms = kar[1:35,]
head(kar_bs)
head(kar_ms, 36)
kar_bs$V6 <- as.numeric(kar_bs$V6)
kar_bs =kar_bs[order(kar_bs$V6, decreasing = T),]
kar_bs <- kar_bs[1:36,]
kar_bs
mykar <- rbind(kar_bs, kar_ms)
mykar[,4] <- gsub('HiC_scaffold_', 'bs_chr', mykar[,4])
mykar[,3] <- gsub('HiC_scaffold_', 'bs_chr', mykar[,3])
mykar

kar <- select(kar, 'V4', 'V5', 'V6', 'V7')
head(kar)
kar <- data.frame(lapply(kar, function(x){
  gsub(" ", "\t", x)
}))
head(kar)
tail(kar)

kar[,1] = gsub('HiC_scaffold_', 'bs_chr', kar[,1])
tail(kar)
colnames(kar) <- c('chr', 'start', 'end', 'col')
head(kar)

links <- read.table('links.tsv', header = F, sep = '\t')
head(links)

kar_bs = kar[37:604,]
kar_ms = kar[1:36,]
head(kar_bs)
head(kar_ms, 36)
kar_bs$end <- as.numeric(kar_bs$end)
kar_bs =kar_bs[order(kar_bs$end, decreasing = T),]
kar_bs <- kar_bs[1:36,]

kar
kar <- rbind(kar_bs, kar_ms)
kar <- kar[!kar$chr == 'chrW',]
head(kar)
tail(kar)


links
links[,4] <- gsub('HiC_scaffold_', 'bs_chr', links[,4])
links[,4]
(newlink <- links[links$V4%in%kar$chr,])
(newlink <- newlink[newlink$V1%in%kar$chr,])

mylinks <- links[links$V4%in%kar$chr,]
mylinks <- mylinks[mylinks$V1%in%kar$chr,]
dim(mylinks)

ms_links <- select(newlink, 'V1', 'V2', 'V3')
bs_links <- select(newlink, 'V4', 'V6', 'V5')

dim(bs_links)
dim(ms_links)
head(kar)

colnames(ms_links) <- c('chr', 'start', 'end')
colnames(bs_links) <- c('chr', 'start', 'end')

graphics.off()
svg('links_with_circlize_ms_bs.svg', width = 6, height = 5)

circos.clear()
circos.par(gap.after = c(rep(1, 34), 5, rep(1, 35), 5))
chromosome.index = c(paste0("chr", c(1:34, "Z")), 
                     rev(paste0("bs_chr", c(1:36))))

circos.initializeWithIdeogram(kar, plotType = NULL, chromosome.index = chromosome.index)


chromosome.index

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
highlight.chromosome(paste0("chr", c(1:34, "Z")), 
                     col = "red", track.index = 1)
highlight.chromosome(paste0("bs_chr", c(1:36)), 
                     col = "blue", track.index = 1)


#head(kar)
circos.track(ylim =c(0.2,0.1))
circos.genomicLink(ms_links, bs_links, col = rand_color(nrow(bs_links), transparency = 0.2), border = NA)
circos.clear()
dev.off()

#kar$chr

write.table(mykar, file = 'ms_bs.karyotype', sep = '\t', quote = F, row.names = F, col.names = F)
newlinks <- cbind(ms_links, bs_links)
#head(newlinks)

write.table(mylinks, file = 'ms_bs_link.tsv', sep = '\t', quote = F, row.names = F, col.names = F)

#dos2unix ms_bs.karyotype
#dos2unix ms_bs_link.tsv
#mv ms_bs.karyotype new.chr.kar
#mv ms_bs_link.tsv new.links.tsv
#perl addCol.pl new.chr.kar new.links.tsv (ref: https://bioinf.cc/misc/2020/08/08/circos-ribbons.html) 




