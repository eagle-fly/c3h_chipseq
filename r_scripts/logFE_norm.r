# C3H chipseq project - normalization of logFE BedGraph data
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Mar 25, 2015

pos_minmax <- function(pos_data) {
  data_max <- max(pos_data)
  pos_data_new <- pos_data / data_max 
}

neg_minmax <- function(neg_data) {
  data_min <- min(neg_data)
  neg_data_new <- neg_data / -data_min
}

infile <- "/Users/ying/Work/BCCRC/Projects/c3h_chipseq/data/results/b6.k4m3.filtered.se/b6.k4m3_logFE.bdg"
data <- read.table (infile, header=FALSE)
logFE_old <- data[,4]
logFE_min <- min(logFE_old)
logFE_max <- max(logFE_old)
pos_sites <- logFE_old >= 0
logFE_pos <- logFE_old[pos_sites]
logFE_pos_new <- pos_minmax(logFE_pos)
logFE_new <- rep(0, length(logFE_old))
logFE_new[pos_sites] <- logFE_pos_new
neg_sites <- logFE_old < 0
logFE_neg <- logFE_old[neg_sites]
logFE_neg_new <- neg_minmax(logFE_neg)
logFE_new[neg_sites] <- logFE_neg_new

data_new <- cbind(data[,1:3], logFE_new)
outfile <- "/Users/ying/Work/BCCRC/Projects/c3h_chipseq/data/results/b6.k4m3.filtered.se/b6.k4m3_logFE.norm.bdg"
write.table(data_new, outfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
