# Monte Carlo simulation of possible presence/absence matrix
# Ying Zhang
# May 4, 2015

row_constraint <- c(2, 1, 3)
col_constraint <- c(3, 2, 2, 1)
row_count <- length(row_constraint)
col_count <- length(col_constraint)
m_matrix <- matrix(rep(0, row_count * col_count), nrow=row_count,ncol=col_count)

for (i in 1:row_count) {
  rand_idx <- sample.int(col_count, row_constraint[i])
  new_matrix <- m_matrix
  new_matrix[i, rand_idx] <- 1
  col_sum <- apply(new_matrix, 2, sum)
  col_chk <- which(col_sum - col_constraint > 0)
  
}