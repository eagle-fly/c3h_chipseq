# Calculate the checkerboard scores for each mutation-chimera pair
# Ying Zhang
# April 27, 2015

# define the c-score matrix
c_matrix <- matrix(rep(0, length(mutations) * length(chimeras)), nrow=length(mutations),ncol=length(chimeras))
colnames(c_matrix) <- chimeras
rownames(c_matrix) <- mutations

# calculate c-values for each Mutation-Chimera pair
for (i in 1:nrow(muta_matrix)) {
  for (j in 1:nrow(chim_matrix)) {
    mu <- as.vector(muta_matrix[i,])
    ch <- as.vector(chim_matrix[j,])
    ri <- sum(mu)
    rj <- sum(ch)
    sij <- sum(mu & ch)
    cij <- (ri - sij) + (rj - sij) / (ri + rj - sij)
    c_matrix[i,j] <- cij
  }
}
