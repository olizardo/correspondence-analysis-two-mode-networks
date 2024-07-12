SimRank1 <- function(A, alpha = 0.8, iter = 10, mode = "in") {
     V <- nrow(A)
     if (mode == "in") {
          d <- colSums(A)
          }
     else if (mode == "out") {
          d <- rowSums(A)
          }
     S <- matrix(0, V, V)
     n <- rownames(A)
     rownames(S) <- n
     colnames(S) <- n
     m <- 1
     while(m < iter) {
          for(i in 1:V) {
               for(j in 1:V) {
                    # similarities based on successors
                    if (d[i] == 0 | d[j] == 0) {
                         S[i, j] <- 0
                         }
                    if (i == j & (d[i] != 0 | d[j] != 0)) {
                         S[i, j] <- 1
                         }
                    if (i != j & d[i] > 0 & d[j] > 0 & mode == "in") {
                         a <- which(A[, i] == 1) # predecessors of i
                         b <- which(A[, j] == 1) # predecessors of j
                         K <- 0
                         for (k in a) {
                              for (l in b) { #avg. similarity between ij in-neighbors
                                   K <- K + (A[k, i]/d[i] * S[k, l] * A[l, j]/d[j])
                              }
                         }
                    S[i, j] <- alpha * K 
                    }
                    else if (i != j & d[i] > 0 & d[j] > 0 & mode == "out") {
                         a <- which(A[i, ] == 1) # successors of i
                         b <- which(A[j, ] == 1) # successors of j
                         K <- 0
                         for (k in a) {
                              for (l in b) { #avg. similarity between ij out-neighbors
                                   K <- K + (A[i, k]/d[i] * S[k, l] * A[j, l]/d[j])
                              }
                         }
                    S[i, j] <- alpha * K 
                    }
               }
          }
     m <- m + 1
     }
     return(list(S = round(S, 3), m = m, mode = mode))
}




