SimRank2 <- function(A, alpha = 0.8, iter = 10) {
     r <- nrow(A)
     c <- ncol(A)
     dr <- rowSums(A)
     dc <- colSums(A)
     Sr <- diag(1, r, r)
     Sc <- diag(1, c, c)
     rn <- rownames(A)
     cn <- colnames(A)
     rownames(Sr) <- rn
     colnames(Sr) <- rn
     rownames(Sc) <- cn
     colnames(Sc) <- cn
     m <- 1
     while(m < iter) {
          for(i in 1:r) {
               for(j in 1:r) {
                    if (i != j) {
                         a <- which(A[i, ] == 1) #objects chosen by i
                         b <- which(A[j, ] == 1) #objects chosen by j
                         Kij <- 0
                         for (k in a) {
                              for (l in b) { #i's similarity to j 
                                   Kij <- Kij + Sc[k, l] 
                              }
                         }
                         Sr[i, j] <- alpha/(dr[i] * dr[j]) * Kij
                    }
               }
          }
          for(i in 1:c) {
               for(j in 1:c) {
                    if (i != j) {
                         a <- which(A[, i] == 1) #people who chose object i
                         b <- which(A[, j] == 1) #people who chose object j
                         Klm <- 0
                         for (k in a) {
                              for (l in b) {
                                   Klm <- Klm + Sr[k, l] #i's similarity to j
                              }
                         }
                         Sc[i, j] <- alpha/(dc[i] * dc[j]) * Klm
                    }
               }
          }
          m <- m + 1
     }
     return(list(Sr = round(Sr, 3), Sc = round(Sc, 3)))
}