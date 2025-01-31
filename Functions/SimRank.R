SimRank <- function(A, C = 0.8, iter = 10) {
     nr <- nrow(A)
     nc <- ncol(A)
     dr <- rowSums(A)
     dc <- colSums(A)
     Sr <- diag(1, nr, nr)
     Sc <- diag(1, nc, nc)
     rn <- rownames(A)
     cn <- colnames(A)
     rownames(Sr) <- rn
     colnames(Sr) <- rn
     rownames(Sc) <- cn
     colnames(Sc) <- cn
     m <- 1
     while(m < iter) {
          Sr.pre <- Sr
          Sc.pre <- Sc
          for(i in 1:nr) {
               for(j in 1:nr) {
                    if (i != j) {
                         a <- names(which(A[i, ] == 1)) #objects chosen by i
                         b <- names(which(A[j, ] == 1)) #objects chosen by j
                         Scij <- 0
                         for (k in a) {
                              for (l in b) {
                                   Scij <- Scij + Sc[k, l] #i's similarity to j
                              }
                         }
                         Sr[i, j] <- C/(dr[i] * dr[j]) * Scij
                    }
               }
          }
          for(i in 1:nc) {
               for(j in 1:nc) {
                    if (i != j) {
                         a <- names(which(A[, i] == 1)) #people who chose object i
                         b <- names(which(A[, j] == 1)) #people who chose object j
                         Srij <- 0
                         for (k in a) {
                              for (l in b) {
                                   Srij <- Srij + Sr[k, l] #i's similarity to j
                              }
                         }
                         Sc[i, j] <- C/(dc[i] * dc[j]) * Srij
                    }
               }
          }
          m <- m + 1
     }
     return(list(Sr = Sr, Sc = Sc))
}