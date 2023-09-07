TwoModeCE <- function(W, 
                      d.1 = 1, d.2 = 2,
                      k.r = 4, k.c = 4, 
                      d.r = 8, d.c = 8) {
     library(igraph)
     
     # Initial matrices and vectors
     A <- W[ , colSums(W) != 0, drop = FALSE] #dropping isolates
     At <- t(A) #affiliation matrix transpose
     r <- nrow(A) # number of rows
     c <- ncol(A) # number of columns
     rn <- rownames(A)
     cn <- colnames(A)
     u.r <- matrix(1, nrow = r, ncol = 1) #unit row vector
     u.c <- matrix(1, nrow = c, ncol = 1) #unit column vector
     
     # Derived matrices and vectors
     AAt <- A %*% At #row projection
     AtA <- At %*% A #column projection
     N <- as.numeric(t(u.r) %*% A %*% u.c) #affiliation matrix size
     D.r <- diag(as.vector(A %*% u.c), r, r) #row diagonal degree matrix
     D.c <- diag(as.vector(At %*% u.r), c, c) #column diagonal degree matrix
     iDr <- solve(D.r) #inverse of row diagonal degree matrix
     iDc <- solve(D.c)  #inverse of column diagonal degree matrix
     S.r <- A %*% iDc %*% At #degree-normalized similarity matrix for rows
     S.c <- At %*% iDr %*% A #degree-normalized similarity matrix for columns
     
     Z.r <- matrix(0, r, r)
     Z.c <- matrix(0, c, c)
     B <- rbind(cbind(Z.r, A), cbind(At, Z.c))
     colnames(B) <- rownames(B)
     B.gr <- graph_from_adjacency_matrix(B)
     clos <- closeness(B.gr)
     iCr <- solve(diag(clos[1:r], r, r))
     a <- r + 1
     iCc <- solve(diag(clos[a:length(clos)], c, c))
     CE.r <- eigen(iCr %*% A %*% iCc %*% At)$vectors #row CE
     CE.c <- eigen(iCc %*% At %*% iCr %*% A)$vectors #column CE
     S.r <- A %*% iCc %*% At #degree-normalized similarity matrix for rows
     S.c <- At %*% iCr %*% A #degree-normalized similarity matrix for columns
     CE.r <- apply(CE.r, 2, as.numeric)
     CE.c <- apply(CE.c, 2, as.numeric)
     eigvec.S.r <- eigen(S.r)$vectors #eigendecomposition of row similarity matrix
     eigvec.S.c <- eigen(S.c)$vectors #eigendecomposition of column similarity matrix
     km.r <- hkmeans(data.frame(eigvec.S.r[, 1:d.r]), k.r)
     km.c <- hkmeans(data.frame(eigvec.S.c[, 1:d.c]), k.c)
     dat.r <- data.frame(cbind(CE.r, cluster = km.r$cluster))
     dat.c <- data.frame(cbind(CE.c, cluster = km.c$cluster))
     dat.r$cluster <- factor(dat.r$cluster)
     dat.c$cluster <- factor(dat.c$cluster)
     p.r <- ggscatter(dat.r, 
                      x = names(dat.r)[d.1], y = names(dat.r)[d.2], 
                      color = "cluster",
                      label = rn, repel = TRUE)
     p.c <- ggscatter(dat.c, 
                      x = names(dat.r)[d.1], y = names(dat.r)[d.2], 
                      color = "cluster",
                      label = cn, repel = TRUE)

     
     return(list(B = B, clos = clos, CE.r = CE.r, CE.c = CE.c, 
                 p.r = p.r, p.c = p.c, dat.r = dat.r))
}