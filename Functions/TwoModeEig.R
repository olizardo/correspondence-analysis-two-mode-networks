TwoModeEig <- function(x, 
                      k.r = 2, 
                      k.c = 2, 
                      b.r = 2,
                      b.c = 2,
                      d.r = 2, 
                      d.c = 2) {
     # Initial matrices and vectors
     A <- x[ , colSums(x) != 0, drop = FALSE] #dropping isolates
     At <- t(A) #affiliation matrix transpose
     r <- nrow(A) #number of rows
     c <- ncol(A) #number of columns
     rn <- rownames(A) #row name vector
     cn <- colnames(A) #column name vector
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
     
     # Spectral matrices
     CA.r <- eigen(iDr %*% A %*% iDc %*% At) #row CA
     CA.c <- eigen(iDc %*% At %*% iDr %*% A) #column CA
     S.r.eig <- eigen(S.r) #eigendecomposition of row similarity matrix 
     S.c.eig <- eigen(S.c) #eigendecomposition of column similarity matrix 
     L.r.eig <- eigen(D.r - S.r) #eigendecomposition of row similarity matrix Laplacian
     L.c.eig <- eigen(D.c - S.c) #eigendecomposition of column similarity matrix Laplacian
     eigvec.r <- Re(CA.r$vectors[, 2:r]) #row coordinates
     eigvec.c <- Re(CA.c$vectors[, 2:c]) #column coordinates
     eigval.r <- CA.r$values[2:r] #row eigenvalues
     eigval.c <- CA.c$values[2:c] #column eigenvalues
     eigvec.S.r <- Re(S.r.eig$vectors) #row similarity eigenvectors
     eigvec.S.c <- Re(S.c.eig$vectors) #column similarity eigenvectors
     eigval.S.r <- S.r.eig$values #row similarity eigenvalues 
     eigval.S.c <- S.c.eig$values #column similarity eigenvalues
     eigvec.L.r <- Re(L.r.eig$vectors) #row Laplacian eigenvectors
     eigvec.L.c <- Re(L.c.eig$vectors) #column Laplacian eigenvectors
     eigval.L.r <- Re(L.r.eig$values) #row Laplacian eigenvectors
     eigval.L.c <- Re(L.c.eig$values) #column Laplacian eigenvectors
     bon.eigvec.r <- as.matrix(eigen(AAt)$vectors) #row eigenvectors
     bon.eigvec.c <- as.matrix(eigen(AtA)$vectors) #column eigenvectors
     bon.eigval.r <- as.matrix(eigen(AAt)$values)  #row eigenvalues
     bon.eigval.c <- as.matrix(eigen(AtA)$values)  #row eigenvalues
     
     #naming objects
     rownames(eigvec.r) <- rn
     rownames(eigvec.c) <- cn
     rownames(bon.eigvec.r) <- rn
     rownames(bon.eigvec.c) <- cn
     a <- r - 1
     b <- c - 1
     colnames(eigvec.r) <- paste("d", 1:a, sep = "")
     colnames(eigvec.c) <- paste("d", 1:b, sep = "")
     colnames(bon.eigvec.r) <- paste("e", 1:r, sep = "")
     colnames(bon.eigvec.c) <- paste("e", 1:c, sep = "")
     rownames(S.r) <- rn
     colnames(S.r) <- rn
     rownames(S.c) <- cn
     colnames(S.c) <- cn

     return(list(
          ca.eigvec.r = eigvec.r,
          ca.eigvec.c = eigvec.c,
          ca.eigval.r = eigval.r,
          ca.eigval.c = eigval.c,
          sim.eigvec.r = eigvec.S.r,
          sim.eigvec.c = eigvec.S.c,
          eigvec.L.r = eigvec.L.r,
          eigvec.L.c = eigvec.L.c,
          eigval.L.r = eigval.L.r,
          eigval.L.c = eigval.L.c,
          D.r = D.r,
          D.c = D.c,
          bon.eigval.r = bon.eigval.r,
          bon.eigval.r = bon.eigval.c,
          bon.eigvec.r = bon.eigvec.r,
          bon.eigvec.c = bon.eigvec.c))
     }