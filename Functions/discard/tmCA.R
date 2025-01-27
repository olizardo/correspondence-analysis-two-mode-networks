tmCA <- function(x) {
     # Initial matrices and vectors
     A <- x[ , colSums(x) != 0, drop = FALSE] #dropping isolates
     At <- t(A) #affiliation matrix transpose
     r <- nrow(A) #number of rows
     c <- ncol(A) #number of columns
     rn <- rownames(A) #row name vector
     cn <- colnames(A) #column name vector
     iDr <- diag(1/rowSums(A)) #inverse of row diagonal degree matrix
     iDc <- diag(1/colSums(A)) #inverse of column diagonal degree matrix
     # Spectral matrices
     eig.r <- eigen(iDr %*% A %*% iDc %*% At) #eigendecomposition of centrality weighted row similarity matrix
     eig.c <- eigen(iDc %*% At %*% iDr %*% A) #eigendecomposition of centrality weighted column similarity matrix
     eigvec.r <- round(Re(eig.r$vectors), 4) # Extracting real row eigenvectors and rounding
     eigvec.c <- round(Re(eig.c$vectors), 4) # Extracting real column eigenvectors and rounding
     eigval.r <- round(Re(eig.r$values), 4) # Extracting real row eigenvalues and rounding
     eigval.c <- round(Re(eig.c$values), 4) # Extracting real column eigenvalues and rounding
     rownames(eigvec.r) <- rn
     rownames(eigvec.c) <- cn
     return(list(eigvec.r = eigvec.r, 
                 eigvec.c = eigvec.c, 
                 eigval.r = eigval.r,
                 eigval.c = eigval.c)
     )
}