
u.r <- matrix(1, nrow = r, ncol = 1) #unit row vector
u.c <- matrix(1, nrow = c, ncol = 1) #unit column vector

# Derived matrices and vectors
AAt <- A %*% At #row projection
AtA <- At %*% A #column projection
N <- as.numeric(t(u.r) %*% A %*% u.c) #affiliation matrix size
D.r <- diag(as.vector(A %*% u.c)) #row diagonal degree matrix
D.c <- diag(as.vector(At %*% u.r)) #column diagonal degree matrix

iDr <- solve(D.r) #inverse of row diagonal degree matrix
iDc <- solve(D.c)  #inverse of column diagonal degree matrix
isDr <- solve(sqrt(D.r)) #inverse of square root of row diagonal degree matrix
isDc <- solve(sqrt(D.c))  #inverse of square root of column diagonal degree matrix
S.r <- A %*% iDc %*% At #similarity matrix for rows
S.c <- At %*% iDr %*% A #similarity matrix for columns
M.r <- iDr %*% S.r #centrality weighted row similarity matrix
M.c <- iDc %*% S.c #centrality weighted column similarity matrix
L.r <- D.r - S.r #row similarity matrix Laplacian
L.c <- D.c - S.c #column similarity matrix Laplacian
NL.r <- isDr %*% L.r %*% isDr #normalized row similarity matrix Laplacian
NL.c <- isDc %*% L.c %*% isDc #normalized row similarity matrix Laplacian

# Spectral matrices
CA.eig.r <- eigen(iDr %*% A %*% iDc %*% At) #eigendecomposition of centrality weighted row similarity matrix
CA.eig.c <- eigen(iDc %*% At %*% iDr %*% A) #eigendecomposition of centrality weighted column similarity matrix
L.eig.r <- eigen(L.r) #eigendecomposition of row similarity matrix Laplacian
L.eig.c <- eigen(L.c) #eigendecomposition of row similarity matrix Laplacian
NC.eig.r <- eigen(NL.r) #eigendecomposition of normalized row similarity matrix Laplacian
NC.eig.c <- eigen(NL.c) #eigendecomposition of normalized row similarity matrix Laplacian

# Extracting real eigenvectors and rounding
CA.vec.r <- round(Re(CA.eig.r$vectors), 4) 
CA.vec.c <- round(Re(CA.eig.c$vectors), 4) 
L.vec.r <- round(Re(L.eig.r$vectors), 4) 
L.vec.c <- round(Re(L.eig.c$vectors), 4) 
NC.vec.r <- round(Re(NC.eig.r$vectors), 4) 
NC.vec.c <- round(Re(NC.eig.c$vectors), 4) 

# Extracting real eigenvalues and rounding
CA.val.r <- round(Re(CA.eig.r$values), 4) 
CA.val.c <- round(Re(CA.eig.c$values), 4) 
L.val.r <- round(Re(L.eig.r$values), 4) 
L.val.c <- round(Re(L.eig.c$values), 4) 
NC.val.r <- round(Re(NC.eig.r$values), 4) 
NC.val.c <- round(Re(NC.eig.c$values), 4) 

name.mat.p <- function(x) {
     rownames(x) <- rn
     colnames(x) <- paste("ev", 1:r, sep = "")
     return(x)
}
name.mat.g <- function(x) {
     rownames(x) <- cn
     colnames(x) <- paste("ev", 1:c, sep = "")
     return(x)
}
return(list(S.r = S.r,
            S.c = S.c,
            M.r = name.mat.p(M.r), 
            M.c = name.mat.g(M.c), 
            CA.vec.r = name.mat.p(CA.vec.r), 
            CA.vec.c = name.mat.g(CA.vec.c), 
            L.vec.r = name.mat.p(L.vec.r), 
            L.vec.c = name.mat.g(L.vec.c), 
            NC.vec.r = name.mat.p(NC.vec.r), 
            NC.vec.r = name.mat.g(NC.vec.c),
            CA.val.r = CA.val.r, 
            CA.val.c = CA.val.c, 
            L.val.r = L.val.r, 
            L.val.c = L.val.c, 
            NC.val.r = NC.val.r, 
            NC.val.r = NC.val.c               
)
)

# Correspondence plot of SimRank similarity matrix
k.r <- hkmeans(data.frame(e.r[, 1:n.eig]), k = k)
k.c <- hkmeans(data.frame(e.c[, 1:n.eig]), k = k)
d.r <- data.frame(lab = rownames(A), 
                  d1 = round(ov.r1, 3),
                  d2 = round(ov.r2, 3),
                  cluster = as.factor(k.r$cluster)
)
d.c <- data.frame(lab = colnames(A), 
                  d1 = round(ov.c1, 3),
                  d2 = round(ov.c2, 3),
                  cluster = as.factor(k.c$cluster)
)
d.rc <- rbind(d.r, d.c)
corr.plot <- corr.scatter(d.rc, s = 14)





# SimRank similarity matrix plot (normalized by matrix maximum value)
Sr <- sr.res$Sr
Sc <- sr.res$Sc
diag(Sr) <- 0
diag(Sc) <- 0
norm.sim.r <- as.matrix(Sr/max(Sr))
norm.sim.c <- as.matrix(Sc/max(Sc))
norm.sim.r <- norm.sim.r[order(ov.r1), order(ov.r1)]
norm.sim.c <- norm.sim.c[order(ov.c1), order(ov.c1)]
sim.plot.r <- sim.plot(norm.sim.r) #row similarity plot
sim.plot.c <- sim.plot(norm.sim.c) #column similarity plot

return(list(Sr = round(sr.res$Sr, 3), 
            Sc = round(sr.res$Sc, 3),
            ov.r1 = ov.r1,
            ov.c1 = ov.c1,
            ov.r2 = ov.r2,
            ov.c2 = ov.c2,
            A.plot = A.plot,
            corr.plot = corr.plot,
            eigval.plot.c = eigval.plot.c,
            eigval.plot.r = eigval.plot.r,
            eigvec.plot.c = eigvec.plot.c,
            eigvec.plot.r = eigvec.plot.r,
            sim.plot.r = sim.plot.r,
            sim.plot.c = sim.plot.c
)
)