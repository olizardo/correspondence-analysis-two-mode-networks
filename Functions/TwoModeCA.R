TwoModeCA <- function(x, 
                      k.r = 2, 
                      k.c = 2, 
                      b.r = 2,
                      b.c = 2,
                      d.r = 2, 
                      d.c = 2) {
        
        #loading functions
        source(here("Functions", "eigval.scatter.R"))
        source(here("Functions", "eigvec.scatter.R"))
        source(here("Functions", "corr.scatter.R"))
        source(here("Functions", "simplot.R"))
        
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
        bon.eigvec.r[, 1] <- bon.eigvec.r[, 1] * -1 #positive values more centrality
        bon.eigvec.c[, 1] <- bon.eigvec.c[, 1] * -1 #positive values more centrality
        bon.eigvec.r[, 2] <- bon.eigvec.r[, 2] * -1 #positive values more centrality
        bon.eigvec.c[, 2] <- bon.eigvec.c[, 2] * -1 #positive values more centrality
        
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
        
        #k-means similarity-based clustering
        km.r <- hkmeans(data.frame(eigvec.S.r[, 1:d.r]), k.r)
        km.c <- hkmeans(data.frame(eigvec.S.c[, 1:d.c]), k.c)
        bon.km.r <- hkmeans(data.frame(bon.eigvec.r[, 1:d.r]), b.r)
        bon.km.c <- hkmeans(data.frame(bon.eigvec.c[, 1:d.c]), b.c)
        
        # Main Eigenvector Plot (CA)
        eigvec.dat.r <- data.frame(rank = rank(eigvec.r[, 1]), 
                                   value = as.numeric(eigvec.r[, 1]),
                                   lab = rn,
                                   cluster = factor(km.r$cluster))
        eigvec.dat.c <- data.frame(rank = rank(eigvec.c[, 1]), 
                                   value = as.numeric(eigvec.c[, 1]),
                                   lab = cn,
                                   cluster = factor(km.c$cluster))
        eigvec.plot.r <- eigvec.scatter(eigvec.dat.r)
        eigvec.plot.c <- eigvec.scatter(eigvec.dat.c)
        
        # Correspondence Plot (CA)
        a <- rbind(eigvec.r[, 1:2], eigvec.c[, 1:2])
        b <- rownames(a)
        a <- apply(a, 2, as.numeric)
        colnames(a) <- c("d1", "d2")
        c <- data.frame(a, 
                        cluster = factor(c(km.r$cluster, km.c$cluster)),
                        lab = b
                        )
        corr.plot <- corr.scatter(c)
        
        # Correspondence Plot (Bonacich)
        a <- rbind(bon.eigvec.r[, 1:2]*-1, bon.eigvec.c[, 1:2]*-1)
        b <- rownames(a)
        a <- apply(a, 2, as.numeric)
        colnames(a) <- c("d1", "d2")
        c <- data.frame(a, 
                        cluster = factor(c(bon.km.r$cluster, bon.km.c$cluster)),
                        lab = b
        )
        bon.corr.plot <- corr.scatter(c)
        
        # Affiliation matrix plot (CA ordering)
        A.ord <- as.matrix(A[order(eigvec.r[, 1]), order(eigvec.c[, 1])])
        A.plot <- ggcorrplot(t(A.ord)) 
        A.plot <- A.plot + scale_x_discrete(position = "top")
        A.plot1 <- A.plot + theme(legend.position = "none",
                                  axis.text.x = element_text(hjust = -0.2))
        
        # Affiliation matrix plot (Bonacich ordering)
        A.ord <- as.matrix(A[order(bon.eigvec.r[, 1]), order(bon.eigvec.c[, 1])])
        A.plot <- ggcorrplot(t(A.ord)) 
        A.plot <- A.plot + scale_x_discrete(position = "top")
        A.plot2 <- A.plot + theme(legend.position = "none",
                                  axis.text.x = element_text(hjust = -0.2)) 
        
        # CA Similarity matrix plot (normalized by matrix maximum value)
        norm.sim.r <- as.matrix(S.r/max(S.r))
        norm.sim.c <- as.matrix(S.c/max(S.c))
        norm.sim.r <- norm.sim.r[order(eigvec.r[, 1]), order(eigvec.r[, 1])]
        norm.sim.c <- norm.sim.c[order(eigvec.c[, 1]), order(eigvec.c[, 1])]
        sim.plot.r <- sim.plot(norm.sim.r) #row similarity plot
        sim.plot.c <- sim.plot(norm.sim.c) #column similarity plot
        
        # Bonacich Similarity matrix plot (normalized by matrix maximum value)
        norm.sim.r <- as.matrix(AAt/max(AAt))
        norm.sim.c <- as.matrix(AtA/max(AtA))
        norm.sim.r <- norm.sim.r[order(bon.eigvec.r[, 1]), order(bon.eigvec.r[, 1])]
        norm.sim.c <- norm.sim.c[order(bon.eigvec.c[, 1]), order(bon.eigvec.c[, 1])]
        bon.sim.plot.r <- sim.plot(norm.sim.r) #row similarity plot
        bon.sim.plot.c <- sim.plot(norm.sim.c) #column similarity plot
        
        # Eigenvalue Plots
        eig.dat.r <- data.frame(k = 1:nrow(A), value = as.numeric(eigval.S.r)) %>% 
                mutate(value = value / max(value))
        eig.dat.c <- data.frame(k = 1:ncol(A), value = as.numeric(eigval.S.c)) %>% 
                mutate(value = value / max(value))
        eigval.plot.r <- eigval.scatter(eig.dat.r) #row eigenvalue plot 
        eigval.plot.c <- eigval.scatter(eig.dat.c) #column eigenvalue plot
        
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
                bon.eigvec.c = bon.eigvec.c,
                eigval.plot.r = eigval.plot.r,
                eigval.plot.c = eigval.plot.c,
                ca.A = A.plot1,
                bon.A = A.plot2,
                eigvec.plot.r = eigvec.plot.r,
                eigvec.plot.c = eigvec.plot.c,
                ca.sim.plot.r = sim.plot.r,
                ca.sim.plot.c = sim.plot.c,
                bon.sim.plot.r = bon.sim.plot.r,
                bon.sim.plot.c = bon.sim.plot.c,                
                ca.corr.plot = corr.plot,
                bon.corr.plot  = bon.corr.plot))
        }