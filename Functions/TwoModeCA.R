TwoModeCA <- function(x, 
                      k.r = 2, 
                      k.c = 2, 
                      b.r = 2,
                      b.c = 2,
                      d.r = 2, 
                      d.c = 2) {
        library(factoextra)
        library(ggcorrplot)
        library(ggpubr)
        library(Matrix)

        # Plot Functions
        eigval.scatter <- function(x) {
                p <- ggscatter(x, 
                               x = "k", 
                               y = "value",
                               size = 2,
                               color = "red")
                p <- p + labs(y = "", x = "k") 
                p <- p + theme(axis.text.x = element_blank(),
                               axis.text.y = element_text(size = 8),
                               axis.line.y = element_blank(),
                               axis.line.x = element_blank())
                return(p)
                }
        eigvec.scatter <- function(x) {
                p <- ggscatter(x, 
                          y = "rank", 
                          x = "value", 
                          point = FALSE,
                          font.label =  12,
                          label = "lab", 
                          color = "cluster",
                          palette = "Dark2",
                          repel = TRUE
                          #ellipse = TRUE,
                          #ellipse.type = "convex",
                          #ellipse.alpha = 0.15,
                          #ellipse.border.remove = TRUE
                          )
                p <- p + theme(axis.text.y = element_blank(),
                               axis.line.x = element_blank(),
                               axis.line.y = element_blank(),
                               legend.position = "none")
                p <- p + labs(y = "Rank", x = "First Axis")
                return(p)
                }
        corr.scatter <- function(x) {
                p <- ggscatter(data.frame(x), 
                          x = "d1", 
                          y = "d2", 
                          point = FALSE,
                          color = "cluster",
                          palette = "uchicago",
                          font.label =  10,
                          label = rownames(x), 
                          repel = TRUE)
                p <- p + theme(legend.position = "none",
                               axis.line.x = element_blank(),
                               axis.line.y = element_blank())
                p <- p + geom_vline(aes(xintercept = 0), 
                                    color = "gray", linetype = 2)
                p <- p + geom_hline(aes(yintercept = 0), 
                                    color = "gray", linetype = 2)
                p <- p + labs(x = "First Axis", y = "Second Axis")
                return(p)
        }
        bon.corr.scatter <- function(x) {
                p <- ggscatter(data.frame(x), 
                               x = "e1", 
                               y = "e2", 
                               point = FALSE,
                               color = "cluster",
                               palette = "uchicago",
                               font.label =  10,
                               label = rownames(x), 
                               repel = TRUE)
                p <- p + theme(legend.position = "none",
                               axis.line.x = element_blank(),
                               axis.line.y = element_blank())
                p <- p + geom_vline(aes(xintercept = 0), 
                                    color = "gray", linetype = 2)
                p <- p + geom_hline(aes(yintercept = 0), 
                                    color = "gray", linetype = 2)
                p <- p + labs(x = "First Eigenvector", y = "Second Eigenvector")
                return(p)
        }
        
        sim.plot <- function(x) {
                p <- ggcorrplot(x)
                p <- p + scale_x_discrete(position = "top")
                p <- p + theme(legend.position = "none",
                               axis.text.x = element_text(hjust = -0.2)) 
                return(p)
        }
        
        # Initial matrices and vectors
        A <- x[ , colSums(x) != 0, drop = FALSE]
        A <- Matrix(A, sparse = TRUE)
        At <- t(A)
        r <- nrow(A) # number of rows
        c <- ncol(A) # number of columns
        rn <- rownames(A)
        cn <- colnames(A)
        u.r <- matrix(1, nrow = r, ncol = 1)
        u.c <- matrix(1, nrow = c, ncol = 1)

        
        # Derived matrices and vectors
        AAt <- A %*% At #row projection
        AtA <- At %*% A #column projection
        N <- as.numeric(t(u.r) %*% A %*% u.c)
        D.r <- diag(as.vector(A %*% u.c), r, r)
        D.c <- diag(as.vector(At %*% u.r), c, c)
        iDr <- solve(D.r)
        iDc <- solve(D.c)
        S.r <- A %*% iDc %*% At #degree-normalized similarity matrix for rows
        S.c <- At %*% iDr %*% A #degree-normalized similarity matrix for columns
        
        # Spectral matrices
        CA.r <- eigen(iDr %*% A %*% iDc %*% At) #row CA
        CA.c <- eigen(iDc %*% At %*% iDr %*% A) #column CA
        S.r.eig <- eigen(S.r) #eigendecomposition of row similarity matrix
        S.c.eig <- eigen(S.c) #eigendecomposition of column similarity matrix
        eigvec.r <- as.matrix(CA.r$vectors[, 2:r]) #row coordinates
        eigvec.c <- as.matrix(CA.c$vectors[, 2:c]) #column coordinates
        eigval.r <- as.matrix(CA.r$values)[2:r]
        eigval.c <- as.matrix(CA.c$values)[2:c]
        eigvec.S.r <- S.r.eig$vectors #row similarity eigenvectors
        eigvec.S.c <- S.c.eig$vectors #column similarity eigenvectors
        eigval.S.r <- S.r.eig$values #row similarity eigenvectors 
        eigval.S.c <- S.c.eig$values #column similarity eigenvectors
        bon.eigvec.r <- as.matrix(eigen(AAt)$vectors) #row eigenvectors
        bon.eigvec.c <- as.matrix(eigen(AtA)$vectors) #column eigenvectors
        bon.eigval.r <- as.matrix(eigen(AAt)$values)  #row eigenvalues
        bon.eigval.c <- as.matrix(eigen(AtA)$values)  #row eigenvalues
        
        # Normalizing row and column vectors
        norm.r <- function(x) {
                den <- t(x) %*% iDr %*% x
                x <- sqrt(N / den) %*% x
                }
        norm.c <- function(x) {
                x <- sqrt((N / (t(x) %*% iDc %*% x))) %*% x
                }
        #eigvec.r <- apply(eigvec.r, 2, norm.r)
        #eigvec.c <- apply(eigvec.c, 2, norm.c)
        #eigvec.r <- sqrt(eigval.r) * eigvec.r
        #eigvec.c <- sqrt(eigval.c) * eigvec.c
        
        #row and column names
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
        
        # Main Eigenvector Plot
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
        rn.a <- rownames(a)
        a <- apply(a, 2, as.numeric)
        b <- c(km.r$cluster, km.c$cluster)
        c <- data.frame(cbind(a, cluster = b))
        c$cluster <- factor(c$cluster)
        rownames(c) <- rn.a
        corr.plot <- corr.scatter(c)
        
        # Correspondence Plot (Bonacich)
        a <- rbind(bon.eigvec.r[, 1:2], bon.eigvec.c[, 1:2])
        rn.a <- rownames(a)
        a <- apply(a, 2, as.numeric)
        b <- c(bon.km.r$cluster, bon.km.c$cluster)
        c <- data.frame(cbind(a, cluster = b))
        c$cluster <- factor(c$cluster)
        rownames(c) <- rn.a
        bon.corr.plot <- bon.corr.scatter(c)
        
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

        # Similarity matrix plot (weighted by degree)
        norm.sim.r <- as.matrix(S.r/max(S.r))
        norm.sim.c <- as.matrix(S.c/max(S.c))
        norm.sim.r <- norm.sim.r[order(eigvec.r[, 1]), order(eigvec.r[, 1])]
        norm.sim.c <- norm.sim.c[order(eigvec.c[, 1]), order(eigvec.c[, 1])]
        sim.plot.r <- sim.plot(norm.sim.r) #row similarity plot
        sim.plot.c <- sim.plot(norm.sim.c) #column similarity plot
        
        # Similarity matrix plot (unweighted)
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