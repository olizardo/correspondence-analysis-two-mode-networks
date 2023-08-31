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
        
        iDr <- matrix(0, nrow(A), nrow(A))
        iDc <- matrix(0, ncol(A), ncol(A))
        diag(iDr) <- 1/rowSums(A)
        diag(iDc) <- 1/colSums(A)
        
        iDs.r <- matrix(0, nrow(A), nrow(A))
        iDs.c <- matrix(0, ncol(A), ncol(A))
        diag(iDs.r) <- 1/sqrt(rowSums(A))
        diag(iDs.c) <- 1/sqrt(colSums(A))
        
        r <- nrow(A)
        c <- ncol(A)
        rn <- rownames(A)
        cn <- colnames(A)
        
        # Derived matrices and vectors
        CA.r <- svd(iDr %*% A %*% iDc %*% At) #row CA
        CA.c <- svd(iDc %*% At %*% iDr %*% A) #column CA
        CA.res <- svd(iDs.r %*% A %*% iDs.c) #row and column CA
        bon.res <- svd(A) #Bonacich eigenvector centrality
        
        S.r <- A %*% iDc %*% At #degree-normalized similarity matrix for rows
        S.c <- At %*% iDr %*% A #degree-normalized similarity matrix for columns
        
        Su.r <- A %*% At #raw similarity matrix for rows
        Su.c <- At %*% A #raw similarity matrix for columns
        
        S.r.eigval <- svd(S.r)$d #eigenvalues of normalized row similarity matrix
        S.c.eigval <- svd(S.c)$d #eigenvalues of normalized column similarity matrix
        
        eigvec.r <- as.matrix(CA.r$u[, 2:r]) #row coordinates
        eigvec.c <- as.matrix(CA.c$u[, 2:c]) #column coordinates
        
        bon.eigvec.r <- bon.res$u #row eigenvectors
        bon.eigvec.c <- bon.res$v #row eigenvectors
        
        #row and column names
        rownames(eigvec.r) <- rn
        rownames(eigvec.c) <- cn
        rownames(bon.eigvec.r) <- rn
        rownames(bon.eigvec.c) <- cn
        a <- r - 1
        b <- c - 1
        colnames(eigvec.r) <- paste("d", 1:a, sep = "")
        colnames(eigvec.c) <- paste("d", 1:b, sep = "")
        colnames(bon.eigvec.r) <- paste("e", 1:c, sep = "")
        colnames(bon.eigvec.c) <- paste("e", 1:c, sep = "")
        rownames(S.r) <- rn
        colnames(S.r) <- rn
        rownames(S.c) <- cn
        colnames(S.c) <- cn
        rownames(Su.r) <- rn
        colnames(Su.r) <- rn
        rownames(Su.c) <- cn
        colnames(Su.c) <- cn
        
        # K-means clustering of weighted similarity matrices
        S.r.dat <- data.frame(svd(S.r[, 1:d.r])$u)
        S.c.dat <- data.frame(svd(S.c[, 1:d.c])$u)
        rownames(S.r.dat) <- rn
        rownames(S.c.dat) <- cn
        km.r <- hkmeans(S.r.dat, k.r)
        km.c <- hkmeans(S.c.dat, k.c)
        km.dat.r <- data.frame(lab = rn, cluster = factor(km.r$cluster))
        km.dat.c <- data.frame(lab = cn, cluster = factor(km.c$cluster))
        
        # K-means clustering of unweighted similarity matrices
        Su.r.dat <- data.frame(svd(Su.r[, 1:d.r])$u)
        Su.c.dat <- data.frame(svd(Su.c[, 1:d.c])$u)
        rownames(Su.r.dat) <- rn
        rownames(Su.c.dat) <- cn
        km.r <- hkmeans(Su.r.dat, b.r)
        km.c <- hkmeans(Su.c.dat, b.c)
        bon.km.dat.r <- data.frame(lab = rn, cluster = factor(km.r$cluster))
        bon.km.dat.c <- data.frame(lab = cn, cluster = factor(km.c$cluster))
        
        # Main Eigenvector Plot
        eigvec.dat.r <- data.frame(rank = rank(eigvec.r[, 1]), 
                                   value = eigvec.r[, 1],
                                   lab = rn) %>% 
                left_join(km.dat.r)
        eigvec.dat.c <- data.frame(rank = rank(eigvec.c[, 1]), 
                                   value = eigvec.c[, 1],
                                   lab = cn)  %>% 
                left_join(km.dat.c)
        eigvec.plot.r <- eigvec.scatter(eigvec.dat.r)
        eigvec.plot.c <- eigvec.scatter(eigvec.dat.c)
        
        # Correspondence Plots (CA)
        corr.plot.r <- corr.scatter(eigvec.r)
        corr.plot.c <- corr.scatter(eigvec.c)
        eigvec.rc <- rbind(eigvec.r[, 1:2], eigvec.c[, 1:2])
        km.dat.rc <- rbind(km.dat.r, km.dat.c)
        eigvec.rc <- cbind(eigvec.rc, km.dat.rc)
        corr.plot <- corr.scatter(eigvec.rc)
        
        # Correspondence Plot (Eigenvector)
        bon.eigvec.rc <- rbind(bon.eigvec.r[, 1:2], bon.eigvec.c[, 1:2])
        bon.eigvec.rc[, 1] <- bon.eigvec.rc[, 1] * -1
        bon.km.dat.rc <- rbind(bon.km.dat.r, bon.km.dat.c)
        bon.eigvec.rc <- cbind(bon.eigvec.rc, bon.km.dat.rc)
        bon.corr.plot <- bon.corr.scatter(bon.eigvec.rc)
        
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
        
        # Eigenvalue Plots
        eig.dat.r <- data.frame(k = 1:nrow(A), value = S.r.eigval) %>% 
                mutate(value = value / max(value))
        eig.dat.c <- data.frame(k = 1:ncol(A), value = S.c.eigval) %>% 
                mutate(value = value / max(value))
        eigval.plot.r <- eigval.scatter(eig.dat.r) #row eigenvalue plot 
        eigval.plot.c <- eigval.scatter(eig.dat.c) #column eigenvalue plot
        
        # Similarity matrix plot (weighted by degree)
        norm.sim.r <- as.matrix(S.r/max(S.r))
        norm.sim.c <- as.matrix(S.c/max(S.c))
        norm.sim.r <- norm.sim.r[order(eigvec.r[, 1]), order(eigvec.r[, 1])]
        norm.sim.c <- norm.sim.c[order(eigvec.c[, 1]), order(eigvec.c[, 1])]
        sim.plot.r <- sim.plot(norm.sim.r) #row similarity plot
        sim.plot.c <- sim.plot(norm.sim.c) #column similarity plot
        
        # Similarity matrix plot (unweighted)
        norm.sim.r <- as.matrix(Su.r/max(Su.r))
        norm.sim.c <- as.matrix(Su.c/max(Su.c))
        norm.sim.r <- norm.sim.r[order(bon.eigvec.r[, 1]), order(bon.eigvec.r[, 1])]
        norm.sim.c <- norm.sim.c[order(bon.eigvec.c[, 1]), order(bon.eigvec.c[, 1])]
        bon.sim.plot.r <- sim.plot(norm.sim.r) #row similarity plot
        bon.sim.plot.c <- sim.plot(norm.sim.c) #column similarity plot
        
        return(list(
                A.plot.ca = A.plot1,
                A.plot.bon = A.plot2,
                eigval.plot.r = eigval.plot.r,
                eigval.plot.c = eigval.plot.c,
                eigvec.plot.r = eigvec.plot.r,
                eigvec.plot.c = eigvec.plot.c,
                corr.plot.r = corr.plot.r,
                corr.plot.c = corr.plot.c,
                corr.plot = corr.plot,
                bon.corr.plot = bon.corr.plot,
                sim.plot.r = sim.plot.r,
                sim.plot.c = sim.plot.c,
                bon.sim.plot.r = bon.sim.plot.r,
                bon.sim.plot.c = bon.sim.plot.c,
                km.r = km.r,
                km.c = km.c,
                eigvec.r = eigvec.r,
                eigvec.c = eigvec.c,
                eigval.r = CA.r$d,
                eigval.c = CA.c$d,
                CA.res = CA.res)
        )
        }