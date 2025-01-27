TwoModeSR <- function(A, k = 4, n.eig = 4) {
   #loading functions
   source(here("Functions", "SimRank.R"))
   source(here("Functions", "eigval.scatter.R"))
   source(here("Functions", "eigvec.scatter.R"))
   source(here("Functions", "corr.scatter.R"))
   
   # SimRank of matrix A
   nr <- nrow(A)
   nc <- ncol(A)
   rn <- rownames(A)
   cn <- colnames(A)
   sr.res <- SimRank(A)
   res.eigen.r <- eigen(sr.res$Sr) #eigenvectors/values of row similarity matrix
   res.eigen.c <- eigen(sr.res$Sc) #eigenvectors/values of column similarity matrix
   e.r <- res.eigen.r$vectors[, 2:nr] #row eigenvectors minus first
   e.c <- res.eigen.c$vectors[, 2:nc] #column eigenvectors minus first
   v.r <- res.eigen.r$values[2:nr] #row eigenvalues minus first
   v.c <- res.eigen.c$values[2:nc] #column eigenvalues minus first
   ev.r1 <- e.r[, 1] #main eigenvector of row similarity matrix
   ev.r2 <- e.r[, 2] #second main eigenvector of row similarity matrix
   ev.c1 <- e.c[, 1] #main eigenvector of column similarity matrix
   ev.c2 <- e.c[, 2] #second main eigenvector of column similarity matrix
   
   # Affiliation matrix plot (SimRank ordering)
   A.ord <- as.matrix(A[order(ev.r1), order(ev.c1*-1)])
   A.plot <- ggcorrplot(t(A.ord)) 
   A.plot <- A.plot + scale_x_discrete(position = "top")
   A.plot <- A.plot + theme(legend.position = "none",
                             axis.text.x = element_text(hjust = -0.2))
   
   # Eigenvalue plot of SimRank similarity matrix
   eig.dat.r <- data.frame(k = 2:nr, value = v.r) %>% 
      mutate(value = value / max(value))
   eig.dat.c <- data.frame(k = 2:nc, value = v.c) %>% 
      mutate(value = value / max(value))
   eigval.plot.r <- eigval.scatter(eig.dat.r)
   eigval.plot.c <- eigval.scatter(eig.dat.c)
   
   # Main eigenvector plot of SimRank similarity matrix
   k.r <- hkmeans(data.frame(ev.r1), k = 2)
   k.c <- hkmeans(data.frame(ev.c1), k = 2)
   eigvec.dat.r <- data.frame(rank = rank(ev.r1), 
                              value = ev.r1,
                              lab = rn,
                              cluster = as.factor(k.r$cluster))
   eigvec.dat.c <- data.frame(rank = rank(ev.c1), 
                              value = ev.c1,
                              lab = cn,
                              cluster = as.factor(k.c$cluster))
   eigvec.plot.r <- eigvec.scatter(eigvec.dat.r, s = 14)
   eigvec.plot.c <- eigvec.scatter(eigvec.dat.c, s = 14)
   
   #k-means similarity-based clustering
   km.r <- hkmeans(data.frame(e.r[, 1:n.eig]), n.eig)
   km.c <- hkmeans(data.frame(e.c[, 1:n.eig]), n.eig)
   
   # Correspondence Plot (SimRank)
   a <- rbind(e.r[, 1:2], e.c[, 1:2]*-1)
   colnames(a) <- c("d1", "d2")
   b <- c(rownames(A), colnames(A))
   c <- factor(c(km.r$cluster, km.c$cluster))
   d <- data.frame(a, cluster = c, lab = b)
   corr.plot <- corr.scatter(d)

return(list(
   ev.r1 = ev.r1,
   ev.c1 = ev.c1, 
   A.plot = A.plot,
   eigval.plot.c = eigval.plot.c,
   eigval.plot.r = eigval.plot.r,
   eigvec.plot.c = eigvec.plot.c,
   eigvec.plot.r = eigvec.plot.r,
   corr.plot = corr.plot)
   )
   }


