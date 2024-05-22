TwoModeGRS <- function(A, k = 4, n.eig = 4) {
   #loading functions
   source(here("Functions", "tm.corr.dist.R"))
   source(here("Functions", "eigval.scatter.R"))
   source(here("Functions", "eigvec.scatter.R"))
   source(here("Functions", "corr.scatter.R"))
   
   # GRS of matrix A
   nr <- nrow(A)
   nc <- ncol(A)
   rn <- rownames(A)
   cn <- colnames(A)
   sr.res <- tm.corr.dist(A)
   res.eigen.r <- eigen(sr.res$row.sims) #eigenvectors/values of row similarity matrix
   res.eigen.c <- eigen(sr.res$col.sims) #eigenvectors/values of column similarity matrix
   e.r <- res.eigen.r$vector #row eigenvectors
   e.c <- res.eigen.c$vectors #column eigenvectors 
   v.r <- res.eigen.r$values #row eigenvalues 
   v.c <- res.eigen.c$values #column eigenvalues
   ev.r1 <- e.r[, 1] #main eigenvector of row similarity matrix
   ev.r2 <- e.r[, 2] #second main eigenvector of row similarity matrix
   ev.c1 <- e.c[, 1] #main eigenvector of column similarity matrix
   ev.c2 <- e.c[, 2] #second main eigenvector of column similarity matrix
   
   # Affiliation matrix plot (GRS ordering)
   A.ord <- as.matrix(A[order(ev.r1), order(ev.c1)])
   A.plot <- ggcorrplot(t(A.ord)) 
   A.plot <- A.plot + scale_x_discrete(position = "top")
   A.plot <- A.plot + theme(legend.position = "none",
                            axis.text.x = element_text(hjust = -0.2))
   
   # Eigenvalue plot of GRS similarity matrix
   eig.dat.r <- data.frame(k = 1:nr, value = v.r) %>% 
      mutate(value = value / max(value))
   eig.dat.c <- data.frame(k = 1:nc, value = v.c) %>% 
      mutate(value = value / max(value))
   eigval.plot.r <- eigval.scatter(eig.dat.r)
   eigval.plot.c <- eigval.scatter(eig.dat.c)
   
   # Main eigenvector plot of GRS similarity matrix
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
   
   return(list(
      ev.r1 = ev.r1,
      ev.c1 = ev.c1,
      A.plot = A.plot,
      eigval.plot.c = eigval.plot.c,
      eigval.plot.r = eigval.plot.r,
      eigvec.plot.c = eigvec.plot.c,
      eigvec.plot.r = eigvec.plot.r))
}


