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
   e.r <- res.eigen.r$vectors #row eigenvectors
   e.c <- res.eigen.c$vectors #column eigenvectors
   v.r <- res.eigen.r$values[2:nr] #row eigenvalues minus first
   v.c <- res.eigen.r$values[2:nc] #column eigenvalues minus first
   e.r <- e.r[, 2:nr] #dropping first row eigenvector
   e.c <- e.c[, 2:nc] #dropping first column eigenvector
   e.r <- e.r * -1 # reversing row eigevectors signs
   e.c <- e.c * -1 # reversing column eigevectors signs
   ov.r1 <- e.r[, 1] #main eigenvector of row similarity matrix
   ov.r2 <- e.r[, 2] #second main eigenvector of row similarity matrix
   ov.c1 <- e.c[, 1] #main eigenvector of column similarity matrix
   ov.c2 <- e.c[, 2] #second main eigenvector of column similarity matrix
   
   # Affiliation matrix plot (SimRank ordering)
   A.ord <- as.matrix(A[order(ov.r1), order(ov.c1)])
   A.plot <- ggcorrplot(t(A.ord)) 
   A.plot <- A.plot + scale_x_discrete(position = "top")
   A.plot <- A.plot + theme(legend.position = "none",
                             axis.text.x = element_text(hjust = -0.2))
   
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
   corr.plot <- corr.scatter(d.rc, s = 15)
   
   # Eigenvalue plot of SimRank similarity matrix
   eig.dat.r <- data.frame(k = 2:nr, value = v.r) %>% 
      mutate(value = value / max(value))
   eig.dat.c <- data.frame(k = 2:nc, value = v.c) %>% 
      mutate(value = value / max(value))
   eigval.plot.r <- eigval.scatter(eig.dat.r)
   eigval.plot.c <- eigval.scatter(eig.dat.c)
   
   # Main eigenvector plot of SimRank similarity matrix
   k.r <- hkmeans(data.frame(ov.r1), k = 2)
   k.c <- hkmeans(data.frame(ov.c1), k = 2)
   eigvec.dat.r <- data.frame(rank = rank(ov.r1), 
                              value = ov.r1,
                              lab = rn,
                              cluster = as.factor(k.r$cluster))
   eigvec.dat.c <- data.frame(rank = rank(ov.c1), 
                              value = ov.c1,
                              lab = cn,
                              cluster = as.factor(k.c$cluster))
   eigvec.plot.r <- eigvec.scatter(eigvec.dat.r)
   eigvec.plot.c <- eigvec.scatter(eigvec.dat.c)
   
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
               eigvec.plot.r = eigvec.plot.r
               )
          )
   }


