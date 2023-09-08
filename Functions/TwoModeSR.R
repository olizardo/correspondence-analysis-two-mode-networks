TwoModeSR <- function(A, C = 0.8, iter = 10, k = 4) {
   nr <- nrow(A)
   nc <- ncol(A)
   dr <- rowSums(A)
   dc <- colSums(A)
   Sr <- diag(1, nr, nr)
   Sc <- diag(1, nc, nc)
   rn <- rownames(A)
   cn <- colnames(A)
   rownames(Sr) <- rn
   colnames(Sr) <- rn
   rownames(Sc) <- cn
   colnames(Sc) <- cn
   m <- 1
   while(m < iter) {
      Sr.pre <- Sr
      Sc.pre <- Sc
      for(i in 1:nr) {
         for(j in 1:nr) {
            if (i != j) {
               a <- names(which(A[i, ] == 1)) #objects chosen by i
               b <- names(which(A[j, ] == 1)) #objects chosen by j
               Scij <- 0
               for (k in a) {
                  for (l in b) {
                     Scij <- Scij + Sc[k, l] #i's similarity to j
                     }
                  }
                  Sr[i, j] <- C/(dr[i] * dr[j]) * Scij
               }
            }
         }
      for(i in 1:nc) {
         for(j in 1:nc) {
            if (i != j) {
               a <- names(which(A[, i] == 1)) #people who chose object i
               b <- names(which(A[, j] == 1)) #people who chose object j
               Srij <- 0
               for (k in a) {
                  for (l in b) {
                     Srij <- Srij + Sr[k, l] #i's similarity to j
                     }
                  }
                  Sc[i, j] <- C/(dc[i] * dc[j]) * Srij
               }
            }
         }
      m <- m + 1
   }

   corr.scatter <- function(x) {
      p <- ggscatter(data.frame(x), 
                     x = "d1", 
                     y = "d2", 
                     point = FALSE,
                     color = "cluster",
                     palette = "uchicago",
                     font.label =  13,
                     label = "lab", 
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
                     font.label =  9,
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
   
   e.r <- eigen(Sr) #eigenvectors/values of row similarity matrix
   e.c <- eigen(Sc) #eigenvectors/values of column similarity matrix
   ov.r <- e.r$vectors[, 2] #main eigenvector of row similarity matrix
   ov.c <- e.c$vectors[, 2] #main eigenvector of column similarity matrix
   
   # Affiliation matrix plot (SimRank ordering)
   A.ord <- as.matrix(A[order(ov.r), order(ov.c)])
   A.plot <- ggcorrplot(t(A.ord)) 
   A.plot <- A.plot + scale_x_discrete(position = "top")
   A.plot <- A.plot + theme(legend.position = "none",
                             axis.text.x = element_text(hjust = -0.2))
   
   # Correspondence plot of similarity matrix
   k.r <- hkmeans(data.frame(e.r$vectors[, 2:4]), k = 4)
   k.c <- hkmeans(data.frame(e.c$vectors[, 2:4]), k = 4)
   d.r <- data.frame(lab = rownames(A), 
                     d1 = round(as.numeric(e.r$vectors[, 2]), 4),
                     d2 = round(as.numeric(e.r$vectors[, 3]), 4),
                     cluster = as.factor(k.r$cluster)
                     )
   d.c <- data.frame(lab = colnames(A), 
                     d1 = round(as.numeric(e.c$vectors[, 2]), 4),
                     d2 = round(as.numeric(e.c$vectors[, 3]), 4),
                     cluster = as.factor(k.c$cluster)
                     )
   d.rc <- rbind(d.r, d.c)
   d.rc$d1 <- d.rc$d1 *-1
   d.rc$d2 <- d.rc$d2 *-1
   corr.plot <- corr.scatter(d.rc)
   eig.dat.r <- data.frame(k = 2:nr, value = e.r$values[2:nr]) %>% 
      mutate(value = value / max(value))
   eig.dat.c <- data.frame(k = 2:nc, value = e.c$values[2:nc]) %>% 
      mutate(value = value / max(value))
   eigval.plot.r <- eigval.scatter(eig.dat.r)
   eigval.plot.c <- eigval.scatter(eig.dat.c)
   
   # Main Eigenvector Plot
   k.r <- hkmeans(data.frame(e.r$vectors[, 2:2]), k = 2)
   k.c <- hkmeans(data.frame(e.c$vectors[, 2:2]), k = 2)
   
   eigvec.dat.r <- data.frame(rank = rank(ov.r*-1), 
                              value = ov.r*-1,
                              lab = rn,
                              cluster = as.factor(k.r$cluster))
   eigvec.dat.c <- data.frame(rank = rank(ov.c*-1), 
                              value = ov.c*-1,
                              lab = cn,
                              cluster = as.factor(k.c$cluster))
   eigvec.plot.r <- eigvec.scatter(eigvec.dat.r)
   eigvec.plot.c <- eigvec.scatter(eigvec.dat.c)
   
   return(list(Sr = round(Sr, 4), 
               Sc = round(Sc, 4),
               ov.r = ov.r,
               ov.c = ov.c,
               A.plot = A.plot,
               corr.plot = corr.plot,
               eigval.plot.c = eigval.plot.c,
               eigval.plot.r = eigval.plot.r,
               eigvec.plot.c = eigvec.plot.c,
               eigvec.plot.r = eigvec.plot.r
               ))
   }


