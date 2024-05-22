tm.corr.dist <- function(x) {
     library(expm)
     r <- nrow(x)
     c <- ncol(x)
     r.c <- diag(r)
     c.c <- diag(c)
     r.m <- rowMeans(x)
     c.m <- colMeans(x)
     
     tm.sim <- function(x) {
          r <- nrow(x)
          c <- ncol(x)
          r.r <- x %*% t(x)
          c.c <- t(x) %*% x
          r.s <- diag(r)
          c.s <- diag(c)
          
          sim <- function(a, b, c) {
               res <- a / sqrt(b*c) #cosine similarity 
               #res <- a (b + c - a) #Jaccard similarity
               #res <- (2*a) / (b + c) #Dice similarity
               return(res)
          }
          
          for (i in 1: r) {
               for (j in 1:r) {
                    if (i != j) {
                         r.s[i, j] <- sim(r.r[i, j], r.r[i, i], r.r[j, j])
                    }
               }
          }
          for (i in 1: c) {
               for (j in 1:c) {
                    if (i != j) {
                         c.s[i, j] <- sim(c.c[i, j], c.c[i, i], c.c[j, j])
                    }
               }
          }
          rownames(r.s) <- rownames(x)
          colnames(r.s) <- rownames(x)
          rownames(c.s) <- colnames(x)
          colnames(c.s) <- colnames(x)
          return(list(row.sims = r.s, col.sims = c.s))
          }
     
     sim.res <- tm.sim(x)
     r.s <- sim.res[[1]]
     c.s <- sim.res[[2]]
     
     for (i in 1: r) {
          for (j in 1:r) {
               if (i != j) {
                    r.x <- x[i, ] - r.m[i]
                    r.y <- x[j, ] - r.m[j]
                    r.xy <- r.x %*% c.s * t(r.y)
                    r.xx <- r.x %*% c.s * t(r.x)
                    r.yy <- r.y %*% c.s * t(r.y) 
                    r.num <- sum(r.xy)
                    r.den <- sqrt(sum(r.xx)) * sqrt(sum(r.yy))
                    r.c[i, j] <- round(r.num / r.den, 2)
                    }
               }
          }
     
     for (i in 1: c) {
          for (j in 1:c) {
               if (i != j) {
                    c.x <- x[, i] - c.m[i]
                    c.y <- x[, j] - c.m[j]
                    c.xy <- c.x %*% r.s * t(c.y) 
                    c.xx <- c.x %*% r.s * t(c.x) 
                    c.yy <- c.y %*% r.s * t(c.y) 
                    c.num <- sum(c.xy)
                    c.den <- sqrt(sum(c.xx)) * sqrt(sum(c.yy))
                    c.c[i, j] <- round(c.num / c.den, 2)
                    }
               }
          }
     
     rownames(r.c) <- rownames(x)
     colnames(r.c) <- rownames(x)
     rownames(c.c) <- colnames(x)
     colnames(c.c) <- colnames(x)
     return(list(row.sims = r.c, col.sims = c.c))
     }