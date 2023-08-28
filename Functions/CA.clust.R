CA.clust <- function(x) {
     library(Matrix)
     library(factoextra)

     A <- x[ , colSums(x) != 0, drop = FALSE]
     A <- Matrix(A, sparse = TRUE)
     At <- t(A)
     iDr <- matrix(0, nrow(A), nrow(A))
     iDc <- matrix(0, ncol(A), ncol(A))
     diag(iDr) <- 1/rowSums(A)
     diag(iDc) <- 1/colSums(A)
     eig.r <- svd(iDr %*% A %*% iDc %*% At)$u[, 2] #row canonical correlation
     eig.c <- svd(iDc %*% At %*% iDr %*% A)$u[, 2] #column canonical correlation
     res.km.r <- hkmeans(data.frame(eig.r), 2, hc.method = "ward.D2") #clustering rows
     res.km.c <- hkmeans(data.frame(eig.c), 2, hc.method = "ward.D2") #clustering cols
     clus.dat.r <- data.frame(actors = rownames(x),
                              cluster = factor(res.km.r$cluster))
     clus.dat.c <- data.frame(events = colnames(x),
                              cluster = factor(res.km.c$cluster))   
     actors.1 <- dplyr::filter(clus.dat.r, cluster == 1)$actors
     actors.2 <- dplyr::filter(clus.dat.r, cluster == 2)$actors
     events.1 <- dplyr::filter(clus.dat.r, cluster == 1)$events
     events.2 <- dplyr::filter(clus.dat.r, cluster == 2)$events
     a.clusters$1 <- list(actors.1 = actors.1, actors.2 = actors.2)
     
    return(list(
         clus.dat.r = clus.dat.r,
         clus.dat.c = clus.dat.c,
         a.clusters = a.clusters)
         )
}