# K-means clustering of weighted similarity matrices
S.r.dat <- data.frame(eigen(S.r[, 1:d.r])$vectors)
S.c.dat <- data.frame(eigen(S.c[, 1:d.c])$vectors)
rownames(S.r.dat) <- rn
rownames(S.c.dat) <- cn
km.r <- hkmeans(S.r.dat, k.r)
km.c <- hkmeans(S.c.dat, k.c)
km.dat.r <- data.frame(lab = rn, cluster = factor(km.r$cluster))
km.dat.c <- data.frame(lab = cn, cluster = factor(km.c$cluster))

# K-means clustering of unweighted similarity matrices
Su.r.dat <- data.frame(eigen(P.r[, 1:d.r])$vectors)
Su.c.dat <- data.frame(eigen(P.c[, 1:d.c])$vectors)
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