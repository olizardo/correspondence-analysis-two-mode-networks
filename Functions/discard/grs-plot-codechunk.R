
source(here("Functions", "TwoModeGRS.R"))
tmSR.res <- TwoModeGRS(A)
a <- tmSR.res$eigvec.plot.r 
b <- tmSR.res$eigvec.plot.c 
c <- tmSR.res$eigval.plot.r
d <- tmSR.res$eigval.plot.c
a <- a + inset_element(c, left = 0.01, right = 0.4, bottom = 0.55, top = 1)
b <- b + inset_element(d, left = 0.01, right = 0.4, bottom = 0.55, top = 1)
e <- tmSR.res$A.plot

dat.r <- data.frame(ca = as.numeric(tm$ca.eigvec.r[,1])*-1,
                    sr = tmSR.res$ev.r1,
                    lab = rownames(A))
dat.c <- data.frame(ca = as.numeric(tm$ca.eigvec.c[,1])*-1,
                    sr = tmSR.res$ev.c1,
                    lab = colnames(A))
cor <- round(cor(dat.r[,c("ca", "sr")])[1, 2], 2)
g <- ggscatter(dat.r, x = "ca", y = "sr", font.label= 8, 
               label = "lab", repel = TRUE, 
               add = "reg.line", 
               add.params = list(color = "gray80", 
                                 linewidth = 0.5, linetype = 2)) +
     annotate("text", x = -.2, y = .3, 
              label = paste("r = ", cor), size = 6) 
cor <- round(cor(dat.c[,c("ca", "sr")])[1, 2], 2)
h <- ggscatter(dat.c, x = "ca", y = "sr", font.label= 8, 
               label = "lab", repel = TRUE, 
               add = "reg.line",
               add.params = list(color = "gray80", 
                                 linewidth = 0.5, linetype = 2)) +
     annotate("text", x = -.2, y = .3, 
              label = paste("r = ", cor), size = 6) 

png(file = here("Plots", "grs-plot-eigen-p.png"), width=600, height=625)
a
dev.off()
png(file = here("Plots", "grs-plot-eigen-g.png"), width=600, height=625)
b
dev.off()
png(file = here("Plots", "grs-plot-reord.png"), width=600, height=625)
e
dev.off()
png(file = here("Plots", "grs-corr-scatter.png"), width=350, height=625)
g / h
dev.off()