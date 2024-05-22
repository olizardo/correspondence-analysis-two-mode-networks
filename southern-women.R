require(kableExtra)
require(here)
require(FactoMineR)

A <- matrix(
     c(
          1, 1, 1, 1, 1,  1, 0, 1, 1,  0, 0, 0, 0, 0, 
          1, 1, 1, 0, 1,  1, 1, 1, 0,  0, 0, 0, 0, 0, 
          0, 1, 1, 1, 1,  1, 1, 1, 1,  0, 0, 0, 0, 0, 
          1, 0, 1, 1, 1,  1, 1, 1, 0,  0, 0, 0, 0, 0, 
          0, 0, 1, 1, 1,  0, 1, 0, 0,  0, 0, 0, 0, 0, 
          0, 0, 1, 1, 1,  1, 0, 1, 0,  0, 0, 0, 0, 0, 
          0, 0, 0, 1, 1,  1, 1, 1, 0,  0, 0, 0, 0, 0, 
          0, 0, 0, 1, 1,  0, 1, 1, 1,  0, 0, 0, 0, 0, 
          
          0, 0, 0, 0, 0,  0, 1, 1, 1,  0, 0, 1, 0, 0, 
          0, 0, 0, 0, 0,  0, 0, 1, 1,  1, 0, 1, 0, 0, 
          0, 0, 0, 0, 0,  0, 0, 1, 1,  1, 0, 1, 1, 1,
          0, 0, 0, 0, 0,  0, 1, 1, 1,  1, 0, 1, 1, 1,
          0, 0, 0, 0, 0,  0, 1, 0, 1,  1, 1, 1, 1, 1,
          0, 0, 0, 0, 0,  0, 1, 1, 0,  1, 1, 1, 0, 0,
          0, 0, 0, 0, 0,  0, 0, 0, 1,  0, 1, 0, 0, 0,
          0, 0, 0, 0, 0,  0, 0, 0, 1,  0, 1, 0, 0, 0, 
          
          0, 0, 0, 0, 0,  1, 0, 1, 1,  0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,  0, 0, 1, 1,  0, 0, 0, 0, 0),
     ncol = 14, byrow = TRUE)
w <- c("EVELYN", "LAURA", "THERESA", "BRENDA", "CHARLOTTE", "FRANCES", "ELEANOR", "RUTH", "VERNE", "MYRNA", "KATHERINE", "SYLVIA", "NORA", "HELEN", "OLIVIA", "FLORA", "PEARL", "DOROTHY")
e <- paste("E", c(1:14), sep = "")
rownames(A) <- w
colnames(A) <- e

A <- A[sort(rownames(A)), ]
kbl(A) %>% 
     kable_minimal(font_size = 18,
                   full_width = FALSE,
                   bootstrap_options = c("condensed")) %>% 
     row_spec(1:18, extra_css = "border-bottom: 0px") %>% 
     save_kable(here("southern-women.html"), bs_theme = "simplex")


ca.res <- CA(A, graph = FALSE)
require(ggrepel)
plot.dat <- data.frame(rbind(ca.res$row$coord, ca.res$col$coord)[, 1:2])
plot.dat <- cbind(plot.dat, mode = c(rep(1, 18), rep(2, 14)), names = rownames(plot.dat))
p <- ggplot(plot.dat, aes(x = Dim.2, y = Dim.1*-1, 
                          color = factor(mode), label = names))
p <- p + geom_vline(xintercept = 0, color = "gray")
p <- p + geom_hline(yintercept = 0, color = "gray")
p <- p + geom_text_repel(size = 5, max.overlaps = 50)
p <- p + theme_minimal() 
p <- p + theme(legend.position = "none",
               axis.text = element_text(size = 20),
               axis.title = element_text(size = 24)
               )
p <- p + labs(x = "Second Dimension", y = "First Dimension")
p <- p + scale_color_manual(values = c('red', 'blue'))
png(file = here("Plots", "ca-southern-women.png"), width=600, height=625)
p
dev.off()

