eigvec.plot <- function(x, s = 12, lab = "Second") {
     # Main Eigenvector Plot (CA)
     km.res <- hkmeans(x, k = 2)
     eigvec.dat <- data.frame(rank = rank(x), 
                           value = as.numeric(x),
                           lab = names(x),
                           cluster = factor(km.res$cluster))
     p <- ggscatter(eigvec.dat, 
                    y = "rank", 
                    x = "value", 
                    point = FALSE,
                    font.label =  s,
                    label = "lab", 
                    color = "cluster",
                    palette = "Dark2",
                    repel = TRUE
     )
     p <- p + theme(axis.text.y = element_blank(),
                    axis.line.x = element_blank(),
                    axis.line.y = element_blank(),
                    axis.title = element_text(size = 12),
                    axis.text.x = element_text(size = 12),
                    legend.position = "none")
     p <- p + labs(y = "Rank", x = paste(lab, "CA Axis", sep = " "))
 return(p)
}