eigvec.scatter <- function(x, s = 10) {
     p <- ggscatter(x, 
                    y = "rank", 
                    x = "value", 
                    point = FALSE,
                    font.label =  s,
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