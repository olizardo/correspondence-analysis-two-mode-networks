corr.scatter <- function(x, s = 12) {
     p <- ggscatter(data.frame(x), 
                    x = "d1", 
                    y = "d2", 
                    point = FALSE,
                    color = "cluster",
                    palette = "uchicago",
                    font.label =  s,
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