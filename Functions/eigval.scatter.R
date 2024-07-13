eigval.scatter <- function(x, s = 2) {
     p <- ggscatter(x, 
                    x = "k", 
                    y = "value",
                    size = s,
                    color = "red")
     p <- p + labs(y = "", x = "k") 
     p <- p + theme(axis.text.x = element_blank(),
                    axis.text.y = element_text(size = 8),
                    axis.line.y = element_blank(),
                    axis.line.x = element_blank())
     return(p)
     }