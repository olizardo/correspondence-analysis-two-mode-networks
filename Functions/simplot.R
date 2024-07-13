sim.plot <- function(x) {
     p <- ggcorrplot(x)
     p <- p + scale_x_discrete(position = "top")
     p <- p + theme(legend.position = "none",
                    axis.text.x = element_text(hjust = -0.2)) 
     return(p)
}