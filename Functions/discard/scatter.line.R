scatter.line <- function(c, r, p.col, t.col, ref.num, t.size = 3) {
     dat <- data.frame(c = c, r = r)
     p <- ggplot(data = dat, aes(x = c, y = r))
     p <- p + geom_point(size = 3, color = p.col) 
     p <- p + geom_smooth(method = "lm", se=FALSE, formula = y ~ x, color = "gray")
     p <- p + geom_text(aes(label = rownames(dat)), color = t.col, size = t.size)
     p <- p + theme_minimal() 
     p <- p + labs(x = "First CA Dimension", 
                   y = paste(ref.num, "Reflection (Standardized)", sep = " "))
     p <- p + theme(axis.text = element_text(size = 16),
                    axis.title = element_text(size = 16))
     return(p)
     }