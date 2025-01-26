bump.plot <- function(x, p.name = "plot.png", w = 600, h = 625, c, rev = 1, n.r = 0.01) {
     p <- ggplot(data = dplyr::filter(x, n.iter < c), 
                 mapping = aes(x = n.iter, y = ref*rev, color = mode))
     p <- p + geom_bump(linewidth = 1.1, smooth = 8)
     p <- p + geom_point(size = 5)
     p <- p + geom_text_repel(data = dplyr::filter(x, n.iter == 2), 
                              aes(x = n.iter, label = mode), 
                              size = 5, nudge_x = -3)
     p <- p + geom_text_repel(data = dplyr::filter(x, n.iter == c),
                              aes(x = n.iter, label = mode),
                              size = 5, nudge_x = n.r)
     p <- p + theme_minimal_grid(font_size = 14, line_size = 0)
     p <- p + scale_x_continuous(limits = c(-2, c + 2), breaks = seq(2, c-2, by = 2))
     p <- p + labs(x = "Reflection")
     p <- p + scale_y_reverse()
     p <- p + theme(legend.position = "none",
                    panel.grid.major = element_blank(),
                    axis.text.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(size = 16)
     )
     set.seed(145)
     c <- sample(cols25(25), length(unique(x$mode)), replace = FALSE)
     p <- p + scale_color_manual(values=as.vector(c))
     png(file = here("Plots", p.name), width = w, height = h)
     p
     dev.off()
return(p)
}