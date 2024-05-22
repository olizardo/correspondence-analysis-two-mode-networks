r <- c("Croatia", "Denmark", "Finland", "France", "Netherlands", "Serbia", "Spain", "Switzerland", "United Kingdom")
c <- c("Arts", "Institutional", "Values and trad.", "Geographical ID", "Group chars.", "Habits and lifestyles", "Human activity", "History and heritage", "Cultivation", "Good manners")

dat <- matrix(c(0.10,0.42,0.02,0.02,0.10,0.03,0.07,0.03,0.03,0.16, 0.16,0.22,0.10,0.06,0.07,0.06,0.18,0.06,0.04,0.05, 0.21,0.19,0.07,0.05,0.08,0.06,0.19,0.05,0.05,0.04, 0.13,0.25,0.03,0.04,0.02,0.03,0.02,0.05,0.42,0.01, 0.13,0.15,0.14,0.12,0.12,0.16,0.08,0.03,0.02,0.05, 0.10,0.27,0.06,0.03,0.08,0.04,0.05,0.03,0.03,0.31, 0.10,0.12,0.08,0.05,0.06,0.07,0.08,0.04,0.38,0.02, 0.15,0.23,0.07,0.06,0.06,0.05,0.11,0.06,0.16,0.05,0.17,0.07,0.23,0.10,0.18,0.11,0.05,0.06,0.02,0.02), nrow = 9, byrow = TRUE)
rownames(dat) <- r
colnames(dat) <- c

library(FactoMineR)
#res.CA <- CA(dat)

dat.f <- data.frame(country = r, dat) %>% 
     pivot_longer(
          cols = 2:11,
          names_to = "topics",
          values_to = "prop"
     ) 

library(ggplot2)
p <- ggplot(aes(y = country, x = prop, fill = country, group = country), data = dat.f)
p <- p + geom_bar(stat = "identity")
p <- p + facet_wrap(~ topics, ncol = 5)
p <- p + theme_minimal() + theme(legend.position = "none") + labs(y = "", x = "")
p