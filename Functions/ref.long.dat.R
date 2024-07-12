ref.long.dat <- function(a, b, max.iter) {
     library(tidyverse)
     pr.long <- data.frame(a)
     pr.long <- pr.long %>% 
          mutate(person = rownames(pr.long)) %>% 
          pivot_longer(
               cols = 1:ncol(a),
               names_to = "iter",
               values_to = "ref"
          ) %>%
       mutate(iter = factor(iter, ordered = TRUE, levels = names(pr.long))) %>%
       mutate(n.iter = as.integer(iter)*2)
     gr.long <- data.frame(b)
     gr.long <- gr.long %>% 
          mutate(group = rownames(gr.long)) %>% 
          pivot_longer(
               cols = 1:ncol(b),
               names_to = "iter",
               values_to = "ref"
          ) %>%
       mutate(iter = factor(iter, ordered = TRUE, levels = names(gr.long))) %>%
       mutate(n.iter = as.integer(iter)*2)
     return(list(person = pr.long, group = gr.long))
}