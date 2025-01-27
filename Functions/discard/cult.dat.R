cult.dat <- function() {
     library(dplyr)
     library(haven)
     taste.dat <- read_dta("C:/Users/Omar Lizardo/Google Drive/MISC DATA SOURCES/SSI-2012/SSI2012.dta")
     taste.dat <- taste.dat %>% 
          dplyr::select("id", ends_with(c("lis", "like")), -starts_with("none")) %>% 
          dplyr::select(c(1:41)) %>% 
          na.omit() %>% 
          mutate(Classical = classicallike * classicallis,
                 Opera = operalike * operalis,
                 Jazz = jazzlike * jazzlis,
                 Broadway = bwaystlike * bwaystlis,
                 Easy = moodezlike * moodezlis, 
                 Bigband = bbandlike * bbandlis,
                 Classic_Rock = croldlike * croldlis,
                 Country = countrylike * countrylis,
                 Bluegrass = blueglike * blueglis,
                 Folk = folklike * folklis,
                 Gospel = hymgoslike * hymgoslis,
                 Salsa = latlpsallike * latspsallis,
                 Rap_Hip_Hop = raphiphoplike * raphiphoplis,
                 Blues_RandB = blurblike * blurblis,
                 Reggae = reggaelike * reggaelis,
                 Pop = toppoplike * toppoplis,
                 Contemp_Rock = controcklike * controcklis,
                 Indie_Alt = indaltlike * indaltlis,
                 Dance_Club = danclublike * danclublis,
                 Metal = hvymtllike * hvymtllis
          ) %>%  #people are linked to genres that the both like and listen to
          dplyr::select(id, Classical:Metal) 
          A <- as.matrix(taste.dat[, 2:ncol(taste.dat)])
          rownames(A) <- taste.dat$id
          A <- A[-which(rowSums(A) == 0), ]
     return(A)
     }