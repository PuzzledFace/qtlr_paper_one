library(tidyverse)
library(lubridate)
library(knitr)
library(qtlr)
source("/home/bceuser/kirkpatj/Papers/qltr/paper1/functions.R", echo=FALSE)


d <- tibble(X=seq(-3, 3, 0.05), 
            Y=dnorm(X, 0, 1),
            Area=cut(X, breaks=c(-Inf, seq(-2.5, 2.5, 0.5), Inf)))

d %>% ggplot() +
        geom_area(aes(x=X, y=Y, fill=Area))
