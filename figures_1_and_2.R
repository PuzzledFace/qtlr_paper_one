library(tidyverse)
library(qtlr)
library(cowplot)
library(knitr)
library(kableExtra)

set.seed(482214)
qtls <- c(0.1, 0.9)

centres <- tibble(Centre=1:15,
                  N=5+floor(20*runif(15)))

generateData <- function(rLabel, vLabel, a, b)
{
  d <- centres %>% 
         mutate(Prob=rbeta(15, a, b),
                Rate=rLabel,
                Var=vLabel,
                R=rbinom(15, N, Prob),
                ObservedRate=R/N)
  return (d)
}
data <- list()
data[[1]] <- generateData("Medium", "High",    2,  2)
data[[2]] <- generateData("Medium", "Medium", 10, 10)
data[[3]] <- generateData("Medium", "Low",    30, 30)
data[[4]] <- generateData("Low",    "High",    2,  6)
data[[5]] <- generateData("Low",    "Medium",  8, 24)
data[[6]] <- generateData("Low",    "Low",    21, 63)
data[[7]] <- generateData("High",   "High",    6,  2)
data[[8]] <- generateData("High",   "Medium", 24,  8)
data[[9]] <- generateData("High",   "Low",    63, 21)
data <- bind_rows(data)

nestedData <- data %>%  
          group_by(Rate, Var) %>%
          nest(.key=Observed)
results <- nestedData %>% 
  mutate(QTL=map(Observed, function(x) fitBinomialModel(n=x$N, r=x$R)),
         Limits=map(QTL, function(x) qtlToQuantile(x, qtls)))
                
tempFunc <- function(i)
{
  createQtlPlot(mcmcData=results$QTL[[i]] %>% filter(Index == 16), 
                groupData=results$Observed[[i]], 
                groupResponse=ObservedRate, qtls=qtls) +
    scale_fill_manual(values=c("grey", "white", "grey")) + 
    labs(x=" ") +
    guides(fill=FALSE)
}

plots <- lapply(1:9, tempFunc)

ggdraw() +
  draw_plot(plots[[9]], x=0.000, y=0.000, width=0.333, height=0.333) +
  draw_plot(plots[[8]], x=0.333, y=0.000, width=0.333, height=0.333) +
  draw_plot(plots[[7]], x=0.667, y=0.000, width=0.333, height=0.333) +
  draw_plot(plots[[3]], x=0.000, y=0.333, width=0.333, height=0.333) +
  draw_plot(plots[[2]], x=0.333, y=0.333, width=0.333, height=0.333) +
  draw_plot(plots[[1]], x=0.667, y=0.333, width=0.333, height=0.333) +
  draw_plot(plots[[6]], x=0.000, y=0.667, width=0.333, height=0.333) +
  draw_plot(plots[[5]], x=0.333, y=0.667, width=0.333, height=0.333) +
  draw_plot(plots[[4]], x=0.667, y=0.667, width=0.333, height=0.333)
