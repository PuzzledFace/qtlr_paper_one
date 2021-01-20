library(tidyverse)
library(qtlr)

generateData <- function(k, n, r1, r2) {
  set.seed(401010)
  
  d <- tibble(
        Centre=1:k,
        CentreRate=ifelse(Centre == 1, r1, r2),
        N=centreSize + round(runif(k, -5, 5))
      ) %>% 
      group_by(Centre) %>% 
      expand(
        nesting(Centre, CentreRate),
        PatNo=1:N
      ) %>% 
      ungroup() %>% 
      mutate(
        Interval=rpois(nrow(.), 20)
      ) %>% 
      group_by(Centre) %>% 
      mutate(
        BaselineDayStudyTime=cumsum(Interval)
      ) %>% 
      ungroup() %>% 
      mutate(
        FirstMissedDose=rpois(group_size(.), CentreRate),
        FirstMissedDoseStudyTime=FirstMissedDose + BaselineDayStudyTime
      )
  return(d)
}

nCentres <- 15
centreSize <- 10
standardRate <- 60
badRate <- 30


d <- generateData(k=nCentres, n=20, badRate, standardRate)

d1 <- d %>% 
        filter(BaselineDayStudyTime <= 180) %>% 
        mutate(
          Censored=FirstMissedDoseStudyTime > 180,
          EventTime=ifelse(Censored, 180 - BaselineDayStudyTime, FirstMissedDose)
        )

d1


d1 %>% 
  group_by(Centre) %>% 
  summarise(
    N=n(),
    Censored=sum(Censored),
    TotalExposure=sum(EventTime),
    MeanExposure=TotalExposure/(N-Censored)
  )

# getwd()
# saveRDS(d1, "temp.Rds")

fitTteModel(d1$EventTime, d1$Censored, TRUE, d1$Centre)
