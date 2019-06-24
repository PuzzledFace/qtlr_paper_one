library(tidyverse)
library(knitr)
library(kableExtra)
library(lubridate)
library(qtlr)

maxPlotRate <- 20
interimMonth <- "09"
interimMonthName <- month.abb[as.integer(interimMonth)]
interimYear <- "2016"
qtls <- c(0.2, 0.8)

readData <- function(folder)
{
  load(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", folder, "/ae.Rdata"))
  load(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", folder, "/dm.Rdata"))
  dm <- dm %>% 
    mutate(COUNTRY=ifelse(str_length(COUNTRY) == 0, "UNK", COUNTRY))
  rv <- createSummaryData(ae, dm)
  rv <- rv %>% 
    add_column(Timepoint=paste0(substr(folder, 1, 4), "-", substr(folder, 5, 6)))
  return (rv)
}

mapCountryToRegion <- function(data)
{
  mapping <- list("ARG"="ROW",  "BEL"="West", "BUL"="East", "CAN"="ROW",
                  "DEN"="West", "ESP"="West", "GBR"="West", "GER"="West",
                  "GRE"="West", "HUN"="East", "ITA"="West", "JPN"="ROW",
                  "LTU"="East", "MEX"="ROW", "NED"="West", "POL"="East",
                  "POR"="West", "ROM"="East", "SUI"="West", "USA"="ROW",
                  "UNK"="Unknown")
  mapping <- mapping[unique(data$Country)]
  data <- data %>%
    mutate(Region=factor(Country,
                         labels=unlist(mapping, use.names=FALSE),
                         levels=unique(Country),
                         exclude=NULL))  #Changde from default of  NA, to allow for missing Country values
  return (data)
}

createSummaryData <- function(aeData, dmData)
{
  ae <- aeData %>% 
    #Create columns for Study, Centre and Subject 
    separate(USUBJID, c("Study", "Centre", "Subject")) %>% 
    #Convert character start dates to Dates
    mutate(AESTDT=ymd(AESTDTC, quiet=TRUE)) %>%
    #Handle incomplete start dates
    mutate(AESTDT=ifelse(is.na(AESTDT), ymd(paste0(AESTDTC, "-01"), quiet=TRUE), AESTDT)) %>% 
    mutate(AESTDT=ifelse(is.na(AESTDT), ymd(paste0(AESTDTC, "-01-01"), quiet=TRUE), AESTDT)) %>% 
    select(-STUDYID)
  dm <- dmData %>%
    separate(USUBJID, c("Study", "Centre", "Subject")) %>%
    mutate(RFENDT=ymd_hm(RFENDTC),
           RFSTDT=ymd_hm(RFSTDTC),
           Exposure=interval(RFSTDT, RFENDT)/dyears(1)) %>%
    replace_na(list(Exposure=0)) %>% 
    select(Study, Centre, Subject, Exposure, COUNTRY, RFENDT) %>% 
    rename(Country=COUNTRY)
  
  dm %>% 
    filter(Centre %in% c("283078","284108","288544")) %>% 
    kable(caption="Centres with dubious country attribution")
  
  #Filter to AEs with an onset no more than 30 days after the last dose of study
  #drug.
  ae <- ae %>% 
    left_join(dm, by=c("Study", "Centre", "Subject")) %>% 
    filter(AESTDT <= (RFENDT + days(30))) %>% 
    select(-RFENDT)
  
  aeCounts <- ae %>% 
    group_by(Study, Centre, Subject) %>% 
    summarise(EventCount=n()) %>% 
    ungroup()
  
  aeSummary <- dm %>% 
    left_join(aeCounts, by=c("Study", "Centre", "Subject")) %>% 
    replace_na(list(EventCount=0)) %>% 
    mutate(Country=ifelse(Country=="UNK" & Centre=="283078", "MEX", Country), 
           Country=ifelse(Country=="UNK" & Centre=="284108", "USA", Country),
           Country=ifelse(Country=="UNK" & Centre=="288544", "JPN", Country)) %>% 
    group_by(Study, Centre, Country) %>% 
    summarise(Subjects=n(),
              Exposure=sum(Exposure),
              EventCount=sum(EventCount),
              Rate=ifelse(Exposure > 0, EventCount/Exposure, 0)) %>% 
    ungroup() %>% 
    filter(Exposure > 0)
  aeSummary <- mapCountryToRegion(aeSummary)
  return (aeSummary)
}

aeBase <- readData(paste0(interimYear, interimMonth)) %>% 
  filter(Exposure > 0.2)

mBase <- fitPoissonModel(aeBase$EventCount, aeBase$Exposure)
fitIndex <-max(unique(mBase$Index), na.rm=TRUE)
mBasePlot <- mBase %>% filter(Parameter=="mu",
                              Index==fitIndex,
                              Value < maxPlotRate)

limits <- qtlToQuantile(mBasePlot, qtls) %>% add_column(Limit=c("Lower", "Upper"))

limits


aeBase %>% 
  ggplot() + 
    geom_point(aes(x=Exposure, y=Rate, size=Subjects, colour=Region)) + 
    geom_area(data=limits, aes(x=c(0, round(max(aeBase$Rate, 1))), y=rep(Quantile[1], 2)), inherit.aes=FALSE, colour="grey", alpha=0.25) + 
                theme_light()
