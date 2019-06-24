library(tidyverse)
library(lubridate)

readData <- function(folder)
{
  ae <- read_sas(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", folder, "/ae.sas7bdat"))
  dm  <- read_sas(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", folder, "/dm.sas7bdat"))
  rv <- createSummaryData(ae, dm)
  rv <- rv %>% 
    add_column(Timepoint=paste0(substr(folder, 1, 4), "-", substr(folder, 5, 6)))
  return (rv)
}

mapCountryToRegion <- function(data)
{
  mapping <- list("ARG"="ROW",  "BEL"="West", "BUL"="East", "CAN"="ROW",
                  "DEN"="West", "ESP"="West", "GBR"="West", "GER"="West",
                  "GRE"="West", "HUN"="East", "ITA"="West", "JPN"="West",
                  "LTU"="East", "MEX"="West", "NED"="West", "POL"="East",
                  "POR"="West", "ROM"="East", "SUI"="West", "USA"="ROW")
  mapping <- mapping[unique(data$Country)]
  data <- data %>%
    mutate(Region=factor(Country,
                         labels=unlist(mapping, use.names=FALSE),
                         levels=sort(unique(Country))))
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

createIncrementalData <- function(base, end)
{
  aeBase <- read_sas(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", base, "/ae.sas7bdat"))
  aeBase <- aeBase %>%  select(STUDYID, USUBJID, AETERM, AESTDTC, AESPID)
  
  aeEnd <- read_sas(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", end, "/ae.sas7bdat"))
  ae <- aeBase %>% 
          full_join(aeEnd, by=c("STUDYID", "USUBJID", "AETERM", "AESTDTC")) %>% 
          filter(is.na(AESPID.x)) %>% 
          select(-AESPID.x, -AESPID.y)
  dm <- read_sas(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", end, "/dm.sas7bdat"))
  ae <- createSummaryData(ae, dm) %>%
          mutate(Baseline=base,
                 Endpoint=end)
  return(ae)
}

aeBase <- readData("201709")
mBase <- fitPoissonModel(aeBase$EventCount, aeBase$Exposure)
createQtlPlot(mBase %>% filter(Value < 50),
              groupData=aeBase,
              groupResponse=Rate,
              groupSize=Subjects,
              groupType=Region,
              qtls=c(0.05, 0.95)) + 
  labs(title="Pre September 2017")

aeInc <- createIncrementalData("201709","201809")
mInc <- fitPoissonModel(aeInc$EventCount, aeInc$Exposure)
createQtlPlot(mInc %>% filter(Value < 50),
              groupData=aeInc,
              groupResponse=Rate,
              groupSize=Subjects,
              groupType=Region,
              qtls=c(0.05, 0.95)) + 
  labs(title="Post September 2017")
