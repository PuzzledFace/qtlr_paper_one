
readData <- function(folder)
{
  load(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", folder, "/ae.Rdata"))
  load(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", folder, "/dm.Rdata"))
  dm <- dm %>% 
    mutate(COUNTRY=ifelse(str_length(COUNTRY) == 0, "UNK", COUNTRY))
  rv <- createSummaryData(ae, dm)
  # rv <- rv %>% 
  #   add_column(Timepoint=paste0(substr(folder, 1, 4), "-", substr(folder, 5, 6)))
  return (rv)
}

mapCountryToRegion <- function(data)
{
  mapping <- list("ARG"="ROW",  "BEL"="ROW", "BUL"="East", "CAN"="ROW",
                  "DEN"="ROW", "ESP"="ROW", "GBR"="ROW", "GER"="ROW",
                  "GRE"="ROW", "HUN"="East", "ITA"="ROW", "JPN"="ROW",
                  "LTU"="East", "MEX"="ROW", "NED"="ROW", "POL"="East",
                  "POR"="ROW", "ROM"="East", "SUI"="ROW", "USA"="ROW",
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
  aeSummary <- aeSummary %>% mutate(Centre=1:nrow(aeSummary)) #Anonymise centres
  return (aeSummary)
}

createIncrementalData <- function(base, end)
{
  load(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", base, "/ae.Rdata"))
  aeBase <- ae %>%  select(STUDYID, USUBJID, AETERM, AESTDTC, AESPID)
  
  load(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", end, "/ae.Rdata"))
  aeEnd <- ae
  ae <- aeBase %>% 
    full_join(aeEnd, by=c("STUDYID", "USUBJID", "AETERM", "AESTDTC")) %>% 
    filter(is.na(AESPID.x)) %>% 
    select(-AESPID.x, -AESPID.y)
  load(paste0("/opt/bee/home_NEW/kirkpatj/R/QTL_test_data/WA29767/", end, "/dm.Rdata"))
  ae <- createSummaryData(ae, dm) %>%
    mutate(Baseline=base,
           Endpoint=end)
  return(ae)
}