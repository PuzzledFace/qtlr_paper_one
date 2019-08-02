
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
  # mapping <- list("ARG"="ROW",  "BEL"="West", "BUL"="East", "CAN"="ROW",
  #                 "DEN"="West", "ESP"="West", "GBR"="West", "GER"="West",
  #                 "GRE"="West", "HUN"="East", "ITA"="West", "JPN"="ROW",
  #                 "LTU"="East", "MEX"="ROW", "NED"="West", "POL"="East",
  #                 "POR"="West", "ROM"="East", "SUI"="West", "USA"="ROW",
  #                 "UNK"="Unknown")
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

createQtlPlot <- function(mcmcData,
                          mcmcVar=Value,
                          groupData,
                          groupResponse=ObservedResponse,
                          groupSize=N,
                          groupType=NULL,
                          basicTheme=ggplot2::theme_light,
                          qtls=c(0.1, 0.2, 0.8, 0.9),
                          postScaleFactor=NULL,
                          mcmcAlpha=0.2,
                          nDensity=512)
{
  #Validate
  if (is.null(groupData)) stop("groupData cannot be null")
  if (is.null(quote(groupResponse))) stop("groupResponse cannot be null")
  if (!tibble::is_tibble(groupData)) stop("groupData must be a tibble")
  if (is.null(mcmcData)) stop("mcmcData cannot be null")
  if (!tibble::is_tibble(mcmcData)) stop("mcmcData must either be null or a tibble.")
  if (is.null(quote(mcmcVar))) stop("mcmcVar cannot be null")
  #Begin
  qResp <- dplyr::enquo(groupResponse)
  qSize <- dplyr::enquo(groupSize)
  if (!is.null(quote(groupType))) qGroup <- dplyr::enquo(groupType)
  qVar <- dplyr::enquo(mcmcVar)
  plot <- groupData %>% ggplot2::ggplot()
  #Plot the observed data
  if (is.null(qGroup)) {
    plot <- plot + ggplot2::geom_linerange(ggplot2::aes_q(x=qResp, ymin=0, ymax=qSize))
  } else {
    plot <- plot + ggplot2::geom_linerange(ggplot2::aes_q(x=qResp, ymin=0, ymax=qSize, linetype=qGroup))
  }
  #Plot the posterior density.  Could use geom_density directly, but shading
  #the action limits requires manipulation...
  if (!is.null(qtls) & length(qtls) > 0)
  {
    #Shading required
    if (is.null(postScaleFactor)) postScaleFactor <- 2 * (groupData %>%  dplyr::summarise(Max=max(!! qSize)))$Max[1]
    qtlTibble <- qtlFromQuantile(mcmcData, qtls)
    d <- ggplot2::ggplot_build(mcmcData %>% 
                                 ggplot2::ggplot() +
                                 ggplot2::geom_density(ggplot2::aes_q(qVar, y=~..scaled..*postScaleFactor), 
                                                       n=nDensity)
    )$data[[1]]
    d <- d %>% dplyr::mutate(Area=cut(x, 
                                      breaks=c(-Inf, qtlTibble$QTL, Inf),
                                      labels=1:(nrow(qtlTibble)+1)))
    plot <- plot +
      ggplot2::geom_line(data=d, ggplot2::aes(x=x, y=y))
    for (a in unique(d$Area))
      plot <- plot + ggplot2::geom_area(data=d %>% dplyr::filter(Area == a),
                                        ggplot2::aes(x=x, y=y, fill=Area),
                                        alpha=mcmcAlpha)
  } else {
    #No shading required
    plot <- plot + ggplot2::geom_density(data=mcmcData,
                                         ggplot2::aes_q(qVar, y=~..scaled..*30), 
                                         n=nDensity)
  }
  plot <- plot + 
    basicTheme() +
    ggplot2::theme(axis.ticks.y=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())
  return(plot)
}