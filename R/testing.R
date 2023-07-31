# library(tidyverse)
#
# # data preprocessing -----------------------------------------------------------
#
# unitzonal <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/unitzonal.rds")
# tsumdatp <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatp.rds")
# pltassign <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/pltassgn.rds")
#
# sampx <- tsumdatp %>%
#   left_join(pltassign, by = "CN") %>%
#   select(tcc16, elev, COUNTYCD) %>%
#   drop_na() %>%
#   mutate(COUNTYCD = as.character(COUNTYCD)) %>%
#   mutate(COUNTYFIPS = case_when(
#     str_length(COUNTYCD) < 2 ~ paste0("4100", COUNTYCD),
#     T ~ paste0("410", COUNTYCD)
#   ))
#
# domains <- sampx$COUNTYFIPS
#
# xsample <- sampx %>%
#   select(tcc16, elev)
#
# xpop <- unitzonal %>%
#   select(tcc16, elev, COUNTYFIPS, npixels) %>%
#   drop_na() %>%
#   rename(N = npixels)
#
# xpop_totals <- xpop %>%
#   mutate(tcc16 = tcc16*N, elev = elev*N) %>%
#   summarise(tcc16 = sum(tcc16), elev = sum(elev))
#
# y <- tsumdatp %>%
#   left_join(pltassign, by = "CN") %>%
#   select(tcc16, elev, COUNTYCD, BA_TPA_live_ADJ) %>%
#   drop_na() %>%
#   select(BA_TPA_live_ADJ) %>%
#   pull()
#
# #
# # ## greg testing
# #
modgreg_est <- modifiedGreg(
  y = y,
  xsample = xsample,
  xpop = xpop,
  domains = domains,
  datatype = "means",
  model = "linear",
  var_est = T,
  var_method = "LinHTSRS",
  messages = F
)


tsumdatc <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatc.rds")
tsumdatp <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatp.rds")

tsumdatp %>% summarise(plots = n_distinct(CN))
tsumdatc %>% summarise(plots = n_distinct(PLT_CN))

tsumdatc %>%
  group_by(PLT_CN) %>%
  sum

