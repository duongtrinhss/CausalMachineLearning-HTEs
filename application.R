library(haven)
stws <- read_stata("data/VietnamSTWS2015_wrangled.dta")
dim(stws)

stws %>% group_by(female) %>% summarise_each(funs(mean, max, min, sd, n()), logWage)

# Make a data.frame containing summary statistics of interest
summ_stats <- fBasics::basicStats(stws)
summ_stats <- as.data.frame(t(summ_stats))
# Rename some of the columns for convenience
summ_stats <- summ_stats %>% select("Mean", "Stdev", "Minimum", "1. Quartile", "Median", "3. Quartile", "Maximum")
summ_stats <- summ_stats %>% rename('Lower quartile'= '1. Quartile', 'Upper quartile' ='3. Quartile')

summ_stats
