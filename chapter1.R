# https://online.stat.psu.edu/stat510/lesson/1

# in time series analysis, data are not necessarily independent and not necessarily
# iid

# first look: trend, seasonality (-> calendar time has effect), long-run cycle, 
# sigma = const, outliers, abrupt breaks in mean or variance

# example: annual earthquake series
library(dplyr)
earthquakes <- read.table(
  file='earthquakes.csv',
  header=TRUE,
  sep=','
  )
earthquakes$Date <- as.Date(earthquakes$Date, format='%m/%d/%Y')
earthquakes$Year <- format(earthquakes$Date, format='%Y')
earthquakes %>%
  group_by(Year) %>%
  count() -> earthquakes_by_year
par(mfcol = c(1, 3))
plot(earthquakes_by_year$Year, earthquakes_by_year$n, type='o', pch=20)
abline(h=mean(earthquakes_by_year$n))

# features of the series:
#   - upward trend
#   - no clear seasonality
#   - maybe 1 outlier in 2011 (712 earthquakes)
#   - more erratic towards right side, perhaps variance is higher there

# the AR(1) model:
# x_t = delta + phi1 * x_tminus1 + w_t
# assumptions:
#   - wt --> iid N(0, sigma_w^2)
#   - wt independent of x

# lag-1 plot
n <- nrow(earthquakes_by_year)
x_t <- earthquakes_by_year$n[2:n]
x_tminus1 <- earthquakes_by_year$n[1:(n-1)]
plot(x_tminus1, x_t)
abline(lm(x_t~x_tminus1))

# assuming we can just use OLS on x_t ~ x_tminus1 (for different reasons this is 
# not correct, discussed below), we'd do this:
ols <- lm(x_t ~ x_tminus1)
summary(ols)

# OLS results 
# model: x_t = 141.5 + 0.672 * x_tminus1, F(1, 50) = 22.81, p<0.001
# residual standard error = 91.4
# slope is significant, so lag1 variable is useful
# R^2 = 0.29
# residual analysis:
plot(ols$fitted.values, ols$residuals)
abline(h=0)


