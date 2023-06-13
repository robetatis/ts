# https://online.stat.psu.edu/stat510/lesson/1

# in time series analysis, data are not necessarily independent and not necessarily
# iid

# first look: trend, seasonality (-> calendar time has effect), long-run cycle, 
# sigma = const, outliers, abrupt breaks in mean or variance

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

earthquakes_by_year


plot(earthquakes_by_year$Year, earthquakes_by_year$n)


trend <- 0
n <- 100
tt <- seq(from=0, to=100, length.out=100)
y <- 10 + trend*tt + rnorm(n, mean=0, sd=0.01)

ytrain <- y[1:50]
ytest <- y[51:length(y)]


# y(t) = phi(0) + phi(1)*y(t-1) + e(t)

yt <- ytrain[2:length(ytrain)]
ytm1 <- ytrain[1:(length(ytrain)-1)]
AR1 <- lm(yt~ytm1)
phi0 <- AR1$coefficients[1]
phi1 <- AR1$coefficients[2]
et <- AR1$residuals
ytpred_train <- AR1$fitted.values

ytm1pred <- ytest[1:(length(ytest)-1)]
ytpred <- phi0 + phi1*ytm1pred

par(mar=c(3, 3, 1, 1))
plot(tt[1:50], ytrain, type='o', pch=20, 
     xlim=c(tt[1], tt[length(tt)]), 
     ylim=c(min(y), max(y)))
lines(tt[2:length(ytrain)], ytpred_train, col='gray', lwd=2)
abline(h=phi0)

lines(tt[51:length(tt)], ytest, pch=20, col='red', type='o')
lines(tt[52:length(tt)], ytpred, col='pink', lwd=2)


