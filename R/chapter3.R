# https://online.stat.psu.edu/stat510/lesson/3

# ****************************************
# Non-seasonal ARIMA models
# ****************************************

# also called Box-Jenkins
# the "I" indicates differencing. No differencing -> ARMA. differencing of 2nd
# order is difference of difference, etc.
# model specification is given through the order of each component -> (p, i, q):
#   - AR order: p
#   - MA order: q
#   - I order: i
#   - example: ARMA(1, 0 1) is AR(1) with no differencing and MA(1)

# identifying the model:
# explore three elements: time series, ACF and PACF

# time series plot:
#   - trend: 
#       - upward or downward linear -> 1st differences
#       - quadratic -> 2nd differences
#       - more complex -> smoothing
#   - seasonality -> seasonal differences
#   - outliers
#   - constant variance -> curved upward trend + increasing variance: log-transform

# ACF and PACF
#   - must be explored TOGETHER!!!
#   - some possible cases:
#                 PACF                   | ACF
#                 -----------------------|---------------------------------
#     - AR(1):    spike at lag 1         | tapers as rho_k = phi_1^k
#     - AR(2):    spikes at lags 1, 2    | sinusoidal, tapering
#     - MA(q):    tapering               | spikes at MA terms
#     - ARMA(p,q) tails off to 0         | tails off to 0
#

# other cases:
#   - if PACF and AF both have significant rho's over many lags, i.e., no tapering,
#     no tailing off to zero -> series may not be stationary -> do differencing!
#   - if ACF and PACF have no significant spikes anywhere -> random noise: 
#     x_t = w_t, w_t -> iid(0, sigma_w^2))
#   - if 1st differences are needed due to trend, but ACF and PACF have no spikes
#     at any lag for the differenced series -> random walk:
#     x_t = x_tminus1 + w_t, w_t -> iid(0, sigma_w^2))


ar1 <- arima.sim(model=list(ar=c(0.7)), n=100, sd=0.2)
ar2 <- arima.sim(model=list(ar=c(1.3, -0.8)), n=100, sd=0.2)
ma <- arima.sim(model=list(ma=c(0.7)), n=100, sd=0.2)
arma <- arima.sim(model=list(ma=c(0.5), ar=c(0.83)), n=100, sd=0.2)

dev.new()
par(mfcol=c(3, 4))

plot(ar1, type='o', pch=20, main='AR(1)')
acf(ar1, main='ACF AR(1)', ylim=c(-1, 1), xlim=c(1, 20))
pacf(ar1, main='PACF AR(1)', ylim=c(-1, 1), xlim=c(1, 20))

plot(ar2, type='o', pch=20, main='AR(2)')
acf(ar2, main='ACF AR(2)', ylim=c(-1, 1), xlim=c(1, 20))
pacf(ar2, main='PACF AR(2)', ylim=c(-1, 1), xlim=c(1, 20))

plot(ma, type='o', pch=20, main='MA(1)')
acf(ma, main='ACF MA(1)', ylim=c(-1, 1), xlim=c(1, 20))
pacf(ma, main='PACF MA(1)', ylim=c(-1, 1), xlim=c(1, 20))

plot(arma, type='o', pch=20, main='ARMA(1, 1)')
acf(arma, main='ACF ARMA(1,1)', ylim=c(-1, 1), xlim=c(1, 20))
pacf(arma, main='PACF ARMA(1,1)', ylim=c(-1, 1), xlim=c(1, 20))

dev.off()


# stochastic processes with no structure

# random noise
w_t <- runif(100, -1, 1)
x_t_random_noise = w_t

# random walk
x_t_random_walk <- numeric(length = 100)
x_t_random_walk[1] <- 0
for(i in 2:length(w_t)){
  x_t_random_walk[i] <- x_t_random_walk[i-1] + w_t[i]
  
}

# random walk with drift
delta <- 1.2
x_t_random_walk_drift <- numeric(length = 100)
x_t_random_walk_drift[1] <- 0
for(i in 2:length(w_t)){
  x_t_random_walk_drift[i] <- delta + x_t_random_walk_drift[i-1] + w_t[i]
  
}

windows(15, 10)
par(mfrow=c(3, 5), mar=c(4, 3, 3, 1))

plot(x_t_random_noise, type='o', pch=20, main='random noise')
acf(x_t_random_noise, xlim=c(1, 20), ylim=c(-1, 1), main='random noise ACF')
pacf(x_t_random_noise, xlim=c(1, 20), ylim=c(-1, 1), main='random noise PACF')
acf(diff(x_t_random_noise), xlim=c(1, 20), ylim=c(-1, 1), main='diff random noise ACF')
pacf(diff(x_t_random_noise), xlim=c(1, 20), ylim=c(-1, 1), main='diff random noise PACF')

plot(x_t_random_walk, type='o', pch=20, main='random walk')
acf(x_t_random_walk, xlim=c(1, 20), ylim=c(-1, 1), main='random walk ACF')
pacf(x_t_random_walk, xlim=c(1, 20), ylim=c(-1, 1), main='random walk PACF')
acf(diff(x_t_random_walk), xlim=c(1, 20), ylim=c(-1, 1), main='random walk diff ACF')
pacf(diff(x_t_random_walk), xlim=c(1, 20), ylim=c(-1, 1), main='random walk diff PACF')

plot(x_t_random_walk_drift, type='o', pch=20, main='random walk drift')
acf(x_t_random_walk_drift, xlim=c(1, 20), ylim=c(-1, 1), main='random walk drift ACF')
pacf(x_t_random_walk_drift, xlim=c(1, 20), ylim=c(-1, 1), main='random walk drift PACF')
acf(diff(x_t_random_walk_drift), xlim=c(1, 20), ylim=c(-1, 1), main='random walk drift diff ACF')
pacf(diff(x_t_random_walk_drift), xlim=c(1, 20), ylim=c(-1, 1), main='random walk drift diff PACF')

dev.off()


# ********************************
# model identification
# ********************************
library(forecast)
library(tseries)

# grab example series
x_t <- ts(scan('data/eriedata.dat'))

# check raw data, ACF and PACF
par(mfcol=c(1, 3))
plot(x_t, type='o', pch=20, main='data')
acf(x_t, xlim=c(1, 20), ylim=c(-1, 1), main='ACF')
pacf(x_t, xlim=c(1, 20), ylim=c(-1, 1), main='PACF')
# seems like AR(1) process

# stationarity test
tseries::kpss.test(x_t)
tseries::adf.test(x_t)

# fit model, check fit results
model_100 <- forecast::Arima(x_t, order=c(1, 0, 0))

# have package 'forecast' find best model automatically
model_auto <- auto.arima(X_train, seasonal=FALSE)

summary(model_100)
summary(model_auto)
# both models have same overall quality
# since 'forecast' implements AR(1) models as x_t - miu = phi_1*(x_tminus1 - miu) + w_t,
# the resulting model is: x_t - 14.23 = 0.617*(x_tminus1 - 14.23) + w_t 
# -> x_t = 5.45 + 0.617*x_tminus1 + w_t, with w_t -> iid(0, 1.722)

# significance of arima coeffs:
n <- length(X_train)
p <- 1
q <- 0
phi_1 <- as.numeric(model_100$coef[1])
tcrit <- qt(0.025, df=n-p-q-1)
tobs <- phi_1/0.1381
print(c(tcrit=tcrit, tobs=tobs))

# obs. vs pred, residuals acf and residuals vs. fitted
ei <- as.numeric(model_100$residuals)
yhat <- as.numeric(model_100$fitted)

par(mfcol=c(1, 3))
plot(x_t, main='x_t and fitted values')
points(yhat, col='red', pch=20)
acf(ei, xlim=c(1, 20), ylim=c(-1, 1), main='ACF residuals')
plot(yhat, ei, pch=20, main='Residuals vs. fitted')


# ********************************
# model diagnostics
# ********************************

# Ljung-Box test -> when we look at PACF or ACF, implicitly, we're doing multiple
# hypothesis testing. The Ljung-Box test is the way to avoid this. It's based on
# this statistic:
# Q(m) = n(n + 2) * sum_{j=1...m}(r_j^2/(n - j))
# -> it measures the 'cumulative' autocorrelation up to lag m
# this statistic can be used to check if the residuals' autocorrelations up to 
# lag m are collectively == 0
# A significant Q(m) suggests non-white-noise residuals (-> problem!)

# when used for checking significance of residuals' ACF, Q(m) has chi-squared
# distribution with df = m - (p + q + 1) -> i.e., Q(m) only defined for 
# m > (p + q + 1)
stats::Box.test(model_100$residuals, lag=20)

# when used for raw data, i.e., for checking 'collective significance' of 
# autocorrelation in data up to lag m, Q(m) has chi-squared distr. with df = m
stats::Box.test(x_t, lag=20)


# ********************************
# forecasting with ARIMA models
# ********************************








