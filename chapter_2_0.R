library(astsa)
library(mgcv)
library(simts)

data(globtemp, package = 'astsa')
globtemp = gts(
  globtemp, 
  start=1880, 
  freq=1, 
  unit_ts='C', 
  name_ts='global temperature deviations',
  data_name="Evolution of Global Temperatures"
  )

plot(globtemp)
