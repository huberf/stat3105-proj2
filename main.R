source("bomb_plan.R")

start_time <- c("2020-08-18T17:50:42.000Z")

start_longitude <- c(-114.003)

start_latitude <- c(46.88678)

prediction <- easy_predict(start_time, start_longitude, start_latitude)

cat(prediction)
