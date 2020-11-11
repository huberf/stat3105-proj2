source("bomb_plan.R")

# DATA TO PREDICT FROM
start_time <- c("2020-08-18T17:50:42.000Z")
start_longitude <- c(-114.003)
start_latitude <- c(46.88678)

# GENERATE PREDICTIONS
prediction <- easy_predict(start_time, start_longitude, start_latitude)

# NOW PRINT TO USER
cat("Day Index:", prediction$day_idx, "\n\n")
cat("PREDICTION 1:\n")
cat("Location (lat,long):", prediction$bomb1$latitude,
    prediction$bomb1$longitude, "\n")
cat("UTM Location (lat,long):", prediction$bomb1$utm_latitude,
    prediction$bomb1$utm_longitude, "\n")
cat("UTC Time:", as.numeric(prediction$bomb1$time), "\n\n")
cat("PREDICTION 2:\n")
cat("Location (lat,long):", prediction$bomb2$latitude,
    prediction$bomb2$longitude, "\n")
cat("UTM Location (lat,long):", prediction$bomb2$utm_latitude,
    prediction$bomb2$utm_longitude, "\n")
cat("UTC Time:", as.numeric(prediction$bomb2$time), "\n\n")
