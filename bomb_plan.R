# GET DEPENDENCIES
library(tidyverse)
library(jsonlite)
library(ggplot2)
library(sp)
library(dlm)

# PREPARES TRAINING DATA
prepare_train_data <- function() {
	# Extract files
	gps_raw <- list()
	for (file in list.files('gps')) {
			gps_raw[[file]] <- fromJSON(paste('gps/', file, sep=''), simplifyVector=TRUE)
	}
  # now pull out the coordinates
	coords <- map(gps_raw, function (data) {t(simplify2array(data$features$geometry$coordinates))})
  # make data frame for each day
	days_dfs <- map(names(gps_raw), function(name) {
			data.frame(time =strptime(gps_raw[[name]]$features$properties$time,format = "%Y-%m-%dT%H:%M:%S", tz="UTC"),
								 longitude = coords[[name]][,1],
								 latitude = coords[[name]][,2])
	})
	names(days_dfs) <- names(gps_raw)

  # Process weather data
	weather_raw <- fromJSON('historical_weather.json', simplifyVector=TRUE)
	weather_obs <- weather_raw$observations
	temperature_all <- weather_obs[, c('valid_time_gmt', 'temp')]

  # PUT TOGETHER IN BIG DF
  # now heap coordinates together with flag for their day
	num_total <- sum(unlist(lapply(days_dfs, function(df) {return(dim(df)[1])})))
	num_total
	reduced_df <- data.frame(
			index=unlist(lapply(days_dfs, function(df) {return(1:length(df$time))})),
			time=unlist(lapply(days_dfs, function(df) {return(df$time)})),
			elapsed=unlist(lapply(days_dfs, function(df) {return(df$time - df$time[1])})),
			longitude=unlist(lapply(days_dfs, function(df) {return(df$longitude)})),
			latitude=unlist(lapply(days_dfs, function(df) {return(df$latitude)})),
			day=factor(unlist(lapply(1:length(days_dfs),
																			function(idx) {return(rep(idx, length(days_dfs[[idx]]$longitude)))}
																		 )))
	)

  # now get UTM coordinates
	spat_df <-SpatialPointsDataFrame(coords=reduced_df[,c("longitude", "latitude")],
																	 data=reduced_df["time"],
																	 proj4string=CRS("+proj=longlat +datum=WGS84 +units=m"))# This step converts the longitude/latitude -> UTM
	utm_df <-spTransform(spat_df, CRSobj = "+proj=utm +zone=12 +datum=WGS84 +units=m")
	utm_coords <- coordinates(utm_df)
	reduced_df$lat_m <- utm_coords[, "latitude"]
	reduced_df$long_m <- utm_coords[, "longitude"]

	# now add temperature
  # now pair data to temperature via linear interpolation
	indices <- unlist(lapply(reduced_df$time, function(time) {which(temperature_all[, 1] > time)[1]}))
	first_temp <- temperature_all[indices-1, 2]
	second_temp <- temperature_all[indices, 2]

  # now interpolate temperatures
	diff <- temperature_all[indices, 1] - temperature_all[indices-1, 1]
	weight <- (reduced_df$time - temperature_all[indices-1, 1])/diff
	final_temps <- first_temp*(1-weight) + second_temp*weight
  # now add temperature to dataframe
	reduced_df$temperature <- final_temps

	# now map "bombable" times
	diff_time <- reduced_df$time[2:length(reduced_df$time)] - reduced_df$time[1:(length(reduced_df$time)-1)]
	stationary <- c(TRUE, diff_time > 2*60)
	reduced_df$stationary <- stationary
	time_to_stationary <- unlist(map(reduced_df$time, function(time) {
			stat_diff <- time - reduced_df$time[reduced_df$stationary]
			ifelse(length(stat_diff[stat_diff >= 0]) != 0, min(stat_diff[stat_diff >= 0]), 0)
	}))
	reduced_df$bombable <- time_to_stationary > 5*60

  # now run Kalman smoother
	gps_variance <- 20^2
	v_mat <- diag(rep(gps_variance, 2))
	avg_walk_speed_m_per_sec <- 1.4  # https://en.wikipedia.org/wiki/Walking
  # now smooth per day
	reduced_df$smooth_lat_m <- NA
	reduced_df$smooth_long_m <- NA
	reduced_df$speed <- NA
	for (day in levels(factor(reduced_df$day))) {
			# estimate dt for every day
			data <- reduced_df[reduced_df$day == day, ]
			dt <- max(data[data$elapsed < 4000, ]$elapsed)/length(data[data$elapsed < 4000, ]$elapsed)
			g_mat <- matrix(c(1, 0, dt, 0,
												0, 1, 0, dt,
												0, 0, 1, 0,
												0, 0, 0, 1), byrow=TRUE, ncol=4)
			dlm_spec <- dlm(
				FF= matrix(c(1, 0, 0, 0,   0, 1, 0, 0), byrow=T, nrow=2),
				GG= g_mat,
				V = v_mat,
				W = diag(c(5, 5, 1, 1)^2),
				m0 = matrix(c(data$long_m[1], data$lat_m[1], rep(avg_walk_speed_m_per_sec / dt, 2)),
										ncol=1), # A vector by R defaults is a k by 1 matrix
				C0 = diag(rep(10^2, 4)))
			dlm_filter_mod <- dlmFilter(cbind(data$long_m, data$lat_m), dlm_spec)
			dlm_smooth_mod <- dlmSmooth(dlm_filter_mod)
			smoothed <- dlm_smooth_mod$s
			reduced_df$smooth_lat_m[reduced_df$day == day] <- smoothed[, 2][2:length(smoothed[, 2])]
			reduced_df$smooth_long_m[reduced_df$day == day] <- smoothed[, 1][2:length(smoothed[, 2])]
			reduced_df$speed[reduced_df$day == day] <- sqrt(dlm_smooth_mod$s[, 3]^2 + dlm_smooth_mod$s[, 4]^2)[2:length(smoothed[, 2])]
	}

  # now return the big dataframe
  return(reduced_df)
}



# MAIN PREDICTION FUNCTION
predict <- function(start_time, start_long, start_lat, train_data) {
    GRID_RESOLUTION = 2 # 2 meters
    #days <- c(1:8,11)
    days <- as.integer(levels(factor(train_data$day)))
    # now get train data slice
    lats <- seq(5193900, 5196700, by=GRID_RESOLUTION)
    train_data <- train_data[train_data$elapsed < 4000, ] # only keep data from first leg
    # only keep bombable points
    train_data <- train_data[train_data$bombable == TRUE, ] # don't train on points that can't be used
    DAYS <- 11
    # now cache certain views
    cache_longitudes <- list()
    cache_elapsed <- list()
    for (day in days) {
        cache_longitudes[[day]] <- train_data$long_m[train_data$day == day]
        cache_elapsed[[day]] <- train_data$elapsed[train_data$day == day]
    }
    # now get densities
    all_uncertainties <- unlist(map(1:length(lats), function(idx) {
        lat <- lats[idx]
        dists <- abs(train_data$lat_m - lat)
        day_idxes <- rep(NA, DAYS)
        latitude_diff <- rep(NA, DAYS)
        longitudes <- rep(NA, DAYS)
        elapsed <- rep(NA, DAYS)
        for (day in days) {
            day_idxes[day] <- which.min(dists[train_data$day == day])
            latitude_diff[day] <- dists[train_data$day == day][day_idxes[day]]
            longitudes[day] <- cache_longitudes[[day]][day_idxes[day]]
            elapsed[day] <- cache_elapsed[[day]][day_idxes[day]]
        }
        # Predict best longitudes and get their location uncertainty
        mean_long <- mean(longitudes, na.rm=TRUE)
        # go through every possible pairing and get the two with minimum uncertainty
        best_uncertainty <- Inf
        final_pred1 <- NA
        final_pred2 <- NA
        final_time1 <- NA
        final_time2 <- NA
        best_groupmask <- NA
        for (i in 1:(length(days)-1)) {
            for (j in (i+1):length(days)) {
                point1 <- longitudes[days[i]]
                point2 <- longitudes[days[j]]
                dists <- matrix(NA, ncol=2, nrow=DAYS)
                #print(dim(dists))
                #print(days)
                #print(length(longitudes))
                #print(length(abs(longitudes - point1)))
                #print(abs(longitudes - point1))
                dists[, 1] <- abs(longitudes - point1)
                dists[, 2] <- abs(longitudes - point2)
                parent <- apply(dists, 1, which.min)
                parent <- unlist(replace(parent, !sapply(parent, length), 0))
                #print(parent)
                # now get predicted point
                pred1 <- mean(longitudes[parent == 1])
                pred2 <- mean(longitudes[parent == 2])
                # now get uncertainty
                uncertainty1 <- sum(sqrt((longitudes[parent == 1] - pred1)^2 + latitude_diff[parent == 1]^2))
                uncertainty2 <- sum(sqrt((longitudes[parent == 2] - pred2)^2 + latitude_diff[parent == 2]^2))
                uncertainty <- (uncertainty1 + uncertainty2)/length(days)
                #uncertainty <- sum(total_uncertainty, na.rm=TRUE)/length(days)
                if (uncertainty < best_uncertainty) {
                    final_pred1 <- pred1
                    final_pred2 <- pred2
                    best_uncertainty <- uncertainty
                    best_groupmask <- parent
                }
            }
        }
        longitude_pred <- mean_long
        # Get position uncertainty
        #total_uncertainty <- sqrt((longitudes - mean_long)^2 + latitude_diff^2)
        #uncertainty <- sum(total_uncertainty, na.rm=TRUE)/length(days)
        # Now get time uncertainty
        final_time1 <- mean(elapsed[parent == 1], na.rm=TRUE)
        time_uncertainty <- sum(abs(elapsed[parent == 1] - final_time1), na.rm=TRUE)
        final_time2 <- mean(elapsed[parent == 2], na.rm=TRUE)
        time_uncertainty <- time_uncertainty +
                sum(abs(elapsed[parent == 2] - final_time2), na.rm=TRUE)
        time_uncertainty <- time_uncertainty/length(days)
        # Now return
        return(c(best_uncertainty, time_uncertainty, final_pred1, final_pred2, final_time1, final_time2))
    }))
    all_uncertainties <- matrix(all_uncertainties, ncol=6, byrow=TRUE)
    lat_uncertainties <- all_uncertainties[, 1]
    time_uncertainties <- all_uncertainties[, 2]
    # now pick the best points
    best_metric <- Inf
    best_lat_uncertainty <- NA
    best_time_uncertainty <- NA
    best_lat <- NA
    best_longs <- NA
    best_times <- NA
    for (i in 1:length(lats)) {
        metric <- time_uncertainties[i] + lat_uncertainties[i]*2
        if (metric < best_metric) {
            best_metric <- metric
            best_lat <- lats[i]
            best_longs <- all_uncertainties[i, 3:4]
            best_times <- all_uncertainties[i, 5:6]
            best_lat_uncertainty <- lat_uncertainties[i]
            best_time_uncertainty <- time_uncertainties[i]
        }
    }
    #cat("BEST:", best_lat_uncertainty, best_time_uncertainty, best_lat, best_longs, "\n")
    return(list(bomb1=list(utm_latitude=best_lat, utm_longitude=best_longs[1],
                  elapsed=best_times[1], time=start_time+best_times[1]),
             bomb2=list(utm_latitude=best_lat, utm_longitude=best_longs[2],
                  elapsed=best_times[2], time=start_time+best_times[2])))
}


## FUNCTIONS TO RUN PREDICTIONS
to_utm <- function(long, lat) {
    spat_df <- SpatialPoints(coords=cbind(long, lat),
                                     proj4string=CRS("+proj=longlat +datum=WGS84 +units=m"))# This step converts the longitude/latitude -> UTM
    utm_df <- spTransform(spat_df, CRSobj = "+proj=utm +zone=12 +datum=WGS84 +units=m")
    utm_coords <- coordinates(utm_df)
    return(list(long=utm_coords[1, "long"], lat=utm_coords[1, "lat"]))
}

easy_predict <- function(start_time, start_latitude, start_longitude) {
  train_data <- prepare_train_data()

  prediction <- NA
	for (day in 1:length(start_time)) {
    # convert provided time to local representation
		st_time <- as.numeric(strptime(start_time[day], format = "%Y-%m-%dT%H:%M:%S", tz="UTC"))
    # generate actual prediction
    prediction <- predict(st_time, start_latitude[day],
                          start_longitude[day], train_data)
	}
  prediction$day_idx <- 1
  return(prediction)
}
