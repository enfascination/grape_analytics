" use a linear model to predict veraison "
library(reshape2)
library(stringr)
library(lubridate)
library(data.table)

env <- Sys.getenv()
FLOWERING_ELSE_VERAISON <- as.logical(env["FLOWERING_ELSE_VERAISON"])
if (is.na(FLOWERING_ELSE_VERAISON)) { stop("is env loading?") }
TAVG_ELSE_MIDRANGE <- TRUE #False is max-min/2 + min .    ave works better for daily, midrange for hourly
if (TAVG_ELSE_MIDRANGE) {
    print("set for AVE TEMP")
} else {
    print("set for MID-RANGE TEMP")
}
if (FLOWERING_ELSE_VERAISON) {
    print("set for FLOWERING")
} else {
    print("set for VERAISON")
}

site_info <- read.csv("../data/01step/site_info.csv")
site_info <- site_info[,-1]

# wrangle GDD
gdd.orig <- read.csv("../data/01step/GDDValidation.csv")
gdd <- melt(gdd.orig, id=c("X"))
colnames(gdd) <- c("date", "identifier", "GDD")
gdd$identifier <- str_replace(gdd$identifier, '\\.', '-')
identifiers <- c( "AO-B1", "AO-B3", "AD-1", "AD-10", "SLV-2", "SLV-4", "CLV-17", "CLV-A1", "CLV-C8", "C-21", "C-2B", "C-1T", "DN-3A", "DN-5C", "DN-3B", "DN-3", "DAR-D3", "DAR-DE1", "DAR-E5", "DAR-I5", "DAR-B2", "DAR-2", "DAR-6", "HE-CS8", "J-3", "K-E2", "K-E2E", "LM-A2B", "LM-A3B", "LM-C8A", "MV-B2", "MV-E2E", "MV-E2W", "NS-3", "NS-7", "O-5H", "O-5I", "PM-3", "PM-4", "Q-DT", "Q-CN", "SE-A5", "SE-B2", "SE-D1E", "SE-D1W", "SE-E2B", "SS-1", "SS-5", "SH-B", "SH-LVL", "SH-SS", "SH-U7", "SO-1", "SO-5", "SO-10", "SO-35", "SW-18", "SW-7", "W-80L", "W-80H", "W-70L", "W-70H", "ST-50", "ST-52", "ST-26", "ST-33", "P-97", "P-98", "SHW-1", "SHW-3", "V-A6")
gdd <- subset(gdd, identifier %in% identifiers)
# imputing two ugly missing values roughly
gdd[is.na(gdd$GDD), "GDD"] <- 200

site_identifiers <- site_info$Identifier
for (i in 1:length(site_identifiers)) {
    gdds_by_id <- subset(gdd, identifier == site_identifiers[i] )
    if (nrow(gdds_by_id) == 0) { next }
	# pull the date of each stage transition at each site
    date_of_phase_change_bb <- date(ymd(site_info[site_info$Identifier == site_identifiers[i],"date_bb"]))
    date_of_phase_change_fl <- date(ymd(site_info[site_info$Identifier == site_identifiers[i],"date_fl"]))
    date_of_phase_change_vr <- date(ymd(site_info[site_info$Identifier == site_identifiers[i],"date_vr"]))
	# for each of the three major stages
	#   pull from gdd spreadsheet the right gdd (the one from the target day)
	#   then add to site_info for that site
	if(!is.na(date_of_phase_change_bb )){ ### this happens if no block observes both bb and veraison
		target_gdd_bb <- gdds_by_id[mdy(gdds_by_id$date) == date_of_phase_change_bb, "GDD"]
		site_info[site_info$Identifier == site_identifiers[i],"GDD_bb"] <- target_gdd_bb
	}
	if(!is.na(date_of_phase_change_fl ) | TRUE){ ### this happens if no block observes both bb and veraison
		# flowering is always observed
		target_gdd_fl <- gdds_by_id[mdy(gdds_by_id$date) == date_of_phase_change_fl, "GDD"]
		site_info[site_info$Identifier == site_identifiers[i],"GDD_fl"] <- target_gdd_fl
	}
	if(!is.na(date_of_phase_change_vr )){ ### this happens if no block observes both bb and veraison
		target_gdd_vr <- gdds_by_id[mdy(gdds_by_id$date) == date_of_phase_change_vr, "GDD"]
		site_info[site_info$Identifier == site_identifiers[i],"GDD_vr"] <- target_gdd_vr
	}
}
write.csv(site_info, "../data/01step/site_info_gdd.csv")

# wrangle action units as well, for comparing to GDD
action_units <- read.csv("../data/01step/actionunits.csv")
action_units <- action_units[,-1]

site_identifiers <- site_info$Identifier
for (i in 1:length(site_identifiers)) {
    units_by_id <- subset(action_units, identifier == site_identifiers[i] )
    if (nrow(units_by_id) == 0) { 
      print("ERROR: FKDSJL.  THis shouldn't happen.") 
      next }
	# pull the date of each stage transition at each site
    date_of_phase_change_bb <- date(ymd(site_info[site_info$Identifier == site_identifiers[i],"date_bb"]))
    date_of_phase_change_fl <- date(ymd(site_info[site_info$Identifier == site_identifiers[i],"date_fl"]))
    date_of_phase_change_vr <- date(ymd(site_info[site_info$Identifier == site_identifiers[i],"date_vr"]))
	# for each of the three major stages
	#   pull from units table the right units (the ones from the target day for each stage's accumulator)
	#   then add to site_info for that site
	if(is.na(date_of_phase_change_bb )){ ### this happens if no block observes both bb and veraison
    date_of_phase_change_bb <- date("2015-01-01")
  }
  target_units_bb <- sum(subset(units_by_id, date(date) <= date(date_of_phase_change_bb) & phase == "BB")$units)
  site_info[site_info$Identifier == site_identifiers[i],"units_bb"] <- target_units_bb
	
	if(!is.na(date_of_phase_change_fl ) | TRUE){ ### this happens if no block observes both bb and veraison
		# flowering is always observed
		target_units_fl <- sum(subset(units_by_id, (date(date) <= date(date_of_phase_change_fl)) & (date(date) >= date(date_of_phase_change_bb)) & phase == "FL")$units)
		target_units_gdd <- sum(subset(units_by_id, (date(date) <= date(date_of_phase_change_fl)) & phase == "GDD")$units)
		site_info[site_info$Identifier == site_identifiers[i],"units_fl"] <- target_units_fl
		site_info[site_info$Identifier == site_identifiers[i],"gdd_fl_alt"] <- target_units_gdd 
	}
	if(!is.na(date_of_phase_change_vr )){ ### this happens if no block observes both bb and veraison
		target_units_vr <- sum(subset(units_by_id, (date(date) <= date(date_of_phase_change_vr)) & (date(date) >= date(date_of_phase_change_fl)) & phase == "VR")$units)
		target_units_gdd <- sum(subset(units_by_id, (date(date) <= date(date_of_phase_change_vr)) & phase == "GDD")$units)
		site_info[site_info$Identifier == site_identifiers[i],"units_vr"] <- target_units_vr
		site_info[site_info$Identifier == site_identifiers[i],"gdd_vr_alt"] <- target_units_gdd
	}
}
write.csv(site_info, "../data/01step/site_info_gdd.csv")

### DAILY TEMP HISTOGRAM AS PREDICTORS
site_info <- read.csv("../data/01step/site_info.csv")
site_info <- site_info[,-1]
climate.data <- read.csv("../data/01step/climate.data.csv")
climate.data <- climate.data[,-1]
### init growth of site_info
### prelims to building histograms
tmin_global <- floor(min(as.integer(climate.data$TMIN)))
tmax_global <- ceiling(max(as.integer(climate.data$TMAX)))
t_names <- paste(tmin_global:tmax_global, 'C', sep='')
t_names <- str_replace(t_names, "-", "_")
## prep dummary data frame and merge it in
temp_hist_template <- data.frame(setNames(rep(list(0), length(tmin_global:tmax_global)), t_names), check.names=FALSE)
site_info <- cbind(site_info, temp_hist_template)
site_info 
### populate hists in site_info
site_siteids <- site_info$Site
for (i in 1:length(site_siteids)) {
	site_idx <- site_info$Site == site_siteids[i]
	if (FLOWERING_ELSE_VERAISON) {
    temp_by_site <- subset(climate.data, Site == site_siteids[i] & 
						                 ymd(date) >= ymd(site_info[site_idx,"date_bb"]) & 
										 ymd(date) <= ymd(site_info[site_idx,"date_fl"]) )
	} else {
    temp_by_site <- subset(climate.data, Site == site_siteids[i] & 
						                 ymd(date) >= ymd(site_info[site_idx,"date_fl"]) & 
										 ymd(date) <= ymd(site_info[site_idx,"date_vr"]) )
	}
	if (TAVG_ELSE_MIDRANGE) {
		specgd <- table(round(temp_by_site$TAV))
	} else {
		specgd <- table(round(( temp_by_site$TMAX - temp_by_site$TMIN ) / 2 + temp_by_site$TMIN))
	}
	new_names <- paste(names(specgd), 'C', sep='')
	site_info[site_idx, new_names] <- specgd
}
write.csv(site_info, "../data/01step/site_info_specgd.csv")

### HOURLY TEMP HISTOGRAM AS PREDICTORS
### TODO
if (FALSE) { ### re-read raw xlsx.  super slow: try to get fread working on excel
	library(readxl)
	col_types <- c(c("numeric", "text"), rep("numeric",15), c("text"), rep("numeric",4))
	hourly1 <- read_xlsx("../data/forrestel/ws_3-95_hourly_data.xlsx", col_types = col_types)
	hourly2 <- read_xlsx("../data/forrestel/ws_101-644_hourly_data.xlsx", col_types = col_types)
	hourly3 <- read_xlsx("../data/forrestel/ws_669-969_hourly_data.xlsx", col_types = col_types)
	ws_hourly <- rbind(hourly1,hourly2, hourly3)
	write.csv(ws_hourly, "../data/01step/ws_hourly.csv")
	names(hourly1)
	names(hourly2)
	names(hourly3)
	names(ws_hourly)
	dim(hourly1)
	dim(hourly2)
	dim(hourly3)
	dim(ws_hourly)
} else {
	ws_hourly <- data.frame(fread("../data/01step/ws_hourly.csv"))
	ws_hourly <- ws_hourly[,-1] 
}
### missing observations at one important weather station
ws_hourly <- ws_hourly[!is.na(ws_hourly[,"temperature_avg"]),]
str(ws_hourly)
unique(ws_hourly$weather_station_id)
unique(ws_hourly$weather_station_name)
summary(ws_hourly[,c("timestamp", "weather_station_id", "weather_station_name", "temperature_avg", "temperature_max", "temperature_min")])

# table I'll merge into
site_info <- read.csv("../data/01step/site_info.csv")
site_info <- site_info[,-1]
# get temp bounds for final histogram
tmin_global <- floor(min(ws_hourly$temperature_min))
tmax_global <- ceiling(max(ws_hourly$temperature_max))
t_names <- paste(tmin_global:tmax_global, 'C', sep='')
t_names <- str_replace(t_names, "-", "_")
## prep dummary data frame and merge it in
temp_hist_template <- data.frame(setNames(rep(list(0), length(tmin_global:tmax_global)), t_names), check.names=FALSE)
site_info <- cbind(site_info, temp_hist_template)
# this implementation is ineffciient: will chug the same ws twice
#   OK though, takes less than a minute
for (i in 1:nrow(site_info)) {
	site_ws <- site_info[i, "nearest_ws_id"]
	if (FLOWERING_ELSE_VERAISON) {
    temp_by_site <- subset(ws_hourly, weather_station_id == site_ws & 
						                 date(ymd_hms(timestamp)) >= ymd(site_info[i,"date_bb"]) & 
										 date(ymd_hms(timestamp)) <= ymd(site_info[i,"date_fl"]) )
	} else {
    temp_by_site <- subset(ws_hourly, weather_station_id == site_ws & 
						                 date(ymd_hms(timestamp)) >= ymd(site_info[i,"date_fl"]) & 
										 date(ymd_hms(timestamp)) <= ymd(site_info[i,"date_vr"]) )
	}
	if (FALSE) {
		# minute scale
		# must figure out the appropriate distribution.  I'm using univorm, but something else like normal or beta or even von mises may be right
		specgd <- table(t(mapply(function(x,y) round(seq(from=x, to=y, length.out=60)), temp_by_site$temperature_min, temp_by_site$temperature_max)))
	} else {
		# hour scale
		if (TAVG_ELSE_MIDRANGE) {
			specgd <- table(round(temp_by_site$temperature_avg))
		} else {
			specgd <- table(round(( temp_by_site$temperature_max-temp_by_site$temperature_min )/2 + temp_by_site$temperature_min))
		}
	}
	new_names <- paste(names(specgd), 'C', sep='')
	site_info[i, new_names] <- specgd
}
write.csv(site_info, "../data/01step/site_info_specgd_hourly.csv")
