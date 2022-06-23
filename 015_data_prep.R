library(lubridate)
library(assertthat)
#library(tidyr)
#library(stringr)

#get phenological models of 
#https://github.com/MoralesCastilla/PhenoDiversity
#working on Forrestel data for one block in one year.
#finding: predictions are very bad


env <- Sys.getenv()
FLOWERING_ELSE_VERAISON <- as.logical(env["FLOWERING_ELSE_VERAISON"])
if (FLOWERING_ELSE_VERAISON) {
    print("set for FLOWERING")
} else {
    print("set for VERAISON")
}

# OBSERVED DATES
obs_phenology_dates_full <- read.csv("../data/winkler2021/winkler_phenology.csv")
pull_stage <- function(ph_data, stage1, stage2) {
    #usage: pull_stage(phenology.csv subset, "Budbreak" or "Bloom", "Bloom" or "Veraison")
    #usage: pull_stage(obs_phenology_dates_full, "Budbreak", "Bloom")
    #usage: pull_stage(obs_phenology_dates_full, "Bloom", "Veraison")
    block_ids_stage1 <- subset(ph_data, phenological_stage %in% c(stage1))$block_id
    block_ids_stage2 <- subset(ph_data, phenological_stage %in% c(stage2))$block_id
    # get blocks that have both stages
    #     intersect() also performs unique()
    block_ids <- intersect(block_ids_stage1, block_ids_stage2)
    # isolate each stage from all shared blocks
    blocks_stage1 <- subset(ph_data, block_id %in% block_ids & phenological_stage == stage1)
    blocks_stage2 <- subset(ph_data, block_id %in% block_ids & phenological_stage == stage2)
    # remove duplicates, picking first (earliest stage change)
    if (TRUE) {
      blocks_stage1 <- blocks_stage1[match(  unique(blocks_stage1$block_id), blocks_stage1$block_id ),] 
      blocks_stage2 <- blocks_stage2[match(  unique(blocks_stage2$block_id), blocks_stage2$block_id ),] 
    } else {
      # remove duplicates, picking last (latest stage change)
      blocks_stage1 <- blocks_stage1[rev(match(  unique(blocks_stage1$block_id), rev(blocks_stage1$block_id))),] 
      blocks_stage2 <- blocks_stage2[rev(match(  unique(blocks_stage2$block_id), rev(blocks_stage2$block_id))),] 
    }
    phenology_dates <- rbind(blocks_stage1, blocks_stage2)
    phenology_dates$X <- NULL
    phenology_dates1 <- phenology_dates[with(phenology_dates, order(block_id, date)),]
    phenology_dates2 <- phenology_dates[with(phenology_dates, order(block_id, phenological_stage)),]
    assert_that( all( phenology_dates1$date == phenology_dates2$date ) | all( phenology_dates1$date != phenology_dates2$date ), msg="ordering problem, stages are happening in the wrong order")
    return(phenology_dates1)
}

### get observed dates of stages in the right format for prediction
if (FLOWERING_ELSE_VERAISON) {
    ### for another day
    obs_phenology_dates_bb_fl <- pull_stage(obs_phenology_dates_full, "Budbreak", "Bloom")
    sites <- unique(obs_phenology_dates_bb_fl$block_id)
    obs_dates_bb_str <- date(mdy_hm(subset(obs_phenology_dates_bb_fl, phenological_stage == "Budbreak")$date))
    obs_dates_bb <- sapply(obs_dates_bb_str, function(x) yday(ymd(x)), simplify=TRUE)
    obs_dates_fl_str <- date(mdy_hm(subset(obs_phenology_dates_bb_fl, phenological_stage == "Bloom")$date))
    obs_dates_fl <- sapply(obs_dates_fl_str, function(x) yday(ymd(x)), simplify=TRUE)
    # dummmy
    obs_dates_vr_str <- rep('', length(obs_dates_fl))
    obs_dates_vr <- rep(0, length(obs_dates_fl))
} else {
		year_min <- min(year(obs_phenology_dates_full$date))
		year_max <- max(year(obs_phenology_dates_full$date))
    obs_phenology_dates_fl_vr <- pull_stage(subset(obs_phenology_dates_full, year_min == year(date)), "Bloom", "Veraison")
		obs_dates_fl_str <- date(ymd(subset(obs_phenology_dates_fl_vr, phenological_stage == "Bloom")$date))
		obs_dates_vr_str <- date(ymd(subset(obs_phenology_dates_fl_vr, phenological_stage == "Veraison")$date))
		for (year in (year_min+1):year_max) {
			obs_phenology_dates_fl_vr_tmp <- pull_stage(subset(obs_phenology_dates_full, year == year(date)), "Bloom", "Veraison") 
			obs_dates_fl_str <- c(obs_dates_fl_str, date(ymd(subset(obs_phenology_dates_fl_vr_tmp, phenological_stage == "Bloom")$date)))
			obs_dates_vr_str <- c(obs_dates_vr_str, date(ymd(subset(obs_phenology_dates_fl_vr_tmp, phenological_stage == "Veraison")$date)))
			obs_phenology_dates_fl_vr <- rbind(obs_phenology_dates_fl_vr, obs_phenology_dates_fl_vr_tmp )
		}
    sites <- unique(obs_phenology_dates_fl_vr$block_id)
    sites_years <- unique(cbind(block_id=obs_phenology_dates_fl_vr$block_id, year=year(obs_phenology_dates_fl_vr$date)))
    obs_dates_fl <- sapply(obs_dates_fl_str, function(x) yday(ymd(x)), simplify=TRUE)
    obs_dates_vr <- sapply(obs_dates_vr_str, function(x) yday(ymd(x)), simplify=TRUE)
    # dummmy
    obs_dates_bb_str <- rep('', length(obs_dates_fl))
    obs_dates_bb <- rep(0, length(obs_dates_fl))
}


### BUILD SITE INFO
site_info <- data.frame(Site = sites_years[,"block_id"], year = as.numeric(sites_years[,"year"]), 
                        date_bb=obs_dates_bb_str, 
                        date_fl=obs_dates_fl_str, 
                        date_vr=obs_dates_vr_str, 
                        doy_bb=obs_dates_bb, 
                        doy_fl=obs_dates_fl, 
                        doy_vr=obs_dates_vr)

### get metadata
###  [this was useful when there were block names, block codes, and block ids.  now there are just block codes]
# various IDs and source info
#block_info <- read.csv("../data/winkler2015/winkler_block_list.csv")
#tmp <- merge(site_info, block_info, by.x="Site", by.y="BDD.ID")
#assert_that( nrow(tmp) == nrow(site_info), msg="merge problem 1")
#site_info <- tmp

# geographical and variety/ag info
block_metadata <- read.csv("../data/winkler2021/winkler_blocks.csv")
tmp <- merge(site_info, block_metadata[,c("block_id", "name", "gps_location_X", "gps_location_Y", "varietal_name")], by.x="Site", by.y="block_id")
assert_that( nrow(tmp) == nrow(site_info), msg="merge problem 2")
site_info <- tmp
# fails:
#assert_that( all(site_info$name == site_info$Block), "metadata mess")
colnames(site_info)[which(names(site_info) == "name")] <- "Block_alt"

### get weather station info
ws_info <- read.csv("../data/winkler2021/winkler_block_weather_station.csv")
### XXXSHORTTERMHACK TO ELIMINATE A DUPLICATE WEATHER STASTAION
tmp <- merge(site_info, ws_info[,c( "identifier", "Winery", "Block", "nearest_ws_id")], by.x=c( "Site"), by.y=c("identifier"), all.x=TRUE, all.y=FALSE)
tmp <- tmp[with(tmp, order(Site, year)),]
assert_that( nrow(tmp) == nrow(site_info), msg="merge problem: names of blocks prob")
site_info <- tmp
print(c("sites:", nrow(site_info)))
write.csv(site_info, "../data/winkler2021/01step/site_info.csv")

# KEY VARIABLES
#' @param climate.data - data frame containing climate data with 7 columns:
#'        Site | Latitude | Year | DOY | TAV | TMIN | TMAX
#  data mainly from https://docs.google.com/spreadsheets/d/1sSXvNZqmkw13YRHl91axMCcMsn_7YaN6/edit#gid=1763009162
#       and https://docs.google.com/spreadsheets/d/1sSXvNZqmkw13YRHl91axMCcMsn_7YaN6/edit#gid=122919871
#       and https://docs.google.com/spreadsheets/d/19ZKCSbgG-jSa7DuOA-Je3rFc_KKPHFKD/edit#gid=1869453609
climate.data.orig <- read.csv("../data/winkler2021/winkler_weather_daily.csv")
climate.data.orig <- subset(climate.data.orig, !is.na(temp_avg))
climate.data <- subset(climate.data.orig, weather_station_id %in% site_info$nearest_ws_id)
climate.data <- data.frame(ws=climate.data[,c("weather_station_id")], Year=year(ymd(climate.data$date)), DOY=yday(ymd(climate.data$date)), date=climate.data[,c("date")], TAV=climate.data[,c("temp_avg")], TMIN=climate.data[,c("temp_min")], TMAX=climate.data[,c("temp_max")])
if (TRUE) {
	print( "converting from Fahrenheit to celsius!!!")
	### convert from Fahrenheit to celsius
	climate.data$TAV <- (climate.data$TAV - 32)*5/9
	climate.data$TMIN <- (climate.data$TMIN - 32)*5/9
	climate.data$TMAX <- (climate.data$TMAX - 32)*5/9
}
tmp <- merge(climate.data, site_info[,c("nearest_ws_id", "Site", "year")], by.x=c("ws", "Year"), by.y=c("nearest_ws_id", "year"))
#assert_that( nrow(tmp) >= nrow(climate.data), msg="some sites have the same nearest weather station, so this will grow the df to give each site it's own weather data") ### except now that I include years this can break and still be right
# Compare:
# table(climate.data$ws)
# table(tmp$Site)
# (also note missing weather data via incomplete days)
# (also note missing weather data via NA temps)
# (also data on phenology dates is surprisingly spotty, and why all the duplicates?)
# (also two missing values in GDD Master: O-5H and O-5I
# (also how were GDD calc'd even?  I tried to figure from the weather sheet and couldn't
#   (and how you predict data from GDD?  like I did: a regression?  from which phase?)
# (also hourly data has many missing observations at station 746 Virtual 691 Vfk, affecting two blocks
# ( phenols: some missing data.  why only flavanoids?  how to use timing? explain timing (why the dates of the 6 visits?)

climate.data <- tmp  ## no longer by weather station, now by site, which is different since sites share weatherstations
									## sites without weather stations have been excluded during merge. alternative is all.y=TRUE which will create NAs instead
# morales code (tenuously) assumes that climate data are ordered chronologically
climate.data <- climate.data[with(climate.data, order(Site, Year, DOY)),]
climate.data.tmp <- climate.data[with(climate.data, order(Site, date)),]
assert_that(all(climate.data$TAV == climate.data.tmp$TAV ))

write.csv(climate.data, "../data/winkler2021/01step/climate.data.csv")

