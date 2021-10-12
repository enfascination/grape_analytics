source("01_stat_functions.R")
source("01_functions_pheno_models.R")
library(lubridate)
print("use phenological model to predict veraison or flowering in forrestel data")

# SANITY CHECKS
aa <- Alphafx(10, 40, 25)
WangEngelfx(10, 40, 25, aa, 30)

# DATA SOURCING
#readRenviron('.') ## rereadthe env file
env <- Sys.getenv()
FLOWERING_ELSE_VERAISON <- as.logical(env["FLOWERING_ELSE_VERAISON"])
#FLOWERING_ELSE_VERAISON <- FALSE #FALSE does veraison, which is where most of the data is
if (FLOWERING_ELSE_VERAISON) {
    print("set for FLOWERING")
} else {
    print("set for VERAISON")
}


# INPUTS
#' @param variety.parameters - fixed parameters for each variety in:
#'        ../data/phenology_parameterization/varieties_parameters_best.csv
variety.parameters <- read.csv("../data/phenology_parameterization/varieties_parameters_best.csv")

site_info <- read.csv("../data/01step/site_info.csv")
site_info <- site_info[,-1]
climate.data <- read.csv("../data/01step/climate.data.csv")
climate.data <- climate.data[,-1]

Variety <- "Cabernet-Sauvignon" # "CABERNET_SAUVIGNON" in the Forrestel data
Year.target <- 2015
Year.previous <- 2014
hemisphere <- "NORTH"

if (FLOWERING_ELSE_VERAISON) {
    obs_dates <- site_info$doy_fl
} else {
    obs_dates <- site_info$doy_vr
}

# PREDICTION
pred_dates <- rep(0, nrow(site_info))
#for (i in 1:2) {
for (i in 1:nrow(site_info)) {
    # DERIVED VARIABLES
    ## into into final
    ## temperature vector for the subsetted year, site and variety
    temps.avg <- climate.data[which(climate.data$Site==site_info[i,"Site"]),]
    temps.avgs <- temps.avg[which(temps.avg$Year==Year.previous | temps.avg$Year==Year.target),]
    parameters <- variety.parameters[which(variety.parameters[,"Variety"]==Variety),2:23]
    ### add crits for gdd model
    parameters$GDD_crit_0_BB <- 52.5
    parameters$GDD_crit_BB_FL <- 350
    parameters$GDD_crit_FL_VER <- 725


    # INTERMEDIATE VARIABLES
    ### calc GDD
    #temps.avgs$GDD <- with(temps.avgs, max(0, TAV - 10))
    temps.avgs$GDD <- with(temps.avgs, ifelse( TAV - 10 > 0, TAV - 10, 0))
    #temps.avgs$TMIDRANGE <- with(temps.avgs, (TMAX - TMIN)/2 + TMIN)
    #temps.avgs$GDD <- with(temps.avgs, ifelse( TMIDRANGE - 10 > 0, TMIDRANGE - 10, 0))
    temps.avgs$GDDsum <- cumsum(temps.avgs$GDD)
    ActionunitsGDD <- temps.avgs[,c("GDD", "GDDsum")]

    ### more sophisticated phenological prediction
    Chillunits <- SmoothUtah(parameters$Tm1,parameters$Topt,parameters$Tn2,parameters$minn,temps.avgs,Year.previous,Year.target,hemisphere)

    alpha.BB <- Alphafx(parameters$Tmin,parameters$Tmax,parameters$Topt.1)
    ActionunitsBB <- WangEngelfx(parameters$Tmin,parameters$Tmax,parameters$Topt.1,alpha.BB,temps.avgs$TAV)

    alpha.Fl <- Alphafx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL)
    ActionunitsFl <- WangEngelfx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL,alpha.Fl,temps.avgs$TAV)

    alpha.Ver <- Alphafx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER)
    ActionunitsVer <- WangEngelfx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER,alpha.Ver,temps.avgs$TAV)

    #FINAL PREDICTED DATES
    # full all-stage sequential predictions
    Cstar.date <- GetCstardate(parameters$Cstar_BB,temps.avgs,Chillunits)
    BB.date <- GetPhenDate(parameters$Fstar_BB,temps.avgs,ActionunitsBB,Cstar.date)
    Flow.date <- GetPhenDate(parameters$Fstar_BB_Fl_obs,temps.avgs,ActionunitsFl,BB.date)
    Verais.date <- GetPhenDate(parameters$Fstar_Fl_Ver_obs,temps.avgs,ActionunitsVer,Flow.date)
    # predictons by GDD 
    #GDD_startdate <- list(date=data.frame(Year=2015, date="2015-01-01"), julian=366)
    Flow.date.GDD <- GetPhenDate(parameters$GDD_crit_BB_FL, temps.avgs,ActionunitsGDD,BB.date)
    #GDD_startdate <- list(date=data.frame(Year=2015, date="2015-03-15"), julian=439)
    Verais.date.GDD <- GetPhenDate(parameters$GDD_crit_FL_VER, temps.avgs,ActionunitsGDD,Flow.date.GDD)
    # predictions by stage
    #    prediction starting at observed budbreak to predicted flower
    phase_dates <- list(date=data.frame(Year=year(ymd(site_info[i,"date_bb"])), date=site_info[i,"date_bb"]), julian=site_info[i,"doy_bb"])
    Flow.date.BBobs <- GetPhenDate(parameters$Fstar_BB_Fl_obs,temps.avgs,ActionunitsFl,phase_dates)
    #    prediction starting at flower to predicted veraison
    phase_dates <- list(date=data.frame(Year=year(ymd(site_info[i,"date_fl"])), date=site_info[i,"date_fl"]), julian=site_info[i,"doy_fl"])
    Verais.date.FLobs <- GetPhenDate(parameters$Fstar_Fl_Ver_obs,temps.avgs,ActionunitsVer,phase_dates)

    if (FLOWERING_ELSE_VERAISON) {
        # Single stage prediction
        #thedoy  <- Flow.date.BBobs$julian
        # all stages prediction
        thedoy <- Flow.date$julian
        #thedoy <- Flow.date.GDD$julian
    } else {
        # Single stage prediction
        #thedoy  <- Verais.date.FLobs$julian
        # all stages prediction
        thedoy <- Verais.date$julian 
        #thedoy <- Verais.date.GDD$julian 
    }

    #   365 because growing cycle is more than a year
    pred_dates[i] <- thedoy %% 365
}

# AMOUNT OF ERROR
print(c("Number of blocks", length(obs_dates)))
print(c(RMSE(obs_dates, pred_dates), MAE(obs_dates, pred_dates)))
### relative RMSE
#print(c(RMSE(obs_dates, pred_dates)/mean(obs_dates), MAE(obs_dates, pred_dates)/mean(obs_dates)))
### model efficiency
print(c(ME(obs_dates, pred_dates)))

# PREDICTIONS, QUALITATIVELY
# DateCstar$Cstardate$date
# DateBB$BBdate$date
# Flow.date$date$date
# Verais.date$date$date
# Flow.date.BBobs$date$date
# Verais.date.FLobs$date$date
#> DateCstar$Cstardate$date
#[1] "2015-01-23"
#> DateBB$BBdate$date
#[1] "2015-03-04"
#> Flow.date$date$date
#[1] "2015-04-30"
#> Verais.date$date$date
#[1] "2015-07-08"
#
#> Flow.date.BBobs$date$date
#[1] "2014-05-16"
#> Verais.date.FLobs$date$date
#[1] "2014-07-16"


# CONCLUSION
# predictions from cstar are within weeks, predictions fromthe previous stage are within days.  acceptable and typical
