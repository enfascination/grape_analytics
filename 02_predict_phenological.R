source("01_stat_functions.R")
source("01_functions_pheno_models.R")
library(lubridate)
print("use fitted phenological model to predict veraison or flowering in forrestel data")

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
Year.target <- 2015
Year.previous <- 2014
hemisphere <- "NORTH"

#' @param variety.parameters - fixed parameters for each variety in:
#'        ../data/phenology_parameterization/varieties_parameters_best.csv
variety.parameters <- read.csv("../data/phenology_parameterization/varieties_parameters_best.csv")

site_info <- read.csv("../data/01step/site_info.csv")
site_info <- site_info[,-1]
climate.data <- read.csv("../data/01step/climate.data.csv")
climate.data <- climate.data[,-1]

Variety <- "Cabernet-Sauvignon" # "CABERNET_SAUVIGNON" in the Forrestel data
parameters <- variety.parameters[which(variety.parameters[,"Variety"]==Variety),2:23]
### add crits for gdd model
parameters$GDD_crit_0_BB <- 52.5
parameters$GDD_crit_BB_FL <- 350
parameters$GDD_crit_FL_VR <- 725
parameters$GDD_start <- 366

if (FLOWERING_ELSE_VERAISON) {
    obs_dates <- site_info$doy_fl
} else {
    obs_dates <- site_info$doy_vr
}

predict_dates_gdd <- function(temps.avgs, parameters) {
    valid <- 1
      ### calc GDD
      #temps.avgs$GDD <- with(temps.avgs, max(0, TAV - 10))
      temps.avgs$GDD <- with(temps.avgs, ifelse( TAV - 10 > 0, TAV - 10, 0))
      #temps.avgs$TMIDRANGE <- with(temps.avgs, (TMAX - TMIN)/2 + TMIN)
      #temps.avgs$GDD <- with(temps.avgs, ifelse( TMIDRANGE - 10 > 0, TMIDRANGE - 10, 0))
      temps.avgs$GDDsum <- cumsum(temps.avgs$GDD)
      ActionunitsGDD <- temps.avgs[,c("GDD", "GDDsum")]
      # predictons by GDD 
      BB.date <- list(date=data.frame(Year=2015, date="2015-01-01"), julian=parameters$GDD_start)
      Flow.date.GDD <- GetPhenDate(parameters$GDD_crit_BB_FL, temps.avgs,ActionunitsGDD,BB.date)
      #GDD_startdate <- list(date=data.frame(Year=2015, date="2015-03-15"), julian=439)
      Verais.date.GDD <- GetPhenDate(parameters$GDD_crit_FL_VR, temps.avgs,ActionunitsGDD,Flow.date.GDD)
    return(list(
                "Flow.date"= Flow.date.GDD,
                "Verais.date"= Verais.date.GDD,
                "ActionunitsVer"= ActionunitsGDD,
                "valid"=valid
                ))
}

predict_dates_phenological <- function(temps.avgs, parameters) {
  # params that are pulled from parameters
    # Tm1,Topt,Topt.1,Tn2,minn, Tmin_BB_FL,Tmax_BB_FL,Topt_BB_FL, Tmin_FL_VER,Tmax_FL_VER,Topt_FL_VER, Cstar_BB, Fstar_BB, Fstar_BB_Fl_obs, Fstar_Fl_Ver_obs,
    # INTERMEDIATE VARIABLES
  valid <- 1

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
    if( is.na(Cstar.date[2])) { 
      valid <- 0 
    } else {
      if( Cstar.date[[2]] < 1) valid <- 0
      if( Cstar.date[[2]] > nrow(temps.avgs)) valid <- 0
    }
    if (valid==1) {
      BB.date <- GetPhenDate(parameters$Fstar_BB,temps.avgs,ActionunitsBB,Cstar.date)
      if( is.na(BB.date$julian)) { 
        valid <- 0 
      } else {
        if( BB.date$julian < 1) valid <- 0
        if( BB.date$julian > nrow(temps.avgs)) valid <- 0
      }
    }
    if (valid==1) {
      Flow.date <- GetPhenDate(parameters$Fstar_BB_Fl_obs,temps.avgs,ActionunitsFl,BB.date)
      if( is.na(Flow.date$julian)) { 
        valid <- 0 
      } else {
        if( Flow.date$julian < 1) valid <- 0
        if( Flow.date$julian > nrow(temps.avgs)) valid <- 0
      }
    }
    if (valid==1) {
      Verais.date <- GetPhenDate(parameters$Fstar_Fl_Ver_obs,temps.avgs,ActionunitsVer,Flow.date)
      #print(c(Verais.date, parameters$Fstar_Fl_Ver_obs,Flow.date))
      if( is.na(Verais.date$julian)) { 
        valid <- 0 
      } else {
        if( Verais.date$julian < 1) { valid <- 0 }
        if( Verais.date$julian > nrow(temps.avgs)) { valid <- 0 }
      }
    }
    # predictions by stage
    #    prediction starting at observed budbreak to predicted flower
    if (FALSE) {
      print("add i as input")
      phase_dates <- list(date=data.frame(Year=year(ymd(site_info[i,"date_bb"])), date=site_info[i,"date_bb"]), julian=site_info[i,"doy_bb"])
      Flow.date.BBobs <- GetPhenDate(parameters$Fstar_BB_Fl_obs,temps.avgs,ActionunitsFl,phase_dates)
      #    prediction starting at flower to predicted veraison
      phase_dates <- list(date=data.frame(Year=year(ymd(site_info[i,"date_fl"])), date=site_info[i,"date_fl"]), julian=site_info[i,"doy_fl"])
      Verais.date.FLobs <- GetPhenDate(parameters$Fstar_Fl_Ver_obs,temps.avgs,ActionunitsVer,phase_dates)
    }
    if (valid == 0) { return(list( "valid"=valid)) }
    out <- list(
                "Cstar.date"= Cstar.date,
                "BB.date"= BB.date,
                "Flow.date"= Flow.date,
                "Verais.date"= Verais.date,
                "Chillunits "= Chillunits ,
                "ActionunitsBB"= ActionunitsBB,
                "ActionunitsFl"= ActionunitsFl,
                "ActionunitsVer"= ActionunitsVer,
                "valid"=valid
                )
    #print(out)
    return(out)
}

fit_all_sites <- function(vals, obs_dates, site_info, climate.data, Year.target, Year.previous, model_type="dryrun") {

  ### parameter processing/checking
  if (model_type == "phenology") {
    valid <- 1
    p <- parameters
    # Tm1,Topt,Tn2,minn,Tmin,Topt.1,Tmax, Tmin_BB_FL,Topt_BB_FL,Tmax_BB_FL, Tmin_FL_VER,Topt_FL_VER,Tmax_FL_VER, Cstar_BB, Fstar_BB, Fstar_BB_Fl_obs, Fstar_Fl_Ver_obs,
    # Cstar_BB, Fstar_BB
    # Cstar_BB, Fstar_BB, Fstar_BB_Fl_obs, Fstar_Fl_Ver_obs,
    if (FALSE) {
    p$Topt.1 <- vals[1]
    p$Cstar_BB <- vals[2]
    p$Fstar_BB <- vals[3]
    p$Fstar_BB_Fl_obs <- vals[4]
    p$Fstar_BB_Fl_obs_min <- vals[4]
    p$Fstar_BB_Fl_obs_max <- vals[4]
    p$Fstar_Fl_Ver_obs <- vals[5]
    p$Fstar_Fl_Ver_obs_min <- vals[5]
    p$Fstar_Fl_Ver_obs_max <- vals[5]
    p$minn <- vals[6]
    }
    if (TRUE) {
      p$Tm1 <- vals[1]
      p$Topt <- vals[2]
      p$Tn2 <- vals[3]
      p$minn <- vals[4]
      p$Tmin <- 0
      p$Topt.1 <- vals[5]
      p$Tmax <- 40
      p$Tmin_BB_FL <- 0
      p$Topt_BB_FL <- vals[5]
      p$Tmax_BB_FL <- 40
      p$Tmin_FL_VER <- 0
      p$Topt_FL_VER <- vals[5]
      p$Tmax_FL_VER <- 40
      p$Cstar_BB <- vals[6]
      p$Fstar_BB <- vals[7]
      p$Fstar_BB_Fl_obs <- vals[8]
      p$Fstar_BB_Fl_obs_min <- vals[8]
      p$Fstar_BB_Fl_obs_max <- vals[8]
      p$Fstar_Fl_Ver_obs <- vals[9]
      p$Fstar_Fl_Ver_obs_min <- vals[9]
      p$Fstar_Fl_Ver_obs_max <- vals[9]
      parameters <- p 
    }
    if (-20 > p$Tm1) valid <- 0
    if (p$Tm1 > p$Topt) valid <- 0
    if (0 > p$Topt) valid <- 0
    if (p$Topt > p$Tn2) valid <- 0
    if (0 > p$Topt.1) valid <- 0
    if (p$Topt.1 > 40) valid <- 0
    if (p$Tmin > p$Topt.1) valid <- 0
    if (p$Topt.1 > p$Tmax) valid <- 0
    if (p$Cstar_BB < 0) valid <- 0
    if (p$Fstar_BB < 0) valid <- 0
    if (p$Fstar_BB_Fl_obs < p$Fstar_BB) valid <- 0
    if (p$Fstar_Fl_Ver_obs < p$Fstar_BB_Fl_obs) valid <- 0
    if (valid == 0) {
      return( 9999999)
    }
    #print(vals)
    #print(parameters)
  } else if (model_type == "GDD") {
    parameters$GDD_crit_BB_FL <- vals[1]
    parameters$GDD_crit_FL_VR <- vals[2]
    parameters$GDD_start <- vals[3]
    #print(c(3, parameters, str(parameters)))
    if (parameters$GDD_crit_BB_FL < 0) {
      return( 9999999)
    } else if ( parameters$GDD_crit_FL_VR < parameters$GDD_crit_BB_FL ) {
      return( 9999999)
    }
    if (parameters$GDD_start < 1) {
      return( 9999999)
    } else if (parameters$GDD_start > 365+365) {
      return( 9999999)
    }
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

    if (model_type == "phenology") {
      rs <- predict_dates_phenological(temps.avgs, parameters)
    } else if (model_type == "GDD") {
      rs <- predict_dates_gdd(temps.avgs, parameters)
    }
    if (rs$valid == 0) return(9999999)

    if (FLOWERING_ELSE_VERAISON) {
      # Single stage prediction
      #thedoy  <- rs$Flow.date.BBobs$julian
      # all stages prediction
      thedoy <- rs$Flow.date$julian
    } else {
      # Single stage prediction
      #thedoy  <- rs$Verais.date.FLobs$julian
      # all stages prediction
      thedoy <- rs$Verais.date$julian 
    }

    #   365 because growing cycle is more than a year
    pred_dates[i] <- thedoy %% 365
    error <- RMSE(obs_dates, pred_dates)
  }
  return(error)
}

pred_dates <- fit_all_sites(parameters, obs_dates, site_info, climate.data, Year.target, Year.previous, model_type="GDD")

if (FALSE) {
  #parameters$GDD_crit_0_BB <- 52.5
  #parameters$GDD_crit_BB_FL <- 350
  #parameters$GDD_crit_FL_VR <- 725
  res_gdd <- optim(c(350,725,366), function(init) fit_all_sites(init, obs_dates, site_info, climate.data, Year.target, Year.previous, model_type="GDD"), method="Nelder-Mead")
  #for GDD: best realistic fits are 530 and 738 (a lot of variety for close fits), for error down to 28.7 (still high)
  #for GDD: best less realistic fit is much better 37, 768, 477 , for error down to 10.17. day 477 is late april
}

if (FALSE) {
init <- c( -19.211955,  0.048543, 32.084374 , 0.534960, 19.841501,178.270852,  25.557653 ,41.114928 , 88.116637)
init <- with(parameters, c(Tm1, Topt, Tn2, minn, Topt.1, Cstar_BB, Fstar_BB,  Fstar_BB_Fl_obs, Fstar_Fl_Ver_obs))
res <- optim(init, function(init) fit_all_sites(init, obs_dates, site_info, climate.data, Year.target, Year.previous, model_type="phenology"), method="Nelder-Mead")
#for process to verasion best fit:
# $par
# [1] -19.211955   0.048543  32.084374   0.534960  19.841501 178.270852  25.557653 41.114928  88.116637
# $value
# [1] 8.0697
}
 
if (FALSE) {
init <- with(parameters, c(Cstar_BB, Fstar_BB))
init <- with(parameters, c(Topt.1, Cstar_BB, Fstar_BB,  Fstar_BB_Fl_obs, Fstar_Fl_Ver_obs, minn))
res <- optim(init, function(x) fit_all_sites(x, obs_dates, site_info, climate.data, Year.target, Year.previous, model_type="phenology"), method="Nelder-Mead")
}
stop()

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
