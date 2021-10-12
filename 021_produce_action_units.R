source("morales/01phenology_models/01_functions_pheno_models.R")
library(lubridate)
print("get datatables of action units by site, date, and phase for lm prediction to compare with GDD")

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

actionunits <- data.frame(date=character(), identifier=character(), units=numeric(), phase=character())
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
    Chillunits <- SmoothUtah(parameters$Tm1,parameters$Topt,parameters$Tn2,parameters$minn,temps.avgs,Year.previous,Year.target,hemisphere)

    alpha.BB <- Alphafx(parameters$Tmin,parameters$Tmax,parameters$Topt.1)
    ActionunitsBB <- WangEngelfx(parameters$Tmin,parameters$Tmax,parameters$Topt.1,alpha.BB,temps.avgs$TAV)

    alpha.Fl <- Alphafx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL)
    ActionunitsFl <- WangEngelfx(parameters$Tmin_BB_FL,parameters$Tmax_BB_FL,parameters$Topt_BB_FL,alpha.Fl,temps.avgs$TAV)

    alpha.Ver <- Alphafx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER)
    ActionunitsVer <- WangEngelfx(parameters$Tmin_FL_VER,parameters$Tmax_FL_VER,parameters$Topt_FL_VER,alpha.Ver,temps.avgs$TAV)

    ### calc GDD
    #temps.avgs$GDD <- with(temps.avgs, max(0, TAV - 10))
    temps.avgs$GDD <- with(temps.avgs, ifelse( TAV - 10 > 0, TAV - 10, 0))
    #temps.avgs$TMIDRANGE <- with(temps.avgs, (TMAX - TMIN)/2 + TMIN)
    #temps.avgs$GDD <- with(temps.avgs, ifelse( TMIDRANGE - 10 > 0, TMIDRANGE - 10, 0))
    temps.avgs$GDDsum <- cumsum(temps.avgs$GDD)
    ActionunitsGDD <- temps.avgs[,c("GDD", "GDDsum")]

    actionunits <- rbind(actionunits,
                         cbind(date=temps.avgs$date, identifier=site_info[i,"Identifier"], units=ActionunitsBB[,1], phase="BB"),
                         cbind(date=temps.avgs$date, identifier=site_info[i,"Identifier"], units=ActionunitsFl[,1], phase="FL"),
                         cbind(date=temps.avgs$date, identifier=site_info[i,"Identifier"], units=ActionunitsVer[,1], phase="VR"),
                         cbind(date=temps.avgs$date, identifier=site_info[i,"Identifier"], units=ActionunitsGDD[,1], phase="GDD")
                         )

}
write.csv(actionunits, "../data/01step/actionunits.csv")

