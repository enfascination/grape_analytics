print(" use a linear model to predict veraison using GDD, very simple")
print("the takeaways of this are the lm out of sample prediction (data-driven)is way better than use of a critical threshold (theory-driven), and also that the simple, data-first GDD is ultimately better in this context than the action units acculumated by more theory-driven models like wang engel")
library(reshape2)
library(stringr)
"%ni%" <- Negate("%in%") # https://stackoverflow.com/questions/5831794/opposite-of-in-exclude-rows-with-values-specified-in-a-vector

source("01_stat_functions.R")

env <- Sys.getenv()
FLOWERING_ELSE_VERAISON <- as.logical(env["FLOWERING_ELSE_VERAISON"])
if (is.na(FLOWERING_ELSE_VERAISON)) { stop("is env loading?") }
if (FLOWERING_ELSE_VERAISON) {
    print("set for FLOWERING")
} else {
    print("set for VERAISON")
}

site_info <- read.csv("../data/01step/site_info_gdd.csv")
site_info <- site_info[,-1]

site_info_fit <- site_info
if (FLOWERING_ELSE_VERAISON) {
    site_info_fit$GDD_prev <- site_info_fit$GDD_bb
    site_info_fit$GDD <- site_info_fit$GDD_fl
    site_info_fit$doy <- site_info_fit$doy_fl
} else {
    site_info_fit$GDD_prev <- site_info_fit$GDD_fl
    site_info_fit$GDD <- site_info_fit$GDD_vr
    site_info_fit$doy <- site_info_fit$doy_vr
}

# remove NAs (sites that have GDD but aren't in the dataset
#   could also remove NAs on GDD_fl, but that's not what I'm predciting at the moment, I just ocmputed it but not using
na.idxs <- !is.na(site_info_fit$GDD)
site_info_fit <- site_info_fit[na.idxs,]

# BUILD MODEL
f_gdd_fl <- as.formula("doy ~ GDD_prev")
f_gdd <- as.formula("doy ~ GDD")
f_gdd_alt <- as.formula("doy ~ gdd_vr_alt")
f_gdd_phen <- as.formula("doy ~ GDD_fl + GDD_vr")
f_actionunits <- as.formula("doy ~ units_bb + units_fl + units_vr")
f_action_vr_fl <- as.formula("doy ~ units_fl + units_vr")
f_action_vr <- as.formula("doy ~ units_vr")
f_hybrid <- as.formula("doy ~ GDD_vr + units_vr")
f_phenmodel <- f_gdd
phen.model <- lm(f_phenmodel, data=site_info_fit)
#summary(phen.model)


# PREDICT WITH MODEL
print("IN SAMPLE")
predict(phen.model)
print( length(predict(phen.model)) )
print(c("r squared", summary(phen.model)$r.squared))
print(c("RMSE", RMSE(site_info_fit$doy, predict(phen.model))))

# PREDICT WITH MODEL: OUT OF SAMPLE CIs WITH OUT OF SAMPLE BOOTSTRAP
l <- nrow(site_info_fit)
n <- 1000
rmses <- rep(0, n)
for (i in 1:n) {
    all_idxs <- 1:l
    ### sample with replacement
    boot_idxs <- sample(all_idxs, replace=TRUE)
    ### use efron's .632 rule to select remaining unsampled observations for out of sample prediction on each bootstrap run
    out_idxs <- all_idxs[all_idxs %ni% boot_idxs]
    # fit model
    phen.model <- lm(f_phenmodel, data=site_info_fit[boot_idxs,])
    # predict out of sample
    doy_preds <- predict(phen.model, newdata=site_info_fit[out_idxs,])
    # rmse
    rmses[i] <- RMSE(site_info_fit[out_idxs,"doy"], doy_preds)
}
print("OUT OF SAMPLE")
print(mean(rmses))
print(quantile(rmses))
