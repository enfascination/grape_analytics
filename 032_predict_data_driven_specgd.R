print(" use a linear model to predict veraison using histogram of temps, either in linear model or ridge")
library(reshape2)
library(stringr)
library(glmnet)
"%ni%" <- Negate("%in%") # https://stackoverflow.com/questions/5831794/opposite-of-in-exclude-rows-with-values-specified-in-a-vector

source("01_stat_functions.R")

env <- Sys.getenv()
FLOWERING_ELSE_VERAISON <- as.logical(env["FLOWERING_ELSE_VERAISON"])
LM_ELSE_RIDGE <- FALSE # false = ridge regression
###FLOWERING_ELSE_VERAISON <- FALSE #FALSE does veraison, which is where most of the data is
HOURLY_ELSE_DAILY <- TRUE #hourly is just better? not yet ...
if (LM_ELSE_RIDGE) {
    print("set for LINEAR MODEL")
} else {
    print("set for RIDGE REGRESSION")
}
if (FLOWERING_ELSE_VERAISON) {
    print("set for FLOWERING")
} else {
    print("set for VERAISON")
}
if (HOURLY_ELSE_DAILY) {
    print("set for HOURLY temp data")
} else {
    print("set for DAILY temp data")
}

if (HOURLY_ELSE_DAILY) {
	site_info <- read.csv("../data/winkler2021/01step/site_info_specgd_hourly.csv")
} else {
	site_info <- read.csv("../data/winkler2021/01step/site_info_specgd.csv")
}
site_info <- site_info[,-1]
site_info_fit <- site_info 

if (FLOWERING_ELSE_VERAISON) {
    site_info_fit$doy <- site_info_fit$doy_fl
} else {
    site_info_fit$doy <- site_info_fit$doy_vr
}

# BUILD MODEL
## Create a formula for a model with a large number of variables:
#   from https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
xnam <- names(site_info_fit)[grep("C$", names(site_info_fit))]
# exclude columns that are all zeros (temp never observed)
xnam <- xnam[(colSums(site_info_fit[,xnam]) != 0)]
xnam <- xnam[!is.na(xnam)]
if (LM_ELSE_RIDGE) {
	fmla <- as.formula(paste("doy ~ ", paste(xnam, collapse= "+")))
	phen.model.lm <- lm(fmla, data=site_info_fit)
	pred <- predict(phen.model.lm)
} else {
	#grid = 10^seq(10, -2, length = 100)
	grid <- 1
	#cv.glmnet(y=site_info_fit[,"doy"], x=as.matrix(site_info_fit[,xnam]))
	# shows that good lambdas are within 0.6 and 1.2, so I'm sticking with 1, which I started with.
phen.model.ridge <- glmnet(y=site_info_fit[,"doy"], x=as.matrix(site_info_fit[,xnam]), alpha = 0, lambda = grid)
	pred <- predict(phen.model.ridge, s=grid, newx=as.matrix(site_info_fit[,xnam]), type="response")
}
#summary(phen.model.lm)

# PREDICT WITH MODEL
print("IN SAMPLE")
print(length(pred))
print(RMSE(site_info_fit$doy, pred))

## INSPECTION
#temp range:
temps_dist <- colSums(site_info_fit[,xnam])
# selected coefficients
temps_penalized <- as.numeric(coef(phen.model.ridge ))
names(temps_penalized) <- row.names(coef(phen.model.ridge ))
cbind(temps_dist , weights=temps_penalized[-1] )

if (FALSE) {

# PREDICT WITH MODEL: OUT OF SAMPLE CIs WITH OUT OF SAMPLE BOOTSTRAP
l <- nrow(site_info_fit)
rmses <- rep(0, 100)
for (i in 1:100) {
    all_idxs <- 1:l
    ### sample with replacement
    boot_idxs <- sample(all_idxs, replace=TRUE)
    ### use efron's .632 rule to select remaining unsampled observations for out of sample prediction on each bootstrap run
    out_idxs <- all_idxs[all_idxs %ni% boot_idxs]
    # fit model
	if (TRUE) { # code block for building formula, not optional; helps me think
		xnam <- names(site_info_fit)[grep("C$", names(site_info_fit))]
		# exclude columns that are all zeros (temp never observed)
		xnam <- xnam[(colSums(site_info_fit[,xnam]) != 0)]
    xnam <- xnam[!is.na(xnam)]
		fmla <- as.formula(paste("doy ~ ", paste(xnam, collapse= "+")))
	}
	if (LM_ELSE_RIDGE) {
		phen.model.lm <- lm(fmla, data=site_info_fit[boot_idxs,])
		# predict out of sample
		doy_preds <- predict(phen.model.lm, newdata=site_info_fit[out_idxs,])
	} else {
		lambda <- 1
		phen.model.ridge.boot <- glmnet(y=site_info_fit[boot_idxs,"doy"], x=as.matrix(site_info_fit[boot_idxs,xnam]), alpha = 0, lambda = lambda)
		doy_preds <- predict(phen.model.ridge.boot, s=lambda, newx=as.matrix(site_info_fit[out_idxs,xnam]), type="response")
	}
    # rmse
    rmses[i] <- RMSE(site_info_fit[out_idxs,"doy"], doy_preds)
}
print("OUT OF SAMPLE")
print(mean(rmses))
print(quantile(rmses))

}

