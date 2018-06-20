#loading in the data and converting to time series
#the directory needs to be set before running the program so that the data can be loaded.

install.packages("car")
library(car)
#installing SSA packages
install.packages("Rssa")
library(Rssa)
install.packages("nortest")
library(nortest)
WH.unemp <- read.csv("WAWHAT5URN.csv")
WH.unemp <- ts(WH.unemp$WAWHAT5URN, start = 1990, frequency = 12)

#using boxcox test to see if transformations are necessary
lambda <- boxCox(WH.unemp ~ 1, family = "bcPower")
#because the 95% CI for the loglikelihood includes zero the 
#transformation will be logistic transformation

#transfroming the time series
mi <- min(WH.unemp); ma <- max(WH.unemp)
log.WH <- log((WH.unemp - mi + .1) / (ma - WH.unemp + .1))
log.WH <- ts(log.WH, start = 1990, frequency = 12)

#checking for normality 
shapiro.test(log_WH) #at a significance level of alpha=.005 the data is not normal. 

n <- length(WH.unemp)
ssa.log.WH <- ssa(log.WH, kind = "1d-ssa")
plot(ssa.log.WH) #eigenvalues
plot(ssa.log.WH, type = "vectors", idx = 1:25) #eigenvectors
plot(ssa.log.WH, type = "paired", idx = 1:25) #pairs of eigenvectors
plot(wcor(ssa.log.WH)) #w-correlation matrix
#pairs are c(8,9) and c(14,15)

#below is a table of the periods, frequencies and yearly period and the
#closest matching whole year as well.
estimates.ssa.column <- matrix(0, nrow = 15, ncol = 3)
for ( i in 1 : 1 : 15 ) {
	#period
		estimates.ssa.column[i,1] <- parestimate(ssa.log.WH, list(c(i,(i+1))), method = "pairs",
		subspace = "column")$period
	#frequency
		frequency <- as.numeric(parestimate(ssa.log.WH, list(c(i,(i+1))), method = "pairs",
		subspace = "column")[3])
		estimates.ssa.column[i,2] <- frequency * 2 * pi
	#period with respect to a year
		estimates.ssa.column[i,3] <- parestimate(ssa.log.WH, list(c(i,(i+1))), method = "pairs",
		subsapce = "column")$period / 12 	
	
}
colnames(estimates.ssa.column) <- c("period", "Arg", "yearly")
rownames(estimates.ssa.column) <- c("1,2","2,3","3,4","4,5","5,6","6,7","7,8","8,9","9,10","10,11",
	"11,12","12,13","13,14","14,15","15,16")
pairs.cycle <- matrix(c("5  6","8  9","14  15"), nrow = 3, ncol = 1)
colnames(pairs.cycle) <- c("Pairs")
rownames(pairs.cycle) <- c("12 Months", "6 Months", "18 Months")

#creating the SSA grouping
group.log.WH <- list(c(1),c(2),c(3),c(4),c(7),c(12),c(13),
	c(5,6),c(8,9),c(10,11),c(14,15))
num.comp.log <- length(group.log.WH)
recon.log.WH <- reconstruct(ssa.log.WH, groups = group.log.WH)

#fully reconstructed time series
dlog.WH <- rep(0,n)
for ( i in 1 : num.comp.log )
{
	dlog.WH <- dlog.WH + recon.log.WH[[i]]
}
plot(cbind(log.WH, dlog.WH), plot.type = "single", col = c("black", "red"))

#########################################
#transforming data back 
#########################################
trans.WH <- (ma * exp(log.WH) + .2 * exp(log.WH) + mi - .2) / (1 + exp(log.WH))
trans.ssa.WH <- (ma * exp(dlog.WH) + .2 * exp(dlog.WH) + mi - .2) / (1 + exp(dlog.WH))
plot(cbind(trans.WH, trans.ssa.WH), plot.type = "single", col = c("black", "red"))

#residuals for reconstruction
res <- residuals(recon.log.WH)

pvalues <- double(10)
for ( lag in 1 : 10 )
{
	pvalues[lag] <- Box.test(res, lag = lag)$p.value
}
plot(pvalues)
#conclusion is that there is autocorrelation present in the residuals

#residual reconstruction
library(forecast)
ar.aicc <- auto.arima(res, d = 0, D = 0, max.p = 5, max.q = 5, max.P = 0, max.Q = 0,
	stationary = TRUE, seasonal = FALSE, ic = "aicc", allowmean = FALSE)
ar1 <- arima.sim(n = (length(res)-1), list(ar = 0.2599), sd = sqrt(.08262))
ad.test(ar1)
#conclusion is that ar1 is normally distributed

#forecasting one year ahead
n.ahead <- 44
for.ssa.log.WH <- rforecast(ssa.log.WH, groups = group.log.WH,
	len = n.ahead, only.new = FALSE)
for.log.WH <- rep(0, (n+n.ahead))
for ( i in 1 : num.comp.log )
{
	for.log.WH <- for.log.WH + for.ssa.log.WH[[i]]
}

#forecasted points using the AR1
forecasts <- predict(ar.aicc, n.ahead = n.ahead, interval = "prediction")
log.point.forecast <- c(forecasts$pred) + for.log.WH[(n+1):(n+n.ahead)]

#degrees of freedom for t-dist using kurtosis
install.packages("fBasics")
library(fBasics)
kurt <- kurtosis(res)
df <- (6 + kurt * 4) / kurt
alpha <- .05
quantile.t <- qt(1-alpha/2, df = df)
log.up <- c(forecasts$pred + quantile.t * forecasts$se) + for.log.WH[(n+1):(n+n.ahead)]
log.lo <- c(forecasts$pred - quantile.t * forecasts$se) + for.log.WH[(n+1):(n+n.ahead)]

upper.lim <- ts((ma * exp(log.up) + .2 * exp(log.up) + mi - .2) /
	(1 + exp(log.up)), start = 2018+1/3, frequency = 4)
lower.lim <- ts((ma * exp(log.lo) + .2 * exp(log.lo) + mi - .2) /
	(1 + exp(log.lo)), start = 2018+1/3, frequency = 4)

#transforming prediction model
for.trans.WH <- ts((ma * exp(for.log.WH) + .2 * exp(for.log.WH) + mi - .2) /
	(1 + exp(for.log.WH)), start = 1990, frequency = 12)

par(cex.lab = 2, cex.axis = 2, cex.main = 4)
plot(for.trans.WH, lwd = 2, xlab = "Year", ylab = "Percent Unemployed",
	main = "Forecasted U-6 Unemployment in Whatcom County to 2022")
abline(v = c(1990,1995,2000,2005,2010,2015,2020), h = c(5,6,7,8,9,10), lty = "dashed")
t.pred <- seq(2018+4/12,2022-1/12,1/12)
lines(upper.lim ~ t.pred)
lines(lower.lim ~ t.pred)
mycol.grey <- rgb(190,190,190, max = 255, alpha = 150, "orange")
polygon(c(t.pred,rev(t.pred)), c(lower.lim,rev(upper.lim)),col = mycol.grey, border = NA)
lines(for.trans.WH[(n+1):(n+n.ahead)] ~ t.pred, col = "red", lwd = 2)

#RMSE: quality of prediction
num.window <- 10
log.error <- matrix(0, ncol = n.ahead, nrow = num.window)
for ( w in 1 : num.window )
{
	log.WH.win <- for.log.WH[w:(w+n-1)]
	log.WH.win.ssa <- ssa(log.WH.win, kind = "1d-ssa")
	recon.log.win <- reconstruct(log.WH.win.ssa, groups = group.log.WH)
	
	d.log.win <- rep(0,n)
	for ( i in 1 : num.comp.log )
	{
		d.log.win <- d.log.win + recon.log.win[[i]]
	}
	
	for.log.win.ssa <- rforecast(log.WH.win.ssa, groups = group.log.WH,
		len = n.ahead, only.new = FALSE)
	for.log.win <- rep(0,(n+n.ahead))
	for ( i in 1 : num.comp.log) 
	{
		for.log.win <- for.log.win + for.log.win.ssa[[i]]
	}

	res.det.win <- log.WH.win - d.log.win
	ar.aicc.win <- Arima(res.det.win, order = c(1,0,0), seasonal = c(0,0,0),
		include.mean = FALSE, method = "CSS", lambda = NULL)

	ar.res.win <- ar.aicc.win$residuals
	for.win <- predict(ar.aicc.win, n.ahead = n.ahead)
	log.point.for.win <- c(for.win$pred) + for.log.win[(n+1):(n+n.ahead)]
	log.actual.win <- for.log.WH[(w+n):(w+n-1+n.ahead)]
	log.error[w,] <- log.point.for.win - log.actual.win
	print(w)
}

rmse <- sqrt(colMeans(log.error, na.rm = TRUE)^2)
plot(rmse, main = "RMSE Values for each Predicted Month", xlab = "Predicted Month",
	ylab = "RMSE", col = "black", pch = 16, cex = 2)
abline(h = c(.2,.4,.6,.8,1), v = c(10,20,30,40), lty = "dashed")
##############
