#################################################
########## Data Import and manipulation #########
#################################################
## Vacicek model on real numbers; Estimation of the risk premium and the fitted yield-curves. 
## Read in (weekly) rate observation
rm(list = ls())
IRdata <- read.csv('F-F_Research_Data_Factors_weekly.csv', skip = 4)
IRdata <- IRdata[ c("X", "RF") ]

## Look at data...
head(IRdata)

## Extract relevant columns
dates <- as.Date(IRdata[,1], "%Y%m%d")
rf_weekly <- as.numeric(IRdata$RF)

## Remove NA's
idxRemove <- is.na(dates) | is.na(rf_weekly)
dates <- dates[!idxRemove]
timepts = (as.numeric(dates) - as.numeric(dates)[1]) / 365
rf_weekly <- rf_weekly[!idxRemove]

## Convert weekly rates of return to yearly continuously compounded
rf <- log((1 + rf_weekly/100)^(52)) # Divide by 100 to get decimals...

## Implement changes in original dataset as we will use it later...
IRdata <- IRdata[!idxRemove, ]
IRdata$RF <- rf

## Plot it...
plot(dates, rf, type = "l", main = "Historical risk-free rate", 
     xlab  = "date", ylab = "risk free rate")

## Assume time steps are equidistant 1-week
dt <- 7 / 365

############## Estimate Analytically ############

## Estimate parameters
r0 <- rf[1]
n <- length(rf) - 1

ri <- rf[2:n]
rim1 <- rf[1:(n-1)]


a <- (sum(ri)*sum(rim1) - n*sum(ri*rim1)) / (sum(rim1)^2 - n*sum(rim1^2))
b <- (sum(ri) - a*sum(rim1))  / n
v <- (sum((ri - rim1*a - b)^2))  / n

## Convert back to original parameters
kappaHat <- -log(a) / dt
thetaHat <- b / (1-a)
sigHat <- sqrt((2*kappaHat*v) / ( 1 - a^2))

PparamHat <- c(thetaHat, kappaHat, sigHat)
PparamHat

##     Plot theta and conditional mean         ##

CondMeanVasicek <- function(xt, u, theta, kappa){
  xt*exp(-kappa*u) + theta*(1 - exp(-kappa*u))
}

u <- seq(0, 20, 0.02)
projDates <- dates[length(dates)] + u*365
rf_projection <- CondMeanVasicek(rf[length(rf)], u, thetaHat, kappaHat)

library(ggplot2)
df <- data.frame(dates = dates, rate = rf, series = "US 3m rate")
dfProc <- data.frame(dates = projDates, rate = rf_projection, series = "Conditional Expectation")
dfTheta <- data.frame(dates = c(dates, projDates), 
                      rate = rep(thetaHat, length(dates) + length(projDates)),
                      series = "thetaHat (long term level)")
dfPlot <- rbind(df, dfProc, dfTheta)
ggplot(dfPlot, aes(x = dates, y = rate, colour = series)) + geom_line() + ggtitle("US 3m rate")



#################################################
############## Estimate Numerically #############
#################################################

# Here we can be exact wrt. time steps

VASICEKnegloglike<-function(param,data,times){
  n.obs<-length(data)
  dum<-sort.list(times)
  data<-data[dum]
  times<-times[dum]
  delta<-diff(times,1)
  mv<-data[1:(n.obs-1)]*exp(-param[2]*delta)+param[1]*(1-exp(-param[2]*delta))
  variance<-param[3]^2*(1-exp(-2*param[2]*delta))/(2*param[2])
  VASICEKnegloglike<--sum(log(dnorm(data[2:n.obs],mv,sqrt(variance))))
}

tstVAS <- optim(par=c(thetaHat, kappaHat, sigHat),fn=VASICEKnegloglike,method = "BFGS",  data=rf,times=timepts, hessian = TRUE)
tstVAS$par
c(thetaHat, kappaHat, sigHat)


# No real change in estimated parameters... 

#################################################################
############## Standard Errors and Confidence Bands #############
#################################################################

## 'Optim' returns the hessian of the objective function evaluated
## at the optimal point. Since we optimize the negative 
## log-likelihood this is exactly -1 times the hessian of the 
## log-likelihood evaluated at the ML estimate.

covMatrix <- solve(tstVAS$hessian)
stdErrs <- sqrt(diag(covMatrix))
stdErrs

## 95% confidence intervals:
CI <- PparamHat + 1.96*stdErrs%*%t(c(-1, 0, 1))

rownames(CI) <- c("thetaHat", "kappaHat", "sigHat")
colnames(CI) <- c("2.5%", "MLE", "97.5%")

CI


############## Investigate Yield Curves #########

# IMPORT: US Treasury Zero-Coupon Yield Curve (from Quandl)

yieldData <- read.csv('FED-SVENY.csv')
head(yieldData)
mat <- 1:30
yieldData[, mat + 1] <- yieldData[, mat + 1] / 100

## Merge with risk-free data set
IRdata$Date <- as.Date(IRdata$X, "%Y%m%d")
yieldData$Date <- as.Date(yieldData$Date, "%Y-%m-%d")
allData <- merge(IRdata,yieldData,by="Date")
head(allData)

## Example yield curve:
plotIdx <- 2000
plot(as.numeric(allData[plotIdx, mat + 3]), xlab = "maturity (in years)", 
     ylab = "yield", main = paste("US Treasury Zero-Coupon Yield Curve (", allData[plotIdx, 1], ")"
     ), type = "l")


###################### Average Historical Yield Curve ##############


# Full 30 year yield curves are not available prior to ~ 1985 thus we limit the dataset:
allDataSub <- allData[allData$Date >= as.Date("1985-11-25","%Y-%m-%d"),]
mat <- 1:30
avgYield <- colMeans(allDataSub[, 3 + mat],na.rm = FALSE)

dfTmp <- data.frame(maturity = mat, yield = avgYield, series = "Historical average")
dfPlotData <- dfTmp 

# Plot
ggplot(dfPlotData, aes(x = maturity, y = yield, color = series)) + geom_line() + ggtitle("Historical average (i.e. 'typical') yield curve")


# if P = Q, ie. riskpremium = tilde{lambda} = 0, the 'typical' vasicek curve
# should match the average observed yield curve well

####################################################################
################## Typical yield curve under Vasicek (rp = 0) ######
####################################################################
VASICEKyield<-function(r,tau,Pparam,riskpremium=0)
{ b<-(Pparam[1]+riskpremium)*Pparam[2]
a<-Pparam[2]
sig<-Pparam[3]
Btau<-(1-exp(-a*tau))/a
Atau<-((Btau-tau)*(a*b-0.5*sig^2)/a^2 - sig^2*Btau^2/(4*a))
VASICEKyield<-r*Btau/tau-Atau/tau
}

# Calculate average historical short rate
avgRF <- mean(allDataSub$RF)

vasicek.yield.avgRF <- VASICEKyield(avgRF, mat, PparamHat, riskpremium=0) #USING rp = 0

dfTmp <- data.frame(maturity = mat, yield = vasicek.yield.avgRF, series = "Typical Vasicek curve (r0 = hist. avg, risk premium = 0)")
dfPlotData <- rbind(dfPlotData, dfTmp)

# Plot
ggplot(dfPlotData, aes(x = maturity, y = yield, color = series)) + geom_line(size = 1.5) + ggtitle("Yield curves") + ylim(low = 0, high = 0.07)

# Conclusion: That does not work so well... Most likely risk premium > 0

####################################################################
################## Finding a reasonable risk premium ###############
####################################################################

## Let us assume that the risk premium tilde{lambda} is a constant  
## and let us estimate it as what gives the best fit of the 'typical' 
## vasicek curve to the 'typical' observed yield curve.

## Minimize sum of absolute differences between hist. avg. curve
## and 'typical' Vasicek curve using different (constant) risk premia (rp):

rpfit <- function(rp){
  rpfit <- sum(abs(avgYield - VASICEKyield(avgRF,mat,PparamHat,riskpremium=rp)))
}

rpFit <- optim(par=c(0),fn=rpfit,method = "BFGS")$par
rpFit

####################################################################
###################### Does it look better now? ####################
####################################################################

vasicek.yield.with.rp <- VASICEKyield(avgRF, mat, PparamHat, riskpremium=rpFit)
dfTmp <- data.frame(maturity = mat, yield = vasicek.yield.with.rp, series = paste(
  "Typical Vasicek curve (risk-free = hist. avg, risk premium = ", round(rpFit, 4), ")"))
dfPlotData <- rbind(dfPlotData, dfTmp)
ggplot(dfPlotData, aes(x = maturity, y = yield, color = series)) + geom_line(size = 1.5) + ggtitle("Yield curves") + ylim(low = 0, high = 0.07)

## Yes, much better (on 'average' at least...)

####################################################################
###################### What about on specific days? ################
####################################################################

# Theory says that the model (with correct parameters) should match 
# the yield curve not just "on average" but EVERY day

# Sample date:
idxPlot <- 1715

df1 <- data.frame(y = VASICEKyield(allDataSub[idxPlot, ]$RF,mat,PparamHat,riskpremium=rpFit), 
                  tau = mat, series = rep(paste("Vasicek with risk premium = ", round(rpFit, 4)), length(mat)))
df2 <- data.frame(y = as.numeric(t(allDataSub[idxPlot, 3 + mat])), tau = mat, series = rep("observed", length(mat)))
dfPlot <- rbind(df1, df2)

# Plot
ggplot(dfPlot, aes(x = tau, y = y, colour = series)) + 
  geom_line() + 
  ggtitle(paste("US Treasury Zero-Coupon Yield Curve (", allDataSub[idxPlot, 1], ")"))




