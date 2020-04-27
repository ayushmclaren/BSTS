library(lubridate)
library(bsts)
library(dplyr)
library(ggplot2)
library(forecast)
library(reshape2)

### Load the data
data("AirPassengers")
Y <- window(AirPassengers, start=c(1949, 1), end=c(1959,12))
y <- log10(Y)
### Fit the ARIMA model
arima <- arima(y, 
               order=c(0, 1, 1), 
               seasonal=list(order=c(0,1,1), period=12))

### Actual versus predicted
d1 <- data.frame(c(10^as.numeric(fitted(arima)), # fitted and predicted
                   10^as.numeric(predict(arima, n.ahead = 12)$pred)),
                 as.numeric(AirPassengers), #actual values
                 as.Date(time(AirPassengers)))
names(d1) <- c("Fitted", "Actual", "Date")


### MAPE (mean absolute percentage error)
MAPE <- filter(d1, year(Date)>1959) %>% summarise(MAPE=mean(abs(Actual-Fitted)/Actual))

### Plot actual versus predicted
ggplot(data=d1, aes(x=Date)) +
  geom_line(aes(y=Actual, colour = "Actual"), size=1.2) +
  geom_line(aes(y=Fitted, colour = "Fitted"), size=1.2, linetype=2) +
  theme_bw() + theme(legend.title = element_blank()) + 
  ylab("") + xlab("") +
  geom_vline(xintercept=as.numeric(as.Date("1959-12-01")), linetype=2) +
  ggtitle(paste0("ARIMA -- Holdout MAPE = ", round(100*MAPE,2), "%")) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))


### Run the bsts model
ss <- AddLocalLinearTrend(list(), y)
ss <- AddSeasonal(ss, y, nseasons = 12)
bsts.model <- bsts(y, state.specification = ss, niter = 1000, ping=0, seed=2016)

### Get a suggested number of burn-ins
burn <- SuggestBurn(0.1, bsts.model)

### Predict
p <- predict.bsts(bsts.model, horizon = 12, burn = burn, quantiles = c(.025, .975))

### Actual versus predicted
d2 <- data.frame(
  # fitted values and predictions
  c(10^as.numeric(-colMeans(bsts.model$one.step.prediction.errors[-(1:burn),])+y),  
    10^as.numeric(p$mean)),
  # actual data and dates 
  as.numeric(AirPassengers),
  as.Date(time(AirPassengers)))
names(d2) <- c("Fitted", "Actual", "Date")

### MAPE (mean absolute percentage error)
MAPE <- filter(d2, year(Date)>1959) %>% summarise(MAPE=mean(abs(Actual-Fitted)/Actual))

### 95% forecast credible interval
posterior.interval <- cbind.data.frame(
  10^as.numeric(p$interval[1,]),
  10^as.numeric(p$interval[2,]), 
  subset(d2, year(Date)>1959)$Date)
names(posterior.interval) <- c("LL", "UL", "Date")

### Join intervals to the forecast
d3 <- left_join(d2, posterior.interval, by="Date")

### Plot actual versus predicted with credible intervals for the holdout period
ggplot(data=d3, aes(x=Date)) +
  geom_line(aes(y=Actual, colour = "Actual"), size=1.2) +
  geom_line(aes(y=Fitted, colour = "Fitted"), size=1.2, linetype=2) +
  theme_bw() + theme(legend.title = element_blank()) + ylab("") + xlab("") +
  geom_vline(xintercept=as.numeric(as.Date("1959-12-01")), linetype=2) + 
  geom_ribbon(aes(ymin=LL, ymax=UL), fill="grey", alpha=0.5) +
  ggtitle(paste0("BSTS -- Holdout MAPE = ", round(100*MAPE,2), "%")) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

credible.interval <- cbind.data.frame(
  10^as.numeric(apply(p$distribution, 2,function(f){quantile(f,0.75)})),
  10^as.numeric(apply(p$distribution, 2,function(f){median(f)})),
  10^as.numeric(apply(p$distribution, 2,function(f){quantile(f,0.25)})),
  subset(d3, year(Date)>1959)$Date)
names(credible.interval) <- c("p75", "Median", "p25", "Date")

components <- cbind.data.frame(
  colMeans(bsts.model$state.contributions[-(1:burn),"trend",]),                               
  colMeans(bsts.model$state.contributions[-(1:burn),"seasonal.12.1",]),
  as.Date(time(Y)))  
names(components) <- c("Trend", "Seasonality", "Date")
components <- melt(components, id="Date")
names(components) <- c("Date", "Component", "Value")

### Plot
ggplot(data=components, aes(x=Date, y=Value,group=1)) + geom_line() + 
  theme_bw() + theme(legend.title = element_blank()) + ylab("") + xlab("") + 
  facet_grid(Component ~ ., scales="free") + guides(colour=FALSE) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

###################################################
###################################################

### Fit the model with regressors
data(iclaims)
ss2 <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss2 <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
bsts.reg <- bsts(iclaimsNSA ~ ., state.specification = ss2, data =
                   initial.claims, niter = 500, ping=0, seed=2016)

### Get the number of burn-ins to discard
burn <- SuggestBurn(0.1, bsts.reg)

### Helper function to get the positive mean of a vector
PositiveMean <- function(b) {
  b <- b[abs(b) > 0]
  if (length(b) > 0) 
    return(mean(b))
  return(0)
}

### Get the average coefficients when variables were selected (non-zero slopes)
coeff <- data.frame(melt(apply(bsts.reg$coefficients[-(1:burn),], 2, PositiveMean)))
coeff$Variable <- as.character(row.names(coeff))
ggplot(data=coeff, aes(x=Variable, y=value)) + 
  geom_bar(stat="identity", position="identity") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("") + ylab("") + ggtitle("Average coefficients")

inclusionprobs <- melt(colMeans(bsts.reg$coefficients[-(1:burn),] != 0))
inclusionprobs$Variable <- as.character(row.names(inclusionprobs))
ggplot(data=inclusionprobs, aes(x=Variable, y=value)) + 
  geom_bar(stat="identity", position="identity") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
  xlab("") + ylab("") + ggtitle("Inclusion probabilities")

P75 <- function(b) {
  b <- b[abs(b) > 0]
  if (length(b) > 0) 
    return(quantile(b, 0.75))
  return(0)
}

p75 <- data.frame(melt(apply(bsts.reg$coefficients[-(1:burn),], 2, P75)))


### Get the components
components.withreg <- cbind.data.frame(
  colMeans(bsts.reg$state.contributions[-(1:burn),"trend",]),
  colMeans(bsts.reg$state.contributions[-(1:burn),"seasonal.52.1",]),
  colMeans(bsts.reg$state.contributions[-(1:burn),"regression",]),
  as.Date(time(initial.claims)))  
names(components.withreg) <- c("Trend", "Seasonality", "Regression", "Date")
components.withreg <- melt(components.withreg, id.vars="Date")
names(components.withreg) <- c("Date", "Component", "Value")

ggplot(data=components.withreg, aes(x=Date, y=Value)) + geom_line() + 
  theme_bw() + theme(legend.title = element_blank()) + ylab("") + xlab("") + 
  facet_grid(Component ~ ., scales="free") + guides(colour=FALSE) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

##########
#Prior Expectations

prior <- SpikeSlabPrior(x=model.matrix(iclaimsNSA ~ ., data=initial.claims), 
                        y=initial.claims$iclaimsNSA, 
                        prior.information.weight = 200,
                        prior.inclusion.probabilities = prior.spikes,
                        optional.coefficient.estimate = prior.mean)

### Run the bsts model with the specified priors
data(iclaims)
ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
bsts.reg.priors <- bsts(iclaimsNSA ~ ., state.specification = ss, 
                        data = initial.claims, 
                        niter = 500, 
                        prior=prior, 
                        ping=0, seed=2016)


### Get the average coefficients when variables were selected (non-zero slopes)
coeff <- data.frame(melt(apply(bsts.reg.priors$coefficients[-(1:burn),], 2, PositiveMean)))
coeff$Variable <- as.character(row.names(coeff))
ggplot(data=coeff, aes(x=Variable, y=value)) + 
  geom_bar(stat="identity", position="identity") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
  xlab("") + ylab("")
