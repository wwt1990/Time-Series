##  Our function for the TRUE SPECTRUM (on the decibel scale) 
##  of an ARMA process.
##  Each time you start R, and want to use this function,
##  you will have to copy it into R.  Or...  Have it load automatically...
my.spectrum <- function(phi.of.b, theta.of.b, variance=1)
{
  p <- length(phi.of.b)
  q <- length(theta.of.b)
  
  omega <- seq(from=0, to=pi, by=.001)
  
  phi.of.e.minus.i.omega <- 1
  phi.of.e.i.omega <- 1
  
  if(p>1)
  {   for(i in 2:p)
  {
    phi.of.e.minus.i.omega <-  phi.of.e.minus.i.omega + phi.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
    phi.of.e.i.omega <-  phi.of.e.i.omega + phi.of.b[i]*exp(complex(imaginary = (i-1))*omega)
  }
  }
  
  theta.of.e.minus.i.omega <- 1
  theta.of.e.i.omega <- 1
  
  if(q>1)
  {
    for(i in 2:q)
    {
      theta.of.e.minus.i.omega <-  theta.of.e.minus.i.omega + theta.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
      theta.of.e.i.omega <-  theta.of.e.i.omega + theta.of.b[i]*exp(complex(imaginary = (i-1))*omega)
    }
  }
  
  my.spectrum <- (variance/(2*pi))*Re(theta.of.e.minus.i.omega*theta.of.e.i.omega/(phi.of.e.minus.i.omega*phi.of.e.i.omega))
  
  plot(omega, 10*log10(my.spectrum), ylab="spectrum (in decibels)", type="l")   
}


#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################

##  Get the data.  Save on your computer or...    Change the path if necessary.
HomePrice <- read.delim("HomePrice.txt")
head(HomePrice)
HomePrice
#n=325=12*27+1-------27 years--------1987/1~2014/1
#index takes 2000/1 to be a "baseline" 
#and on that date the index is given a numerical value of 1, 
#and then compares the value of homes at other times

#Most people believe that there was a housing bubble in the 2000's
#where prices increased too rapidly until prices became unsustainable at which point the bubble burst and prices fell
ts.HomePrice<-ts(HomePrice$index, start=c(1987,1), end=c(2014,1), frequency =12)
index<-ts.HomePrice[1:325]



## 1. raw data plot
plot(index, xlab="observations", ylab="index")
#house bubble begins 130-150, Oct-97-Jun-99
#take a liitle closer
plot(1:120, index[1:120], ylim=c(0.6,2.3), xlim=c(0,350), xlab = "time", ylab="index",type="b")
lines(121:156, index[121:156], ylim=c(0.6,2.3), xlim=c(0,350), type="b", col=2)
lines(157, index[157], ylim=c(0.6,2.3), xlim=c(0,350), type="b", col=4, pch=19)
lines(158:234, index[158:234], ylim=c(0.6,2.3), xlim=c(0,350), type="b", col=2)
lines(235:263, index[235:263], ylim=c(0.6,2.3), xlim=c(0,350), type="b", col=3)
lines(264:325, index[264:325], ylim=c(0.6,2.3), xlim=c(0,350), type="b", col=6)
## cutoff at 234=peak
##    Hmmmm   stationary...  need to difference
##  Should always plot the data first!!!
##  Should have known from the acf

my.log.data<-log(index)
plot(my.log.data, type="b")

# 2. differencing
my.df1.data <-diff(index, lag = 1, differences = 1)
plot(my.df1.data, type="b", xlab="time")

my.df2.data <-diff(index, lag = 1, differences = 2)
plot(my.df2.data, type="b", xlab="time")

my.df3.data <-diff(index, lag = 1, differences = 3)
plot(my.df3.data, type="b", xlab="time")

my.df4.data <-diff(index, lag = 1, differences = 4)
plot(my.df4.data, type="b", xlab="time")



## split data into two parts and fit separately
#first part
data.1 <- my.df2.data[1:232]
plot(data.1, type="b")
par(mfrow=c(2,1))
acf(data.1, type = c("correlation"), ylab="Sample ACF")
acf(data.1, type = c("partial"), ylab="Sample PACF")
par(mfrow=c(1,1)) 
#like ARMA(11,1)
my.df2.data <-diff(index, lag = 1, differences = 2)
arima(x=my.df2.data[1:232] , order = c(12, 0, 0), include.mean=FALSE)

##  look only at first 90% 
par(mfrow=c(2,1))
acf(my.df3.data[1:210], type = c("correlation"), ylab="90% Sample ACF")
acf(my.df3.data[1:210], type= c("partial"), ylab="90% Sample PACF")
par(mfrow=c(1,1))
#
my.fit.1 <- arima(x=my.df3.data[1:210] , order = c(11, 0, 0), include.mean=FALSE)
my.real.preds.1 <- predict(my.fit.1, n.ahead = 22, se.fit = TRUE)
my.real.preds.1

range(my.df3.data[1:210])

plot(my.df3.data[1:210], ylim=c(-0.01,0.012), xlim=c(0,236), type="b")

lines(my.df3.data[211:232], my.real.preds.1$pred, type="b", col=2)
lines(my.df3.data[211:232], my.real.preds.1$pred + 2*my.real.preds.1$se, type="l", col=3)
lines(my.df3.data[211:232], my.real.preds.1$pred - 2*my.real.preds.1$se, type="l", col=3)
#take a liitle closer
plot(190:210, my.df3.data[190:210], ylim=c(-0.01,0.012), xlim=c(190,232), type="b")

lines(211:232, my.real.preds.1$pred, type="b", col=2)
lines(211:232, my.real.preds.1$pred + 2*my.real.preds.1$se, type="l", col=3)
lines(211:232, my.real.preds.1$pred - 2*my.real.preds.1$se, type="l", col=3)

cbind(my.real.preds.1$pred, my.df3.data[211:232])

points(211:232, my.df3.data[211:232], col=4)

## periodogram
fit.1<-arima(index[1:234], order=c(12,2,0))
fit.1
AIC(fit.1)
fit.sample.1<-arima.sim(model=list(ar=c(-0.7160,-0.4470,-0.5338,-0.6973,-0.6917,-0.6776,-0.4867,-0.5026,-0.5444,-0.4849,-0.2856)), n=234)
spec.pgram(fit.sample.1, taper=.1)

#Plot the TRUE SPECTRUM for this process (measured in decibels)
my.spectrum(phi.of.b=c(1,0.7160,0.4470,0.5338,0.6973,0.6917,0.6776,0.4867,0.5026,0.5444,0.4849,0.2856), theta.of.b=c(1), variance=1)
#Does the periodogram match up with the theoretical spectrum?
par(mfrow=c(2,1))
spec.pgram(fit.sample.1)  ##  tapr = .1 by default
my.spectrum(phi.of.b=c(1,0.7160,0.4470,0.5338,0.6973,0.6917,0.6776,0.4867,0.5026,0.5444,0.4849,0.2856), theta.of.b=c(1), variance=1)
par(mfrow=c(1,1))

## 6. residuals diagnostics
#try ARIMA(3,1,0)
tsdiag(my.fit.1)   ## can't do this if fit was using "ar"
qqnorm(my.fit.1$residuals)

plot(my.fit.1$residuals)
#Close up of a portion of the data, to see it better.
plot(my.fit.1$residuals[1:100], type="b")






## cut point 234
data.2 <- my.df2.data[233:323]
plot(data.2, type="b")
par(mfrow=c(2,1))
acf(data.2, type = c("correlation"), ylab="Sample ACF")
acf(data.2, type = c("partial"), ylab="Sample PACF")
par(mfrow=c(1,1)) 
#like AR(11)
my.df2.data <-diff(index, lag = 1, differences = 2)
arima(x=my.df2.data[233:323] , order = c(11, 0, 1), include.mean=FALSE)


##  look only at first 90% 
par(mfrow=c(2,1))
acf(my.df2.data[233:314], type = c("correlation"), ylab="90% Sample ACF")
acf(my.df2.data[233:314], type= c("partial"), ylab="90% Sample PACF")
par(mfrow=c(1,1))
#
my.fit.2 <- arima(x=my.df2.data[233:314] , order = c(11, 0, 1), include.mean=FALSE)
my.real.preds.2 <- predict(my.fit.2, n.ahead = 9, se.fit = TRUE)
my.real.preds.2

range(my.df2.data[233:314])

plot(my.df2.data[233:314], ylim=c(-0.0153,0.0234), xlim=c(0,95), type="b")

lines(my.df2.data[315:323], my.real.preds.2$pred, type="b", col=2)
lines(my.df2.data[315:323], my.real.preds.2$pred + 2*my.real.preds.2$se, type="l", col=3)
lines(my.df2.data[315:323], my.real.preds.2$pred - 2*my.real.preds.2$se, type="l", col=3)
#take a liitle closer
plot(300:314, my.df2.data[300:314], ylim=c(-0.0229,0.0234), xlim=c(300,323), type="b")

lines(315:323, my.real.preds.2$pred, type="b", col=2)
lines(315:323, my.real.preds.2$pred + 2*my.real.preds.2$se, type="l", col=3)
lines(315:323, my.real.preds.2$pred - 2*my.real.preds.2$se, type="l", col=3)

cbind(my.real.preds.2$pred, my.df2.data[315:323])

points(315:323, my.df2.data[315:323], col=4)



## periodogram
fit.2<-arima(index[235:325], order=c(11,2,1))
fit.2
fit.sample.2<-arima.sim(model=list(ar=c(0.6320,-0.1647,-0.1301,0.0539,-0.0799,-0.0835,-0.0388,0.0832,-0.2559,0.1012,0.3375),ma=c(-0.3343)), n=91)
spec.pgram(fit.sample.2, taper=.1)

#Plot the TRUE SPECTRUM for this process (measured in decibels)
my.spectrum(phi.of.b=c(1,-0.6320,0.1647,0.1301,-0.0539,0.0799,0.0835,0.0388,-0.0832,0.2559,-0.1012,-0.3375), theta.of.b=c(1-0.3343), variance=1)
#Does the periodogram match up with the theoretical spectrum?
par(mfrow=c(2,1))
spec.pgram(fit.sample.2)  ##  tapr = .1 by default
my.spectrum(phi.of.b=c(1,-0.6320,0.1647,0.1301,-0.0539,0.0799,0.0835,0.0388,-0.0832,0.2559,-0.1012,-0.3375), theta.of.b=c(1-0.3343), variance=1)
par(mfrow=c(1,1))

## 6. residuals diagnostics
#try ARIMA(3,1,0)
tsdiag(my.fit.2)   ## can't do this if fit was using "ar"
qqnorm(my.fit.2$residuals)

plot(my.fit.2$residuals)
#Close up of a portion of the data, to see it better.
plot(my.fit.2$residuals[1:81], type="b")




##########################################
## 7. Make predictions....
##########################################

#choose (3,1,0)
my.final.fit <- arima(index[235:325], order=c(11,2,1))
my.preds <- predict(my.final.fit, n.ahead = 12, se.fit = TRUE)
names(my.preds)
my.preds

my.lower.limits  <- my.preds$pred - 2*my.preds$se
my.upper.limits <- my.preds$pred + 2*my.preds$se

cbind(my.lower.limits, my.preds$pred, my.upper.limits)


plot(300:325, index[300:325], ylim=c(1.42,2.3), xlim=c(300,340), type="b")

lines(326:337, my.preds$pred, type="b", col=2)
lines(326:337, my.upper.limits, type="l", col=3)
lines(326:337, my.lower.limits, type="l", col=3)






my.final.fit <- arima(index[235:325], order=c(11,2,1))
my.preds <- predict(my.final.fit, n.ahead = 100, se.fit = TRUE)
names(my.preds)
my.preds

my.lower.limits  <- my.preds$pred - 2*my.preds$se
my.upper.limits <- my.preds$pred + 2*my.preds$se

cbind(my.lower.limits, my.preds$pred, my.upper.limits)


plot(1:500, index[1:500], ylim=c(0,5), xlim=c(1,500), type="b")

lines(326:425, my.preds$pred, type="b", col=2)
lines(326:425, my.upper.limits, type="l", col=3)
lines(326:425, my.lower.limits, type="l", col=3)
